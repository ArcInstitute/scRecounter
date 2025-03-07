#!/usr/bin/env python3
# import
## batteries
import os
import logging
import argparse 
from pathlib import Path
from itertools import chain
from typing import List, Set, Tuple, Optional
## 3rd party
import pandas as pd
from pypika import Query, Table, Criterion
## package
from db_utils import db_connect

# format logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

# classes
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass

# functions
def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
    desc = 'Find scRNA-seq count matrix files.'
    epi = """DESCRIPTION:
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=CustomFormatter)
    parser.add_argument(
        'base_dir',  type=str, help='Base directory to search for input data files'
    )
    parser.add_argument(
        '--feature-type', default='GeneFull_Ex50pAS', 
        choices=['Gene', 'GeneFull', 'GeneFull_Ex50pAS', 'GeneFull_ExonOverIntron', 'Velocyto'], 
        help='Feature type to process'
    )
    parser.add_argument(
        '--raw', action='store_true', default=False,
        help='Use raw count matrix files instead of filtered'
    )
    parser.add_argument(
        '--max-datasets', type=int, default=0,
        help='Maximum number of datasets to process'
    )
    parser.add_argument( 
        '--multi-mapper', default=['Unique'], nargs='+',
        choices=['Unique', 'EM', 'Uniform'],
        help='Multi-mapper strategy to use' 
    )
    parser.add_argument(
        '--organisms', type=str, default=None,
        help="Comma-separated list of organisms to process; if none, process all"
    )
    parser.add_argument(
        '--redo-processed', action='store_true', default=False,
        help="Do not skip already procssed SRX IDs, but instead include them in the output"
    )
    return parser.parse_args()

def load_srx_metadata(organisms: str) -> Set[str]:
    """
    Load metadata from scBasecamp database.
    Args:
        organisms: Comma-separated list of organisms to process; if none, process all
    """
    logging.info("Obtaining srx metadata...")

    # get metadata from scRecounter postgresql database
    srx_metadata = Table("srx_metadata")
    stmt = (
        Query
        .from_(srx_metadata)
        .select(
            srx_metadata.srx_accession,
            srx_metadata.organism,
        ).where(
            Criterion.all([
                ~srx_metadata.is_illumina.isnull(),
                ~srx_metadata.organism.isnull(),
                ~srx_metadata.organism.isin(['NA', 'None', 'NaN', 'other', 'metagenome']),
            ])
        )
    )
    
    # load metadata
    with db_connect() as conn:
        metadata = pd.read_sql(str(stmt), conn)

    # filter by organism
    if organisms:
        organisms = organisms.split(',')
        metadata = metadata[metadata['organism'].isin(organisms)]
        
    # return set of SRX accessions
    srx = set(metadata['srx_accession'].tolist())
    logging.info(f"Found {len(srx)} SRX accessions with metadata.")
    return srx

def load_scbasecamp_metadata(feature_type: str) -> Set[str]:
    """
    Load metadata from scBasecamp database.
    """
    logging.info("Obtaining scbasecamp metadata...")

    # get metadata from scRecounter postgresql database
    scbc_metadata = Table("scbasecamp_metadata")
    stmt = (
        Query
        .from_(scbc_metadata)
        .select(
            scbc_metadata.srx_accession,
        ).where(
            scbc_metadata.feature_type == feature_type
        )
    )
    with db_connect() as conn:
        metadata = pd.read_sql(str(stmt), conn)
    srx = set(metadata['srx_accession'].tolist())
    logging.info(f"Found {len(srx)} SRX accessions already processed.")
    return srx

def find_matrix_files(
        base_dir: str, 
        feature_type: str, 
        has_srx_metadata: Set[str],
        processed_srx: Set[str],
        multi_mapper: List[str],
        raw: bool=False, 
        max_datasets: Optional[int]=0
    ) -> List[tuple]:
    """
    Recursively find *.mtx.gz files and extract SRX/ERX IDs.
    Args:
        base_dir: Base directory to search
        feature_type: 'Gene' or 'GeneFull'
        has_srx_metadata: Set of SRX IDs with metadata; records skipped if no metadata
        processed_srx: Set of existing SRX IDs
        multi_mapper: 'EM', 'uniform', and/or 'Unique'
        raw: Use raw count matrix files instead of filtered
        max_datasets: Maximum number of datasets to process
    Returns:
        List of tuples (matrix_path, srx_id)
    """
    logging.info(f"Searching for new data files in {base_dir}...")
    base_path = Path(base_dir)
    subdir = 'raw' if raw else 'filtered'
    results = []
    stats = {
        'found': 0, 
        'has_metadata': 0, 
        'no_metadata': 0, 
        'already_processed': 0, 
        'permissions': 0,
        'mtx_file_missing': 0, 
        'novel': 0
    }
    
    # Determine which matrix file to look for based on multi_mapper
    matrix_filenames = []
    if "Velocyto" in feature_type: 
        matrix_filenames += ["ambiguous.mtx.gz", "spliced.mtx.gz", "unspliced.mtx.gz"]
    else:
        if 'Unique' in multi_mapper:
            matrix_filenames = ['matrix.mtx.gz']
        if 'EM' in multi_mapper:
            matrix_filenames += ['UniqueAndMult-EM.mtx.gz']
        if 'Uniform' in multi_mapper:
            matrix_filenames += ['UniqueAndMult-Uniform.mtx.gz']

    # set max datasets
    max_datasets *= len(matrix_filenames)

    # Walk through directory structure
    num_dirs = 0
    for srx_dir in chain(base_path.glob('**/SRX*'), base_path.glob('**/ERX*')):
        # skip files
        if not srx_dir.is_dir():
            continue
        else:
            stats['found'] += 1

        # status
        num_dirs += 1
        if num_dirs % 1000 == 0:
            logging.info(f"  Searched {num_dirs} SRX directories so far...")

        # Check if SRX directory exists in database
        if srx_dir.name in processed_srx:
            stats['already_processed'] += 1
            continue

        # Check if SRX directory exists in srx_metadata
        if srx_dir.name in has_srx_metadata:
            stats['has_metadata'] += 1
        else:
            stats['no_metadata'] += 1
            continue

        # Find target matrix file in SRX directory
        for matrix_filename in matrix_filenames:
            # EM and Uniform are in the raw directory
            if matrix_filename in ['UniqueAndMult-EM.mtx.gz', 'UniqueAndMult-Uniform.mtx.gz']:
                subdir_tmp = 'raw'
            else:
                subdir_tmp = subdir
            # iterate through the glob results
            for mtx_file in srx_dir.glob(f'**/{matrix_filename}'):
                # check for `feature_type/subdir` in file path
                hit = None
                for i,x in enumerate(mtx_file.parts):
                    try:
                        if feature_type == x and mtx_file.parts[i+1] == subdir_tmp:
                            hit = True
                            break
                    except IndexError:
                        continue
                # if target file found, check if it exists, and add to results
                if hit:
                    features_file = mtx_file.parent / "features.tsv.gz"
                    barcodes_file = mtx_file.parent / "barcodes.tsv.gz"
                else:
                    continue
                try:
                    if not mtx_file.exists() or not features_file.exists() or not barcodes_file.exists():
                        stats['mtx_file_missing'] += 1
                    else:
                        stats['novel'] += 1
                        results.append([srx_dir.name, mtx_file, features_file, barcodes_file])  
                except PermissionError:
                    logging.warning(f"Permission denied for {mtx_file}. Skipping.")
                    stats['permissions'] += 1
                #break
        
        # Check max datasets
        if max_datasets > 0 and len(results) >= max_datasets:
            logging.info(f"  Found --max-datasets datasets. Stopping search.")
            break

    # Status
    logging.info(f"  {stats['found']} total SRX directories found (total).")
    logging.info(f"  {stats['has_metadata']} has srx metadata (kept).")
    logging.info(f"  {stats['no_metadata']} lacks srx metadata (skipped).")
    logging.info(f"  {stats['already_processed']} existing SRX directories found (skipped).")
    logging.info(f"  {stats['mtx_file_missing']} missing matrix files (skipped).")
    logging.info(f"  {stats['permissions']} directories with permission errors (skipped).")
    logging.info(f"  {stats['novel']} novel SRX directories found (final).")
    return results

def main():
    """Main function to run the TileDB loader workflow."""
    args = parse_arguments()

    # Load metadata
    has_srx_metadata = load_srx_metadata(args.organisms)
    
    # Load metadata
    if args.redo_processed:
        logging.warning("Redoing processed SRX IDs.")
        processed_srx = set()
    else:
        processed_srx = load_scbasecamp_metadata(args.feature_type)

    # Find all matrix files and their corresponding SRX IDs
    matrix_files = find_matrix_files(
        args.base_dir, args.feature_type, 
        has_srx_metadata = has_srx_metadata, 
        processed_srx = processed_srx,
        multi_mapper=args.multi_mapper,
        raw=args.raw, 
        max_datasets=args.max_datasets,
    )

    # convert to dataframe
    df = pd.DataFrame(
        matrix_files, columns=['srx', 'matrix_path', 'features_path', 'barcodes_path']
    ).sort_values(['srx'])

    # sort by srx and matrix_path and drop duplicate of the same srx+path
    df = df.sort_values(by=['srx', 'matrix_path'])
    df["basename"] = df["matrix_path"].apply(lambda x: x.name)
    df = df.drop_duplicates(subset=['srx', 'basename'], keep='last').drop(columns=['basename'])

    # if feature_type is Velocyto, check for 3 per SRX and filter incomplete records
    if args.feature_type == "Velocyto":
        # identify SRX with complete records (3 files)
        complete_srx = df.groupby('srx').filter(lambda x: len(x) == 3)
        # filter to keep only complete records
        df = complete_srx.copy()
        num_filtered = len(set(df['srx'])) - len(set(complete_srx['srx']))
        if num_filtered > 0:
            logging.warning(f"Filtered {num_filtered} SRX records that did not have all 3 Velocyto matrix files")

    # write as csv
    df.to_csv('mtx_files.csv', index=False)
    logging.info(f"File written: mtx_files.csv")

if __name__ == "__main__":
    #from dotenv import load_dotenv
    #load_dotenv(override=True)
    main()
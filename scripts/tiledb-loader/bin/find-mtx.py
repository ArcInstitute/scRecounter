#!/usr/bin/env python3
# import
## batteries
import os
import logging
import argparse 
from pathlib import Path
from itertools import chain, repeat
from typing import List, Set, Tuple, Optional
## 3rd party
import pandas as pd
import tiledbsoma
import tiledbsoma.io
from pypika import Query, Table, Criterion
## package
from db_utils import db_connect

# format logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
logging.getLogger("tiledbsoma").setLevel(logging.WARNING)
logging.getLogger("tiledbsoma.io").setLevel(logging.WARNING)
logging.getLogger("tiledb").setLevel(logging.WARNING) 

# classes
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass

# functions
def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
    desc = 'Find scRNA-seq count matrix files for TileDB loader.'
    epi = """DESCRIPTION:
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=CustomFormatter)
    parser.add_argument(
        'base_dir',  type=str, help='Base directory to search for input data files'
    )
    parser.add_argument(
        '--feature-type', default='GeneFull_Ex50pAS', 
        choices=['Gene', 'GeneFull', 'GeneFull_Ex50pAS', 'GeneFull_ExonOverIntron', 'Velocyto', None], 
        help='Feature type to process'
    )
    parser.add_argument(
        '--raw', action='store_true', default=False,
        help='Use raw count matrix files instead of filtered'
    )
    parser.add_argument(
        '--db-uri', type=str, default="tiledb_exp", 
        help='URI of existing TileDB database, or it will be created if it does not exist'
    )
    parser.add_argument(
        '--batch-size', type=int, default=8, help='batch size for downstream processing'
    )
    parser.add_argument(
        '--max-datasets', type=int, default=None,
        help='Maximum number of datasets to process'
    )
    parser.add_argument(
        '--organisms', type=str, default=None,
        help="Comma-separated list of organisms to process; if none, process all"
    )
    parser.add_argument( 
        '--multi-mapper', default='None', choices=['None', 'EM', 'uniform'],
        help='Multi-mapper strategy to use' 
    )
    parser.add_argument(
        '--redo-processed', action='store_true', default=False,
        help="Do not skip already procssed SRX IDs, but instead include them in the output"
    )
    return parser.parse_args()

def get_tiledb_srx_ids(db_uri: str) -> Set[str]:
    """
    Read metadata from existing database and return set of SRX IDs.
    Args:
        db_uri: URI of the TileDB database
    Returns:
        Set of SRX IDs already in the database
    """
    logging.info(f"Checking for existing SRX accessions in {db_uri}...")

    srx = set()
    if not os.path.exists(db_uri):
        logging.info("Database does not exist yet. No SRX/ERX accessions to obtain.")
    else:
        with tiledbsoma.open(db_uri) as exp:
            try:
                metadata = (exp.obs.read(column_names=["SRX_accession"])
                    .concat()
                    .group_by(["SRX_accession"])
                    .aggregate([
                        ([], 'count_all'),
                    ])
                    .to_pandas())
                srx = set(metadata["SRX_accession"].unique())
            except tiledbsoma._exception.DoesNotExistError:
                metadata = (exp.obs.read(column_names=["SRX_accession"])
                    .concat()
                    .to_pandas())
                srx = set(metadata["SRX_accession"].unique())
    # status
    logging.info(f"  Found {len(srx)} SRX/ERX accessions in the tiledb database.")
    return srx

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
    logging.info(f"  Found {len(srx)} SRX/ERX accessions with metadata.")
    return srx

# def load_scbasecamp_metadata(feature_type: str) -> Set[str]:
#     """
#     Load metadata from scBasecamp database.
#     """
#     logging.info("Obtaining scbasecamp metadata...")

#     # get metadata from scRecounter postgresql database
#     scbc_metadata = Table("scbasecamp_metadata")
#     stmt = (
#         Query
#         .from_(scbc_metadata)
#         .select(
#             scbc_metadata.srx_accession,
#         ).where(
#             scbc_metadata.feature_type == feature_type
#         )
#     )
#     with db_connect() as conn:
#         metadata = pd.read_sql(str(stmt), conn)
#     return set(metadata['srx_accession'].tolist())

def find_matrix_files(
        base_dir: str, 
        feature_type: str, 
        has_srx_metadata: Set[str],
        processed_srx: Set[str],
        multi_mapper: str='None',
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
        multi_mapper: 'EM', 'uniform', or 'None'
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

    # account for all 3 matrix files if feature_type is Velocyto
    if feature_type == "Velocyto":
        if max_datasets > 0:
            max_datasets = max_datasets * 4
    
    # Determine which matrix file to look for based on multi_mapper
    if multi_mapper == 'None':
        matrix_filename = ['matrix.mtx.gz']
        if "Velocyto" in feature_type:
            matrix_filename += ["ambiguous.mtx.gz", "spliced.mtx.gz", "unspliced.mtx.gz"]
    elif multi_mapper == 'EM':
        matrix_filename = ['UniqueAndMult-EM.mtx.gz']
    elif multi_mapper == 'uniform':
        matrix_filename = ['UniqueAndMult-Uniform.mtx.gz']
    else:
        raise ValueError(f"Invalid multi-mapper strategy: {multi_mapper}")

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
        mtx_files = []
        for f in matrix_filename:
            mtx_files.extend(srx_dir.glob(f'**/{f}'))
        for mtx_file in mtx_files:
            hit = None
            # check for `feature_type/subdir` in file path
            for i,x in enumerate(mtx_file.parts):
                try:
                    if feature_type == x and mtx_file.parts[i+1] == subdir:
                        hit = True
                        break
                except IndexError:
                    continue
            # if target file found, check if it exists, and add to results
            if hit:
                features_file = mtx_file.parent / "features.tsv.gz"
                barcodes_file = mtx_file.parent / "barcodes.tsv.gz"
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

# def make_batch(num_repeats: int, total_numbers: int, feature_type: str) -> List[int]:
#     """
#     Bin numbers into batches of num_repeats.
#     Args:
#         num_repeats: Number of repeats per unique number
#         total_numbers: Total number of unique numbers
#     Returns:
#         List of batch numbers
#     """
#     batch_counts = []
#     unique_count = int(round(total_numbers / num_repeats + 0.5))
#     for i in range(1, unique_count + 1):
#         batch_counts.extend(repeat(i, num_repeats))
#     return batch_counts[:total_numbers]

def make_batch(df: pd.DataFrame, batch_size: int) -> pd.DataFrame:
    unique_srx = sorted(df['srx'].drop_duplicates().tolist())
    batch_mapping = {}
    batch_num = 1
    count = 0
    for s in unique_srx:
        if count >= batch_size:
            batch_num += 1
            count = 0
        batch_mapping[s] = batch_num
        count += 1
    df["batch"] = df["srx"].map(batch_mapping)
    return df

def main():
    """Main function to run the TileDB loader workflow."""
    args = parse_arguments()

    # Load scRecounter SQL db records
    has_srx_metadata = load_srx_metadata(args.organisms)
    
    # Load tiledb records
    if args.redo_processed:
        processed_srx = set()
    else:
        processed_srx = get_tiledb_srx_ids(args.db_uri)

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

    # assign batches ensuring all records for the same SRX are in the same batch
    df = make_batch(df, args.batch_size)

    # write as csv
    df.to_csv('mtx_files.csv', index=False)
    logging.info(f"File written: mtx_files.csv")

if __name__ == "__main__":
    main()
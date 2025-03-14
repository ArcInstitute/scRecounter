#!/usr/bin/env python3
# import
## batteries
import os
import gzip
import logging
import argparse
from typing import List, Set, Tuple, Dict, Optional, Union, TextIO
## 3rd party
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from scipy import sparse
from pypika import Query, Table
## package
from db_utils import db_connect, db_upsert

# format logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
logging.getLogger("psycopg2").setLevel(logging.CRITICAL)
logging.getLogger("google.auth.transport.requests").setLevel(logging.CRITICAL)
logging.getLogger("urllib3").setLevel(logging.CRITICAL)
logging.getLogger("google.auth").setLevel(logging.CRITICAL)

# classes
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass

# functions
def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
    desc = 'Convert mtx files to h5ad.'
    epi = """DESCRIPTION:
    Convert mtx files to h5ad in parallel.
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=CustomFormatter)
    parser.add_argument(
        '--srx', type=str, help="SRX accessions", required=True
    )
    parser.add_argument(
        '--matrix-paths', type=str, nargs="+", help="Path to >=1 *.mtx.gz file", required=True
    )
    parser.add_argument(
        '--feature-paths', type=str, nargs="+", help="Path to >=1 feature file", required=True
    )
    parser.add_argument(
        '--barcode-paths', type=str, nargs="+", help="Path to >=1 barcode file", required=True
    )
    parser.add_argument(
        '--matrix-types', type=str, nargs="+", help="Matrix types", required=True
    )
    parser.add_argument(
        '--publish-path', type=str, help="Publishing path", required=True
    )
    parser.add_argument(
        '--tissue-categories', type=str, help="Tissue category csv file", required=True
    )
    parser.add_argument(
        '--missing-metadata', type=str, default="error", 
        choices=["error", "skip", "allow"],
        help="How do handle missing metadata?"
    )
    parser.add_argument(
        '--feature-type', default='GeneFull_Ex50pAS', 
        choices=['Gene', 'GeneFull', 'GeneFull_Ex50pAS', 'GeneFull_ExonOverIntron', 'Velocyto'], 
        help='Feature type to process'
    )
    parser.add_argument(
        '--update-database', action="store_true", default=False, 
        help="Update the database?"
    )
    return parser.parse_args()

def get_metadata(srx_id: str, missing_metadata: str="error") -> Optional[pd.DataFrame]:
    """
    Get metadata for an SRX accession.
    """
    logging.info("Obtaining srx_metadata...")

    # get metadata from scRecounter postgresql database
    srx_metadata = Table("srx_metadata")
    stmt = (
        Query
        .from_(srx_metadata)
        .select(
            srx_metadata.entrez_id,
            srx_metadata.srx_accession,
            srx_metadata.lib_prep, 
            srx_metadata.tech_10x,
            srx_metadata.cell_prep,
            srx_metadata.organism,
            srx_metadata.tissue,
            srx_metadata.disease,
            srx_metadata.perturbation,
            srx_metadata.cell_line, 
            srx_metadata.czi_collection_id,
            srx_metadata.czi_collection_name,
        )
        .where(srx_metadata.srx_accession == srx_id) 
    )
    metadata = None
    with db_connect() as conn:
        metadata = pd.read_sql(str(stmt), conn)

    ## if metadata is not found, return None
    if metadata is None or metadata.shape[0] == 0:
        if missing_metadata == "allow":
            logging.warning(
                f"    Metadata not found for SRX accession {srx_id}, but `--missing-metadata allow` used"
            )
            pass
        elif missing_metadata == "skip":
            logging.warning(
                f"    Metadata not found for SRX accession {srx_id}, but `--missing-metadata skip` used"
            )
            return None
        elif missing_metadata == "error":
            raise ValueError(f"    Metadata not found for SRX accession {srx_id}")
        else:
            raise ValueError(f"    Invalid value for `--missing-metadata`")
    if metadata.shape[0] > 1:
        raise ValueError(f"Multiple metadata entries found for SRX accession {srx_id}")
    elif metadata.shape[0] == 1:
        # lib_prep should be "10x_Genomics"
        if metadata["lib_prep"].values[0] != "10x_Genomics":
            metadata["lib_prep"] = "10x_Genomics"
            metadata["tech_10x"] = "other"
    return metadata

def add_tissue_category(metadata: pd.DataFrame, tissue_categories_path: str) -> None:
    """
    Add tissue category to metadata.
    """
    # if tissue_categories_path exists, update tissue to category
    if os.path.exists(tissue_categories_path):
        logging.info(f"Loading tissue categories...")
        df = pd.read_csv(tissue_categories_path)  
        tissue_category = df[df["tissue"] == metadata["tissue"].values[0]]
        if tissue_category.shape[0] > 0:
            metadata["tissue"] = tissue_category["category"].values[0]
    return metadata

def rename_files(files: List[str], matrix_types: List[str], prefix: str) -> List[str]:
    """
    Rename files based on matrix_types
    """
    new_files = {}
    for file, matrix_type in zip(files, matrix_types):
        if matrix_type == "matrix":
            new_file = f"{prefix}.tsv.gz" 
        else:    
            new_file = f"{prefix}_{matrix_type}.tsv.gz" 
        if os.path.exists(file):
            os.rename(file, new_file)
        new_files[matrix_type] = new_file
    return new_files

def build_velocyto_anndata(matrix_paths: Dict[str,str], feature_paths: Dict[str,str], barcode_paths: Dict[str,str]) -> sc.AnnData:
    """
    Generate an anndata object from the STAR aligner output folder
    """
    # check exists
    for _,matrix_path in matrix_paths.items():
        if not os.path.exists(matrix_path):
            raise FileNotFoundError(f"{matrix_path} not found")

    # Transpose counts matrix to have Cells as rows and Genes as cols as expected by AnnData objects
    ## Using spliced.mtx.gz as the reference matrix
    X = sc.read_mtx(matrix_paths['spliced']).X.transpose()

    # Load the 3 matrices containing Spliced, Unspliced and Ambigous reads
    mtxU = np.loadtxt(matrix_paths['unspliced'], skiprows=3, delimiter=' ')
    mtxS = np.loadtxt(matrix_paths['spliced'], skiprows=3, delimiter=' ')
    mtxA = np.loadtxt(matrix_paths['ambiguous'], skiprows=3, delimiter=' ')

    # Extract sparse matrix shape informations from the third row
    shapeU = np.loadtxt(matrix_paths['unspliced'], skiprows=2, max_rows = 1 ,delimiter=' ')[0:2].astype(int)
    shapeS = np.loadtxt(matrix_paths['spliced'], skiprows=2, max_rows = 1 ,delimiter=' ')[0:2].astype(int)
    shapeA = np.loadtxt(matrix_paths['ambiguous'], skiprows=2, max_rows = 1 ,delimiter=' ')[0:2].astype(int)

    # Read the sparse matrix with csr_matrix((data, (row_ind, col_ind)), shape=(M, N))
    # Subract -1 to rows and cols index because csr_matrix expects a 0 based index
    # Traspose counts matrix to have Cells as rows and Genes as cols as expected by AnnData objects
    spliced = sparse.csr_matrix((mtxS[:,2], (mtxS[:,0]-1, mtxS[:,1]-1)), shape = shapeS).transpose()
    unspliced = sparse.csr_matrix((mtxU[:,2], (mtxU[:,0]-1, mtxU[:,1]-1)), shape = shapeU).transpose()
    ambiguous = sparse.csr_matrix((mtxA[:,2], (mtxA[:,0]-1, mtxA[:,1]-1)), shape = shapeA).transpose()

    # Load Genes and Cells identifiers
    obs = pd.read_csv(barcode_paths['spliced'], header = None, index_col = 0)

    # Remove index column name to make it compliant with the anndata format
    obs.index.name = None
    var = pd.read_csv(feature_paths['spliced'], sep='\t', names = ('gene_ids', 'feature_types'), index_col = 1)
  
    # Build AnnData object to be used with ScanPy and ScVelo
    adata = anndata.AnnData(
        X = X, obs = obs, var = var,
        layers = {'spliced': spliced, 'unspliced': unspliced, 'ambiguous': ambiguous}
    )
    adata.var_names_make_unique()

    # Subset Cells based on STAR filtering
    selected_barcodes = pd.read_csv(barcode_paths['spliced'], header = None)
    return adata[selected_barcodes[0]]

def open_file(filename: str, mode: str = 'r') -> Union[TextIO, gzip.GzipFile]:
    """Open a file, handling gzip if the filename ends with .gz"""
    if filename.endswith('.gz'):
        return gzip.open(filename, mode + 't')  # Text mode for gzip
    else:
        return open(filename, mode)

def match_barcodes(
        barcode_file: str, matrix_file: str, target_barcode_file: str, 
        out_cb_file: str, out_mat_file: str
    ) -> None:
    """
    Filter matrix to keep only cells that match target barcodes
    Args:
        barcode_file: Path to input barcode file
        matrix_file: Path to input matrix file
        target_barcode_file: Path to target barcodes file
        out_cb_file: Path to output barcodes file
        out_mat_file: Path to output matrix file
    """
    
    logging.info(f"Starting barcode matching")
    logging.info(f"  Input barcode file: {barcode_file}")
    logging.info(f"  Input matrix file: {matrix_file}")
    logging.info(f"  Target barcodes file: {target_barcode_file}")
    
    # Read original barcodes and create mapping
    barcode_to_index = {}
    with open_file(barcode_file) as f:
        for i, line in enumerate(f, 1):
            barcode = line.strip()
            barcode_to_index[barcode] = i
    
    logging.info(f"  Found {len(barcode_to_index)} original barcodes")
    
    # Read target barcodes
    target_barcodes = set()
    with open_file(target_barcode_file) as f:
        for line in f:
            target_barcodes.add(line.strip())
    
    logging.info(f"  Found {len(target_barcodes)} target barcodes (used to filter)")
    
    # Find matching barcodes and create new index mapping
    matched_indices = {}  # Maps original index to new index
    matched_barcodes = []  # List of matched barcodes
    
    for barcode, orig_idx in sorted(barcode_to_index.items(), key=lambda x: x[1]):
        if barcode in target_barcodes:
            new_idx = len(matched_indices) + 1
            matched_indices[orig_idx] = new_idx
            matched_barcodes.append(barcode)
    
    logging.info(f"  Found {len(matched_barcodes)} matching barcodes")
    
    # Write matched barcodes
    with open(out_cb_file, 'w') as out_f:
        for barcode in matched_barcodes:
            out_f.write(f"{barcode}\n")
    
    # Process matrix file
    header_lines = []
    filtered_entries = []
    n_features = 0
    
    with open_file(matrix_file) as f:
        # Process header lines (comments)
        line = f.readline().strip()
        while line.startswith('%'):
            header_lines.append(line)
            line = f.readline().strip()
        
        # Process dimensions line
        dimensions = line.split()
        n_features = int(dimensions[0])
        
        # Process data lines
        for i,line in enumerate(f,1):
            cols = line.strip().split()
            gene_idx = cols[0]
            try:
                cell_idx = int(cols[1])
            except (IndexError, ValueError) as e:
                raise ValueError(f"Line {i}: invalid matrix value: \"{line.strip()}\"")
            umi_count = cols[2]
            if cell_idx in matched_indices:
                filtered_entries.append((gene_idx, matched_indices[cell_idx], umi_count))
    
    # Write filtered matrix
    with open(out_mat_file, 'w') as out_f:
        # Write header comments
        for line in header_lines:
            out_f.write(f"{line}\n")
        
        # Write new dimensions
        out_f.write(f"{n_features} {len(matched_barcodes)} {len(filtered_entries)}\n")
        
        # Write filtered entries
        for gene_idx, cell_idx, umi_count in filtered_entries:
            out_f.write(f"{gene_idx} {cell_idx} {umi_count}\n")
    
    logging.info(f"  Filtered barcodes saved to {out_cb_file}")
    logging.info(f"  Filtered matrix saved to {out_mat_file}")
    logging.info(f"  Barcode matching complete: {len(matched_barcodes)} cells, {len(filtered_entries)} matrix entries")

def build_gene_anndata(matrix_paths: Dict[str, str], feature_paths: Dict[str,str], barcode_paths: Dict[str, str]) -> sc.AnnData:
    """
    Generate an anndata object from the STAR aligner output folder
    Args:
        matrix_paths: Path to matrix files, {matrix_type: path}
        feature_paths: Path to feature files, {matrix_type: path}
        barcode_paths: Path to barcode files, {matrix_type: path}
    Returns:
        Anndata object
    """    
    # primary load count matrix
    logging.info("Loading primary count matrix...")
    adata = sc.read_10x_mtx(
        os.path.dirname(matrix_paths['matrix']),
        var_names="gene_ids",
        make_unique=True
    )
    
    # filtering multi-mapper count matrices
    logging.info("Adding multi-mapper count matrices as layers...")
    for matrix_type in ['UniqueAndMult-Uniform', 'UniqueAndMult-EM']:
        logging.info(f"Filtering {matrix_type} matrix...")
        match_barcodes(
            barcode_paths[matrix_type],
            matrix_paths[matrix_type],
            barcode_paths['matrix'],
            f'barcodes_{matrix_type}_filtered.tsv',
            f'{matrix_type}_filtered.mtx'
        )
        adata.layers[matrix_type] = sc.read_mtx(f'{matrix_type}_filtered.mtx').X.transpose()
    return adata


def load_matrix_as_anndata(
        srx_id: str, 
        metadata: pd.DataFrame,
        matrix_paths: Dict[str, str],
        feature_paths: Dict[str, str],
        barcode_paths: Dict[str, str],
        publish_path: str,
        feature_type: str="GeneFull_Ex50pAS",
        update_database: bool=False
    ) -> sc.AnnData:
    """
    Load a matrix.mtx.gz file as an AnnData object.
    Args:
        srx_id: SRX accession
        metadata: Metadata for the SRX accession
        matrix_paths: Path to matrix files
        feature_paths: Path to feature files
        barcode_paths: Path to barcode files
        publish_path: Path to publish directory
        feature_type: Feature type
        update_database: Update the database?
    Returns:
        AnnData object
    """
    # add publish path
    metadata["file_path"] = metadata["organism"].apply(
        lambda org: os.path.join(publish_path, "h5ad", feature_type, str(org).replace(" ", "_"), f"{srx_id}.h5ad")
    )

    # build anndata
    if all(x in matrix_paths.keys() for x in ["spliced", "unspliced", "ambiguous"]):
        logging.info("Building Velocyto anndata...")
        adata = build_velocyto_anndata(matrix_paths, feature_paths, barcode_paths)
    elif 'matrix' in matrix_paths.keys():
        logging.info("Building gene anndata...")
        adata = build_gene_anndata(matrix_paths, feature_paths, barcode_paths)
    else:
        raise ValueError("Invalid matrix_paths")
    
    # calculate total counts
    if sparse.issparse(adata.X):
        adata.obs["gene_count"] = (adata.X > 0).sum(axis=1).A1
        adata.obs["umi_count"] = adata.X.sum(axis=1).A1
    else:
        adata.obs["gene_count"] = (adata.X > 0).sum(axis=1)
        adata.obs["umi_count"] = adata.X.sum(axis=1)

    # add metadata to adata
    adata.obs["SRX_accession"] = srx_id

    # add obs_count to metadata
    metadata["obs_count"] = adata.shape[0]

    ## write to h5ad
    outdir = os.path.join("h5ad", feature_type, metadata["organism"].values[0].replace(" ", "_"))
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, f"{srx_id}.h5ad")
    logging.info(f"Writing to {outfile}...")
    adata.write_h5ad(outfile, compression="gzip")

    # write out obs dataframe as csv
    os.makedirs("metadata", exist_ok=True)
    outfile = os.path.join("metadata", f"{srx_id}.csv.gz")
    adata.obs["cell_barcode"] = adata.obs.index
    adata.obs["organism"] = metadata["organism"].values[0]
    adata.obs.to_csv(outfile, index=False, compression="gzip")

    # add feature type
    metadata["feature_type"] = feature_type

    # upsert metadata to postgresql database
    if update_database:
        logging.info(f"Upserting metadata for SRX accession {srx_id}...")
        with db_connect() as conn:
            db_upsert(metadata, "scbasecamp_metadata_tmp", conn)
    else:
        logging.info(f"Skipping upserting metadata for SRX accession {srx_id}")

def main():
    # parse args
    args = parse_arguments()
    # get metadata
    metadata = get_metadata(args.srx, args.missing_metadata)
    if metadata is None:
        return None
    ## add tissue category to metadata
    metadata = add_tissue_category(metadata, args.tissue_categories)

    # rename feature and barcode files
    args.matrix_paths = dict(zip(args.matrix_types, args.matrix_paths))
    args.feature_paths = rename_files(args.feature_paths, args.matrix_types, prefix="features")
    args.barcode_paths = rename_files(args.barcode_paths, args.matrix_types, prefix="barcodes")

    # Load mtx file
    load_matrix_as_anndata(
        srx_id = args.srx, 
        metadata = metadata,
        matrix_paths = args.matrix_paths, 
        feature_paths = args.feature_paths,
        barcode_paths = args.barcode_paths,
        publish_path = args.publish_path, 
        feature_type = args.feature_type,
        update_database = args.update_database
    )

if __name__ == "__main__":
    #from dotenv import load_dotenv
    #load_dotenv(override=True)
    main()
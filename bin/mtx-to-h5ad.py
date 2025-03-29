#!/usr/bin/env python
# import
## batteries
import os
import sys
import gzip
import json
import argparse
import logging
from glob import glob
from collections import defaultdict
from typing import Optional, List, Tuple, Dict, TextIO, Union
## 3rd party
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from scipy import sparse
from pypika import Query, Table
## package
from db_utils import db_connect, db_upsert


# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
logging.getLogger("psycopg2").setLevel(logging.CRITICAL)
logging.getLogger("google.auth.transport.requests").setLevel(logging.CRITICAL)
logging.getLogger("urllib3").setLevel(logging.CRITICAL)
logging.getLogger("google.auth").setLevel(logging.CRITICAL)


# argparse
FEATURE_TYPES = ["Gene", "GeneFull", "GeneFull_Ex50pAS", "GeneFull_ExonOverIntron", "Velocyto"]
#MTX_TYPES = ["Unique", "EM", "Uniform", "Veloc"

def parse_args():
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                          argparse.RawDescriptionHelpFormatter):    
        pass

    desc = 'Convert STARsolo mtx output to h5ad'
    epi = """DESCRIPTION:
    Convert STARsolo mtx output to h5ad.
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                    formatter_class=CustomFormatter)
    parser.add_argument('star_output', type=str, help='STARsolo output directory')
    parser.add_argument('--sample', type=str, required=True,
                        help='Sample name')
    parser.add_argument('--output-dir', type=str, default='mtx2h5ad_out',
                        help='Output directory')
    parser.add_argument('--feature-types', type=str, nargs="+",
                        default=FEATURE_TYPES, choices=FEATURE_TYPES, 
                        help='Feature types to include')
    parser.add_argument(
        '--missing-metadata', type=str, default="error", 
        choices=["error", "skip", "allow"],
        help="How do handle missing metadata?"
    )
    parser.add_argument(
        '--update-database', action="store_true", default=False, 
        help="Update the database?"
    )
    return parser.parse_args()

# functions
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

def build_velocyto_anndata(matrix_paths: Dict[str,str], feature_paths: str, barcode_paths: str) -> sc.AnnData:
    """
    Generate an anndata object from the STAR aligner output folder
    Args:
        matrix_paths: Path to matrix files, {matrix_type: path}
        feature_paths: Path to feature file
        barcode_paths: Path to barcode file
    Returns:
        Anndata object
    """
    # check exists
    for matrix_path in matrix_paths.values():
        if not os.path.exists(matrix_path):
            raise FileNotFoundError(f"{matrix_path} not found")

    # Transpose counts matrix to have Cells as rows and Genes as cols as expected by AnnData objects
    ## Using spliced.mtx as the reference matrix
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
    obs = pd.read_csv(barcode_paths, header = None, index_col = 0)

    # Remove index column name to make it compliant with the anndata format
    obs.index.name = None
    var = pd.read_csv(feature_paths, sep='\t', names = ('gene_ids', 'feature_types'), index_col = 1)
  
    # Build AnnData object to be used with ScanPy and ScVelo
    adata = anndata.AnnData(
        X = X, obs = obs, var = var,
        layers = {'spliced': spliced, 'unspliced': unspliced, 'ambiguous': ambiguous}
    )
    adata.var_names_make_unique()

    # Subset Cells based on STAR filtering
    selected_barcodes = pd.read_csv(barcode_paths, header = None)
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

def build_gene_anndata(
    mtx_filt: str, feat_filt: str, barcode_filt: str, 
    mtx_raw: Dict[str,str], barcode_raw: str
    ) -> sc.AnnData:
    """
    Generate an anndata object from the STAR aligner output folder
    Args:
        mtx_filt: Path to filtered matrix file
        feat_filt: Path to filtered feature file
        barcode_filt: Path to filtered barcode file
        mtx_raw: Path to raw matrix files
        barcode_raw: Path to raw barcode file
    Returns:
        Anndata object
    """
    # primary load count matrix
    logging.info("Loading primary count matrix...")
    adata = sc.read_10x_mtx(
        os.path.dirname(mtx_filt),
        var_names="gene_ids",
        make_unique=True,
        cache=False,
        gex_only=True
    )  
    # filtering multi-mapper count matrices
    logging.info("Adding multi-mapper count matrices as layers...")
    for matrix_type in ['UniqueAndMult-Uniform', 'UniqueAndMult-EM']:
        logging.info(f"Filtering {matrix_type} matrix...")
        match_barcodes(
            barcode_file = barcode_raw,
            matrix_file = mtx_raw[matrix_type],
            target_barcode_file = barcode_filt,
            out_cb_file = f'barcodes_{matrix_type}_filtered.tsv',
            out_mat_file = f'{matrix_type}_filtered.mtx'
        )
        adata.layers[matrix_type] = sc.read_mtx(f'{matrix_type}_filtered.mtx').X.transpose()
    return adata

def load_matrix_as_anndata(
        srx_id: str, 
        metadata: pd.DataFrame,
        feature_type: str,
        mtx_raw: Dict[str, str],
        barcode_raw: str,
        mtx_filt: Dict[str, str],
        feat_filt: str,
        barcode_filt: str,
        update_database: bool = False
    ) -> None:
    """
    Load a matrix.mtx file as an AnnData object.
    Args:
        srx_id: SRX accession
        metadata: Metadata for the SRX accession
        feature_type: Feature type
        mtx_raw: Path to raw matrix file
        barcode_raw: Path to raw barcode file
        mtx_filt: Path to filtered matrix file
        feat_filt: Path to filtered feature file
        barcode_filt: Path to filtered barcode file
        update_database: bool
    """
    # build anndata
    if all(x in mtx_filt for x in ["spliced", "unspliced", "ambiguous"]):
        logging.info("Building Velocyto anndata...")
        adata = build_velocyto_anndata(mtx_filt, feat_filt, barcode_filt)
    elif 'matrix' in mtx_filt:
        logging.info("Building gene anndata...")
        adata = build_gene_anndata(mtx_filt['matrix'], feat_filt, barcode_filt, mtx_raw, barcode_raw)
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
    outfile = os.path.join("metadata", f"{srx_id}.csv")
    adata.obs["cell_barcode"] = adata.obs.index
    adata.obs["organism"] = metadata["organism"].values[0]
    adata.obs.to_csv(outfile, index=False)

    # add feature type
    metadata["feature_type"] = feature_type

    # upsert metadata to postgresql database
    if update_database:
        logging.info(f"Upserting metadata for SRX accession {srx_id}...")
        with db_connect() as conn:
            db_upsert(metadata, "scbasecamp_metadata_tmp", conn)
    else:
        logging.info(f"Skipping upserting metadata for SRX accession {srx_id}")

def get_basename(path: str) -> str:
    # remove trailing .gz
    if path.endswith(".gz"):
        path = path[:-3]
    return os.path.basename(os.path.splitext(path)[0])

def main(args: argparse.Namespace, log_df: pd.DataFrame) -> Optional[None]:
    # get metadata
    metadata = get_metadata(args.sample, args.missing_metadata)
    if metadata is None:
        return None
    ## add tissue category to metadata
    #metadata = add_tissue_category(metadata, args.tissue_categories)

    # find target files
    for feat_type in args.feature_types:
        p = os.path.join(args.star_output, f"{feat_type}", "raw", "*.mtx.gz")
        mtx_raw = {get_basename(x): x for x in glob(p)}
        p = os.path.join(args.star_output, f"{feat_type}", "filtered", "*.mtx.gz")
        mtx_filt = {get_basename(x): x for x in glob(p)}

        #barcode_raw = {get_basename(x): x for x in glob(p)}
        barcode_raw = os.path.join(args.star_output, f"{feat_type}", "raw", "barcodes.tsv.gz")
        feat_filt = os.path.join(args.star_output, f"{feat_type}", "filtered", "features.tsv.gz")
        #feat_filt = {get_basename(x): x for x in glob(p)}
        barcode_filt = os.path.join(args.star_output, f"{feat_type}", "filtered", "barcodes.tsv.gz")
        #barcode_filt = {get_basename(x): x for x in glob(p)}
        

        # Convert to h5ad
        load_matrix_as_anndata(
            srx_id = args.sample, 
            metadata = metadata,
            feature_type = feat_type,
            mtx_raw = mtx_raw,
            barcode_raw = barcode_raw,
            mtx_filt = mtx_filt,
            feat_filt = feat_filt,
            barcode_filt = barcode_filt,
            update_database = args.update_database
        )


## script main
if __name__ == '__main__':
    args = parse_args()

    # setup
    os.makedirs(args.output_dir, exist_ok=True)
    log_df = pd.DataFrame(
        columns=["sample", "accession", "process", "step", "status", "message"]
    )

    # run main
    main(args, log_df)
    
    # upsert log to database
    with db_connect() as conn:
       db_upsert(log_df, "screcounter_log", conn)   
#!/usr/bin/env python
# import
## batteries
import os
import gc
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
from scipy.io import mmread
from pypika import Query, Table
## package
from db_utils import db_connect, db_upsert


# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
logging.getLogger("psycopg2").setLevel(logging.CRITICAL)
logging.getLogger("google.auth.transport.requests").setLevel(logging.CRITICAL)
logging.getLogger("urllib3").setLevel(logging.CRITICAL)
logging.getLogger("google.auth").setLevel(logging.CRITICAL)
logging.getLogger("scanpy").setLevel(logging.CRITICAL)
logging.getLogger("anndata").setLevel(logging.CRITICAL)
logging.getLogger("scipy").setLevel(logging.CRITICAL)


# argparse
FEATURE_TYPES = ["Gene", "GeneFull", "GeneFull_Ex50pAS", "GeneFull_ExonOverIntron", "Velocyto"]
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
    parser.add_argument(
        '--sample', type=str, required=True,
        help='Sample name'
    )
    parser.add_argument(
        '--output-dir', type=str, default='h5ad',
        help='Output directory'
    )
    parser.add_argument(
        '--keep-raw-h5ad', action="store_true", default=False,
        help='Keep raw h5ad files instead of just the filtered ones'
    )
    parser.add_argument(
        '--keep-input', action="store_true", default=False,
        help='Keep input files'
    )
    parser.add_argument(
        '--feature-types', type=str, nargs="+",
        default=FEATURE_TYPES, choices=FEATURE_TYPES,
        help='Feature types to include'
    )
    parser.add_argument(
        '--use-database', action="store_true", default=False, 
        help="Use the scRecounter SQL database?"
    )
    return parser.parse_args()

# functions
def get_metadata(srx_id: str) -> Optional[pd.DataFrame]:
    """
    Get metadata for an SRX accession.
    Args:
        srx_id: SRX accession
        use_database: Use the scRecounter SQL database?
    Returns:
        pd.DataFrame: Metadata for the SRX accession
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

    ## check metadata
    if metadata is None or metadata.shape[0] == 0:
        raise ValueError(f"Metadata not found for SRX accession {srx_id}")
    if metadata.shape[0] > 1:
        raise ValueError(f"Multiple metadata entries found for SRX accession {srx_id}")
    return metadata

def build_velocyto_anndata(matrix_paths: Dict[str,str], feature_paths: str, barcode_paths: str) -> sc.AnnData:
    """
    Generate an anndata object from the STAR aligner output folder using memory-efficient sparse matrix loading.
    Args:
        matrix_paths: Path to matrix files, {matrix_type: path}
        feature_paths: Path to feature file
        barcode_paths: Path to barcode file
    Returns:
        Anndata object
    """
    # check exists
    required_matrices = ['spliced', 'unspliced', 'ambiguous']
    for matrix_type in required_matrices:
        matrix_path = matrix_paths.get(matrix_type)
        if not matrix_path or not os.path.exists(matrix_path):
            raise FileNotFoundError(f"{matrix_type} matrix path not found or file missing: {matrix_path}")
    if not os.path.exists(feature_paths):
        raise FileNotFoundError(f"Feature file not found: {feature_paths}")
    if not os.path.exists(barcode_paths):
        raise FileNotFoundError(f"Barcode file not found: {barcode_paths}")

    # Load the 3 matrices containing Spliced, Unspliced and Ambigous reads
    # mmread reads as (features, cells), transpose to (cells, features) for AnnData
    try:
        logging.info("Reading spliced matrix...")
        spliced = mmread(matrix_paths['spliced']).astype('float32').transpose().tocsr()
        logging.info("Reading unspliced matrix...")
        unspliced = mmread(matrix_paths['unspliced']).astype('float32').transpose().tocsr()
        logging.info("Reading ambiguous matrix...")
        ambiguous = mmread(matrix_paths['ambiguous']).astype('float32').transpose().tocsr()
    except Exception as e:
        logging.error(f"Error reading MTX files: {e}")
        raise

    # Use spliced counts as the primary matrix X
    X = spliced

    # Load Genes and Cells identifiers
    try:
        obs = pd.read_csv(barcode_paths, header=None, index_col=0, names=['barcode'])
        obs.index.name = None # AnnData expects unnamed index for obs

        var = pd.read_csv(
            feature_paths, sep='\t', header=None, names=('gene_ids', 'feature_types'), index_col=0
        )
        var.index.name = None # AnnData expects unnamed index for var
    except Exception as e:
        logging.error(f"Error reading feature/barcode files: {e}")
        raise

    # Ensure var index matches matrix shape
    if X.shape[1] != len(var):
         raise ValueError(f"Shape mismatch: Matrix columns ({X.shape[1]}) != Feature count ({len(var)})")
    # Ensure obs index matches matrix shape
    if X.shape[0] != len(obs):
         raise ValueError(f"Shape mismatch: Matrix rows ({X.shape[0]}) != Barcode count ({len(obs)})")


    # Build AnnData object to be used with ScanPy and ScVelo
    try:
        adata = anndata.AnnData(
            X = X, obs = obs, var = var,
            layers = {'spliced': spliced, 'unspliced': unspliced, 'ambiguous': ambiguous}
        )
        adata.var_names_make_unique()
    except Exception as e:
        logging.error(f"Error creating AnnData object: {e}")
        raise

    # # Subset Cells based on STAR filtering (This seems redundant if barcode_paths already points to filtered barcodes)
    #selected_barcodes = pd.read_csv(barcode_paths, header = None)
    #return adata[selected_barcodes[0]]

    # Assuming the barcode_paths file provided already corresponds to the cells in the matrices
    return adata

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

def build_gene_anndata_filt(
    mtx_filt: str, barcode_filt: str, 
    mtx_raw: Dict[str,str], barcode_raw: str,
    keep_input: bool = False,
    ) -> sc.AnnData:
    """
    Generate an anndata object from the STAR aligner output folder.
    This function is used for filtered matrices.
    Args:
        mtx_filt: Path to filtered matrix file
        feat_filt: Path to filtered feature file
        barcode_filt: Path to filtered barcode file
        mtx_raw: Path to raw matrix files
        barcode_raw: Path to raw barcode file
        keep_input: bool
    Returns:
        Anndata object
    """
    # primary load count matrix
    logging.info("Loading primary count matrix...")
    adata = sc.read_10x_mtx(
        os.path.dirname(mtx_filt),
        var_names="gene_ids",
        make_unique=True
    )  
    # filtering multi-mapper count matrices
    logging.info("Adding multi-mapper count matrices as layers...")
    for matrix_type in ['UniqueAndMult-Uniform', 'UniqueAndMult-EM']:
        logging.info(f"Filtering {matrix_type} matrix...")
        # filter barcodes
        match_barcodes(
            barcode_file = barcode_raw,
            matrix_file = mtx_raw[matrix_type],
            target_barcode_file = barcode_filt,
            out_cb_file = f'barcodes_{matrix_type}_filtered.tsv',
            out_mat_file = f'{matrix_type}_filtered.mtx'
        )
        # add to anndata
        X = mmread(f'{matrix_type}_filtered.mtx').astype('float32')
        adata.layers[matrix_type] = sparse.csr_matrix(X).transpose()
        # delete filtered files
        if not keep_input:
            logging.info(f"  Deleting temporary filtered files for {matrix_type}...")
            os.remove(f'{matrix_type}_filtered.mtx')
            os.remove(f'barcodes_{matrix_type}_filtered.tsv')
        
    return adata

def build_gene_anndata_raw(
    mtx_raw: Dict[str,str], feat_raw: str, barcode_raw: str
    ) -> sc.AnnData:
    """
    Generate an anndata object from the STAR aligner output folder.
    This function is used for raw matrices.
    Args:
        mtx_raw: Path to raw matrix file
        feat_raw: Path to raw feature file
        barcode_raw: Path to raw barcode file
    Returns:
        Anndata object
    """
    # primary load count matrix
    logging.info("Loading primary count matrix...")
    adata = sc.read_10x_mtx(
        os.path.dirname(mtx_raw['matrix']),
        var_names="gene_ids",
        make_unique=True
    )  
    # Adding multi-mapper count matrices as layers
    logging.info("Adding multi-mapper count matrices as layers...")
    for matrix_type in ['UniqueAndMult-Uniform', 'UniqueAndMult-EM']:
        X = mmread(mtx_raw[matrix_type]).astype('float32')
        adata.layers[matrix_type] = sparse.csr_matrix(X).transpose()
    return adata

def load_matrix_as_anndata(
        srx_id: str, 
        metadata: pd.DataFrame,
        feature_type: str,
        output_dir: str,
        mtx_raw: Dict[str, str],
        barcode_raw: str,
        keep_input: bool = False,
        feat_raw: Optional[str] = None,
        mtx_filt: Optional[Dict[str, str]] = None,
        feat_filt: Optional[str] = None,
        barcode_filt: Optional[str] = None,
        use_database: bool = False,
    ) -> None:
    """
    Load a matrix.mtx file as an AnnData object.
    Args:
        srx_id: SRX accession
        metadata: Metadata for the SRX accession
        feature_type: Feature type
        mtx_raw: Path to raw matrix file
        barcode_raw: Path to raw barcode file
        feat_raw: Path to raw feature file
        mtx_filt: Path to filtered matrix file
        feat_filt: Path to filtered feature file
        barcode_filt: Path to filtered barcode file
        use_database: bool
        keep_input: bool
    """
    # build anndata
    if mtx_filt is None:
        logging.info("Processing raw matrix...")
        if all(x in mtx_raw for x in ["spliced", "unspliced", "ambiguous"]):
            logging.info("Building Velocyto anndata...")
            adata = build_velocyto_anndata(mtx_raw, feat_raw, barcode_raw)
        elif 'matrix' in mtx_raw:
            logging.info("Building gene anndata...")
            adata = build_gene_anndata_raw(mtx_raw, feat_raw, barcode_raw)
        else:
            x = ','.join(mtx_raw.keys())
            raise ValueError(f"Invalid matrix_paths. Available keys: {x}")
        out_prefix = "raw"
    else:
        logging.info("Processing filtered matrix...")
        if all(x in mtx_filt for x in ["spliced", "unspliced", "ambiguous"]):
            logging.info("Building Velocyto anndata...")
            adata = build_velocyto_anndata(mtx_filt, feat_filt, barcode_filt)
        elif 'matrix' in mtx_filt:
            logging.info("Building gene anndata...")
            adata = build_gene_anndata_filt(
                mtx_filt['matrix'], barcode_filt, 
                mtx_raw, barcode_raw, keep_input
            )
        else:
            x = ','.join(mtx_filt.keys())
            raise ValueError(f"Invalid matrix_paths. Available keys: {x}")
        out_prefix = "filtered"

    # drop 'feature_types' column in var
    if 'feature_types' in adata.var.columns:
        adata.var.drop(columns=['feature_types'], inplace=True)

    # list layers
    layers_str = ", ".join(list(adata.layers.keys()))
    logging.info(f"Layers: {layers_str}")

    # calculate total counts
    ## primary matrix
    if sparse.issparse(adata.X):
        adata.obs["gene_count"] = (adata.X > 0).sum(axis=1).A1
        adata.obs["umi_count"] = adata.X.sum(axis=1).A1
    else:
        adata.obs["gene_count"] = (adata.X > 0).sum(axis=1)
        adata.obs["umi_count"] = adata.X.sum(axis=1)
    ## layer matrices
    for layer in adata.layers: 
        if sparse.issparse(adata.layers[layer]):
            adata.obs[f"gene_count_{layer}"] = (adata.layers[layer] > 0).sum(axis=1).A1
            adata.obs[f"umi_count_{layer}"] = adata.layers[layer].sum(axis=1).A1
        else:
            adata.obs[f"gene_count_{layer}"] = (adata.layers[layer] > 0).sum(axis=1)
            adata.obs[f"umi_count_{layer}"] = adata.layers[layer].sum(axis=1)

    # add metadata to adata
    adata.obs["SRX_accession"] = srx_id

    ## write to h5ad
    h5ad_outdir = os.path.join(output_dir, out_prefix)
    os.makedirs(h5ad_outdir, exist_ok=True)
    outfile = os.path.join(h5ad_outdir, f"{feature_type}.h5ad")
    logging.info(f"Writing to {outfile}...")
    adata.write_h5ad(outfile, compression="gzip")

    # write out obs dataframe as csv
    # metadata_outdir = os.path.join(output_dir, out_prefix, "metadata")
    # os.makedirs(metadata_outdir, exist_ok=True)
    # outfile = os.path.join(metadata_outdir, f"{srx_id}.csv")
    # adata.obs["cell_barcode"] = adata.obs.index
    # if use_database:
    #     adata.obs["organism"] = metadata["organism"].values[0]
    # adata.obs.to_csv(outfile, index=False)

def get_basename(path: str) -> str:
    # remove trailing .gz
    if path.endswith(".gz"):
        path = path[:-3]
    return os.path.basename(os.path.splitext(path)[0])

def main(args: argparse.Namespace, log_df: pd.DataFrame) -> Optional[None]:
    # get metadata
    if args.use_database:
        metadata = get_metadata(args.sample)
    else:
        metadata = None

    # find target files
    for feat_type in args.feature_types:
        # matrices
        p = os.path.join(args.star_output, f"{feat_type}", "raw", "*.mtx.gz")
        mtx_raw = {get_basename(x): x for x in glob(p)}
        p = os.path.join(args.star_output, f"{feat_type}", "filtered", "*.mtx.gz")
        mtx_filt = {get_basename(x): x for x in glob(p)}
        # barcodes
        barcode_raw = os.path.join(args.star_output, f"{feat_type}", "raw", "barcodes.tsv.gz")
        barcode_filt = os.path.join(args.star_output, f"{feat_type}", "filtered", "barcodes.tsv.gz")
        # features
        feat_raw = os.path.join(args.star_output, f"{feat_type}", "raw", "features.tsv.gz")
        feat_filt = os.path.join(args.star_output, f"{feat_type}", "filtered", "features.tsv.gz")

        # Convert filtered matrices to h5ad
        load_matrix_as_anndata(
            srx_id = args.sample, 
            metadata = metadata,
            output_dir = args.output_dir,
            feature_type = feat_type,
            mtx_raw = mtx_raw,
            barcode_raw = barcode_raw,
            mtx_filt = mtx_filt,
            feat_filt = feat_filt,
            barcode_filt = barcode_filt,
            use_database = args.use_database,
            keep_input = args.keep_input
        )

        # delete the filtered matrix files
        if not args.keep_input:
            logging.info(f"Deleting input files for {feat_type} (filtered)...")
            for mtx in mtx_filt.values():
                if os.path.exists(mtx):
                    os.remove(mtx)
            if os.path.exists(barcode_filt):
                os.remove(barcode_filt)
            if os.path.exists(feat_filt):
                os.remove(feat_filt)

        # garbage collect
        del mtx_filt, barcode_filt, feat_filt
        gc.collect()
        
        # convert raw matrices to h5ad
        if args.keep_raw_h5ad:
            load_matrix_as_anndata(
                srx_id = args.sample, 
                metadata = metadata,
                output_dir = args.output_dir,
                feature_type = feat_type,
                mtx_raw = mtx_raw,
                barcode_raw = barcode_raw,
                feat_raw = feat_raw,
                use_database = args.use_database,
                keep_input = args.keep_input
            )

        # delete the raw mtx files
        if not args.keep_input:
            logging.info(f"Deleting input files for {feat_type} (raw)...")
            for mtx in mtx_raw.values():
                if os.path.exists(mtx):
                    os.remove(mtx)
            if os.path.exists(barcode_raw):
                os.remove(barcode_raw)
            # feat_raw might not exist for Velocyto
            if feat_raw and os.path.exists(feat_raw):  
                os.remove(feat_raw)

        # garbage collect
        del mtx_raw, barcode_raw, feat_raw
        gc.collect()


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
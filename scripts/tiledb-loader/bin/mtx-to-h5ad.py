#!/usr/bin/env python3
# import
## batteries
import os
import logging
import argparse 
import concurrent.futures
from pathlib import Path
from itertools import chain, repeat
from typing import List, Set, Tuple, Optional
## 3rd party
import numpy as np
import scipy.sparse
import pandas as pd
import tiledbsoma
import tiledbsoma.io
import scanpy as sc
from pypika import Query, Table
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
    desc = 'Convert mtx files to h5ad.'
    epi = """DESCRIPTION:
    Convert mtx files to h5ad in parallel.
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=CustomFormatter)
    parser.add_argument(
        '--srx', type=str, help="SRX accessions", required=True
    )
    parser.add_argument(
        '--mtx-path', type=str, help="Path to matrix.mtx.gz files", nargs='+', required=True
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
        '--threads', type=int, default=8, help="Number of threads to use"
    )
    return parser.parse_args()


def load_matrix_as_anndata(
        srx_id: str, 
        matrix_path: str, 
        missing_metadata: str="error",
    ) -> sc.AnnData:
    """
    Load a matrix.mtx.gz file as an AnnData object.
    Args:
        srx_id: SRX accession
        matrix_path: Path to matrix.mtx.gz file
        missing_metadata: How to handle missing metadata
    Returns:
        AnnData object
    """
    # get metadata from scRecounter postgresql database
    srx_metadata = Table("srx_metadata")
    stmt = (
        Query
        .from_(srx_metadata)
        .select(
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

    # load count matrix
    adata = sc.read_10x_mtx(
        os.path.dirname(matrix_path),
        var_names="gene_ids",
        make_unique=True,
        prefix=os.path.basename(matrix_path).split("_")[0] + "_"
    )

    # calculate total counts
    if scipy.sparse.issparse(adata.X):
        adata.obs["gene_count"] = (adata.X > 0).sum(axis=1).A1
        adata.obs["umi_count"] = adata.X.sum(axis=1).A1
    else:
        adata.obs["gene_count"] = (adata.X > 0).sum(axis=1)
        adata.obs["umi_count"] = adata.X.sum(axis=1)
    adata.obs["barcode"] = adata.obs.index

    # append SRX to barcode to create a global-unique index
    adata.obs.index = adata.obs.index + f"_{srx_id}"

    # add metadata to adata
    adata.obs["SRX_accession"] = srx_id
    for col in metadata.columns:
        try:
            adata.obs[col] = str(metadata[col].values[0])
        except IndexError:
            adata.obs[col] = None

    return adata

def mtx_to_h5ad(
    matrix_files: str, 
    missing_metadata: str="error",
    threads: int=8
    ) -> sc.AnnData:
    """
    Convert a list of matrix.mtx.gz files to a single h5ad file.
    Args:
        matrix_files: DataFrame with columns "srx" and "mtx_path"
        missing_metadata: How to handle missing metadata
        threads: Number of threads
    Returns:
        AnnData object
    """
    logging.info("Loading mtx files to h5ad...")

    # paralle load mtx files
    if threads == 1:
        adata = [load_matrix_as_anndata(x["srx"], x["mtx_path"], missing_metadata=missing_metadata) for _,x in matrix_files.iterrows()]
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            adata = list(executor.map(
                lambda x: load_matrix_as_anndata(
                    x[1]["srx"], x[1]["mtx_path"], missing_metadata=missing_metadata
                ), 
                matrix_files.iterrows()
        )   )
        ## filter out empty objects
        adata = [a for a in adata if a is not None]

    ## concat
    adata = sc.concat(adata, join="outer")

    ## write to h5ad
    adata.write_h5ad(f"data.h5ad")
    logging.info(f"Saved h5ad file to data.h5ad")

def parse_arg(arg: str) -> List[str]:
    """Parse a comma-separated argument into a list."""
    return [x.strip() for x in arg.lstrip("[").rstrip("]").split(",")]

def main():
    """Main function to run the TileDB loader workflow."""
    args = parse_arguments()

    # parse args
    srx_ids = parse_arg(args.srx)

    # combine srx and path
    srx_mtx = []
    for i in range(len(srx_ids)):
        srx_mtx.append([srx_ids[i], f"{i+1}_matrix.mtx.gz"])
    srx_mtx = pd.DataFrame(srx_mtx, columns=["srx", "mtx_path"])

    # create h5ad files
    mtx_to_h5ad(
        srx_mtx, 
        threads=args.threads,
        missing_metadata=args.missing_metadata
    )

if __name__ == "__main__":
    main()
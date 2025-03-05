#!/usr/bin/env python3
# import
## batteries
import os
import logging
import argparse 
from typing import List, Set, Tuple, Optional
## 3rd party
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
    desc = 'Convert mtx file to h5ad.'
    epi = """DESCRIPTION:
    Convert mtx files to h5ad
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=CustomFormatter)
    parser.add_argument(
        '--srx', type=str, help="SRX accessions", required=True
    )
    parser.add_argument(
        '--mtx-path', type=str, help="Path to matrix.mtx.gz file", required=True
    )
    parser.add_argument(
        '--missing-metadata', type=str, default="error", 
        choices=["error", "skip", "allow"],
        help="How do handle missing metadata?"
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

def main():
    """Main function to run the TileDB loader workflow."""
    args = parse_arguments()

    adata = load_matrix_as_anndata(
        args.srx,
        args.mtx_path,
        missing_metadata=args.missing_metadata
    ) 
    adata.write_h5ad(f"{args.srx}.h5ad")


if __name__ == "__main__":
    main()
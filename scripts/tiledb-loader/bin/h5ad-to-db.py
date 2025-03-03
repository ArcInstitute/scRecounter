#!/usr/bin/env python3
# import
## batteries
import os
import pickle
import logging
import argparse
from typing import List, Set, Tuple, Optional
## 3rd party
import pandas as pd
import tiledbsoma
import tiledbsoma.io
import scanpy as sc

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
    desc = 'Add scRNA-seq data to a TileDB database.'
    epi = """DESCRIPTION:
    If the database does not exist, it will be created. Otherwise, the data will be appended.
    The database collection structure is {feature_type}/{organism}/{experiment}.
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=CustomFormatter)
    parser.add_argument(
        'h5ad_path', type=str, help='Path to the h5ad file to load.'
    )
    parser.add_argument(
        '--registration-plan', type=str, help='Path to the TileDB registration plan (pkl) file', required=True
    )
    parser.add_argument(
        '--db-uri', type=str, help='URI of the TileDB database.', required=True
    )
    parser.add_argument(
        '--organism', type=str, help='Organism name.', required=True
    )
    parser.add_argument(
        '--feature-type', default='GeneFull_Ex50pAS', 
        choices=['Gene', 'GeneFull', 'GeneFull_Ex50pAS', 'GeneFull_ExonOverIntron', 'Velocyto'], 
        help='Feature type to process'
    )
    return parser.parse_args()

def append_to_database_from_mem(adata: sc.AnnData, db_uri: str, organism: str) -> None:
    """
    Append an AnnData object to the TileDB database.
    Args:
        db_uri: URI of the TileDB database
        adata: AnnData object to append
    """
    logging.info(f"  Appending data to {organism} collection...")

    # add organism to db_uri
    db_uri = os.path.join(db_uri, organism)

    # register AnnData objects
    rd = tiledbsoma.io.register_anndatas(
        db_uri,
        [adata],
        measurement_name="RNA",
        obs_field_name="obs_id",
        var_field_name="var_id",
    )

    # pickle the rd object
    with open("registration_data.pkl", "wb") as f:
        dump(rd, f)
        print(f"  Pickled registration data to {os.path.join(db_uri, 'registration_data.pkl')}")
    exit();

    ## resize the experiment
    with tiledbsoma.Experiment.open(db_uri) as exp:
        tiledbsoma.io.resize_experiment(
            exp.uri,
            nobs=rd.get_obs_shape(),
            nvars=rd.get_var_shapes()
        )

    # ingest new data into the db
    tiledbsoma.io.from_anndata(
        db_uri,
        adata,
        measurement_name="RNA",
        registration_mapping=rd,
    )


def create_tiledb_from_mem(adata: sc.AnnData, db_uri: str, organism: str) -> None:
    """
    Create a new tiledb database.
    Args:
        db_uri: URI of the TileDB database
        adata: AnnData object to append
    """
    # create collection
    with tiledbsoma.Collection.create(db_uri) as base_collection:
        print(f"Created base collection at {base_collection.uri}")
    db_uri = os.path.join(db_uri, organism)

    # create organism collection and experiment/measurement
    logging.info(f"  Creating collection for {organism}...")
    tiledbsoma.io.from_anndata(
        experiment_uri=db_uri, 
        anndata=adata,
        measurement_name="RNA",
    )

def load_tiledb_from_mem(h5ad_path: str, db_uri: str) -> None:
    """
    Load an AnnData object from memory and append to the TileDB database.
    Args:
        h5ad_path: Path to the h5ad file to load
        db_uri: URI of the TileDB database
    """
    # load anndata objects in parallel
    adata = sc.read_h5ad(h5ad_path)
    organism = adata.obs['organism'].unique()[0].replace(" ", "_")
    org_db_uri = os.path.join(db_uri, organism)
    
    # append to database
    if not os.path.exists(org_db_uri):
        create_tiledb_from_mem(adata, db_uri, organism)
    else:
        append_to_database_from_mem(adata, db_uri, organism)

def main():
    """Main function to run the TileDB loader workflow."""
    args = parse_arguments()

    # unpickle the registration plan
    with open(args.registration_plan, "rb") as inF:
        registration_plan = pickle.load(inF)

    # load h5ad file to db
    experiment_uri = os.path.join(args.db_uri, args.feature_type, args.organism)
    tiledbsoma.io.from_h5ad(
        experiment_uri,
        args.h5ad_path,
        measurement_name="RNA",
        obs_id_name="obs_id",
        var_id_name="feature_name",
        registration_mapping=registration_plan,
    )

if __name__ == "__main__":
    main()


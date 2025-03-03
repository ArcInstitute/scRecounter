#!/usr/bin/env python3
# import
## batteries
import os
import gc
import logging
import argparse
from pickle import dump
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
    desc = 'Register h5ad files for loading into a TileDB database.'
    epi = """DESCRIPTION:
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=CustomFormatter)
    parser.add_argument(
        'h5ad_files', type=str, help='Path to the h5ad file(s) to register.', nargs='+'
    )
    parser.add_argument(
        '--db-uri', type=str, help='URI of the TileDB database.', required=True
    )
    parser.add_argument(
        '--organism', type=str, help='Organism name.', required=True
    )
    parser.add_argument(
        '--matrix-type', type=str, help='Matrix type', required=True
    )
    parser.add_argument(
        '--feature-type', default='GeneFull_Ex50pAS', 
        choices=['Gene', 'GeneFull', 'GeneFull_Ex50pAS', 'GeneFull_ExonOverIntron', 'Velocyto'], 
        help='Feature type to process'
    )
    return parser.parse_args()

def create_db_OLD(db_uri, feature_type, organism, h5ad_path):
    # Check if the base collection exists, if not create it
    if not os.path.exists(db_uri):
        with tiledbsoma.Collection.create(db_uri) as base_collection:
            print(f"Created base collection at {db_uri}")
    else:
        base_collection = tiledbsoma.open(db_uri, mode="w")
        print(f"Opening existing base collection at {db_uri}")

    # Check if feature type collection exists, if not create it
    ft_uri = f"{db_uri}/{feature_type}"
    try:
        with base_collection.add_new_collection(feature_type) as ft_collection:
            print(f"Created feature type collection at {ft_collection.uri}")
    except tiledbsoma._exception.AlreadyExistsError:
        ft_collection = base_collection[feature_type]
        print(f"Opening existing feature type collection at {ft_uri}")

    # Using the AnnData schema, create an empty Experiment if it doesn't exist
    experiment_uri = f"{ft_uri}/{organism}"
    if organism not in ft_collection:
        try:
            tiledbsoma.io.from_h5ad(
                experiment_uri,
                h5ad_path,
                measurement_name="RNA",
                obs_id_name="obs_id",
                var_id_name="feature_name",
                ingest_mode="schema_only",
            )
            # Add the new Experiment to the Collection
            ft_collection[organism] = tiledbsoma.open(experiment_uri)
            print(f"Created Experiment at {experiment_uri}")
        except tiledbsoma.AlreadyExistsError:
            print(f"Experiment at {experiment_uri} already exists")
    else:
        print(f"Experiment at {experiment_uri} already exists")

def create_db_OLD(db_uri, feature_type, organism, h5ad_path):
    with tiledbsoma.Collection.create(db_uri) as base_collection:
        print(f"Created base collection at {db_uri}")

        with base_collection.add_new_collection(feature_type) as ft_collection:
            print(f"Created feature type collection at {ft_collection.uri}")

            # Using the AnnData schema, create an empty Experiment and add it to the feature type collection
            experiment_uri = f"{ft_collection.uri}/{organism}"
            tiledbsoma.io.from_h5ad(
                experiment_uri,
                h5ad_path,
                measurement_name="RNA",
                obs_id_name="obs_id",
                var_id_name="feature_name",
                ingest_mode="schema_only",
            )
            # Add the new Experiment, which was created by tiledbsoma.io.from_anndata, to the
            # Collection.
            ft_collection[organism] = tiledbsoma.open(experiment_uri)
            print(f"Created Experiment at {experiment_uri}")


def create_db(db_uri, feature_type, organism, h5ad_path):
    # create base collection
    try:
        base_collection = tiledbsoma.Collection.create(db_uri)
        print(f"Created base collection at {db_uri}")
    except tiledbsoma.AlreadyExistsError:
        base_collection = tiledbsoma.Collection.open(db_uri)
        print(f"Base collection exists, opened {db_uri}")

    # create feature type collection
    try:
        ft_collection = base_collection.add_new_collection(feature_type)
        print(f"Created feature type collection at {ft_collection.uri}")
    except (tiledbsoma.AlreadyExistsError, KeyError) as e:
        ft_collection = base_collection[feature_type]
        print(f"Feature type collection exists, opened {ft_collection.uri}")

    # Using the AnnData schema, create an empty Experiment for target organism and add it to the feature type collection
    experiment_uri = f"{ft_collection.uri}/{organism}"
    try:
        tiledbsoma.io.from_h5ad(
            experiment_uri,
            h5ad_path,
            measurement_name="RNA",
            obs_id_name="obs_id",
            var_id_name="feature_name",
            ingest_mode="schema_only",
        )

        ft_collection[organism] = tiledbsoma.open(experiment_uri)
        print(f"Created Experiment at {experiment_uri}") 
    except tiledbsoma._exception.SOMAError:
        print(f"Experiment at {experiment_uri} already exists")

    # return experiment uri
    return experiment_uri

def main():
    """Main function to run the TileDB loader workflow."""
    args = parse_arguments()

    # create database
    experiment_uri = create_db(args.db_uri, args.feature_type, args.organism, args.h5ad_files[0])

    # register h5ad files
    registration_plan = tiledbsoma.io.register_h5ads(
        experiment_uri,
        args.h5ad_files,
        measurement_name="RNA",
        obs_field_name="obs_id",
        var_field_name="feature_name",
    )

    # resize experiment
    tiledbsoma.io.resize_experiment(
        experiment_uri,
        nobs=registration_plan.get_obs_shape(),
        nvars=registration_plan.get_var_shapes(),
    )
    
    # pickle the registration plan
    with open("registration-plan.pkl", "wb") as outF:
        dump(registration_plan, outF)
    print("Pickled registration plan: registration-plan.pkl")


if __name__ == "__main__":
    main()


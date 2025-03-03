#!/usr/bin/env python3
# import
## batteries
import os
import logging
import argparse
import pickle
## 3rd party
import pandas as pd
import tiledbsoma
import tiledbsoma.io
import scanpy as sc

# Configure logging
logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.DEBUG)
logging.getLogger("tiledbsoma").setLevel(logging.WARNING)
logging.getLogger("tiledbsoma.io").setLevel(logging.WARNING)
logging.getLogger("tiledb").setLevel(logging.WARNING)


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    """Custom formatter for argparse to allow both default values and raw descriptions."""
    pass


def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.
    
    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    desc = "Register h5ad files for loading into a TileDB database."
    epi = "DESCRIPTION: Registers AnnData (.h5ad) files to a TileDB database."

    parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=CustomFormatter)
    parser.add_argument("h5ad_files", type=str, nargs="+", help="Path to the h5ad file(s) to register.")
    parser.add_argument("--db-uri", type=str, required=True, help="URI of the TileDB database.")
    parser.add_argument("--organism", type=str, required=True, help="Organism name.")
    parser.add_argument("--matrix-type", type=str, required=True, help="Matrix type.")
    parser.add_argument(
        "--feature-type",
        default="GeneFull_Ex50pAS",
        choices=["Gene", "GeneFull", "GeneFull_Ex50pAS", "GeneFull_ExonOverIntron", "Velocyto"],
        help="Feature type to process.",
    )
    
    return parser.parse_args()


def create_db(db_uri: str, feature_type: str, organism: str, h5ad_path: str) -> str:
    """
    Creates or retrieves a TileDB collection and registers an h5ad dataset.
    
    Args:
        db_uri (str): URI of the TileDB database.
        feature_type (str): The type of feature collection.
        organism (str): The organism for which the experiment is being created.
        h5ad_path (str): Path to the h5ad file.

    Returns:
        str: The URI of the created experiment.
    """

    # Create or open the base collection
    try:
        base_collection = tiledbsoma.Collection.create(db_uri)
        print(f"Created base collection at {db_uri}")
    except tiledbsoma.AlreadyExistsError:
        base_collection = tiledbsoma.Collection.open(db_uri)
        print(f"Base collection exists, opened {db_uri}")

    # Create or retrieve the feature type collection
    try:
        ft_collection = base_collection.add_new_collection(feature_type)
        print(f"Created feature type collection at {ft_collection.uri}")
    except (tiledbsoma.AlreadyExistsError, KeyError):
        ft_collection = base_collection[feature_type]
        ft_collection = tiledbsoma.open(ft_collection.uri, "w")  # Ensure write mode
        print(f"Feature type collection exists, opened {ft_collection.uri} in write mode")

    # Define the experiment URI
    experiment_uri = f"{ft_collection.uri}/{organism}"

    # Create the experiment if it does not already exist
    try:
        tiledbsoma.io.from_h5ad(
            experiment_uri,
            h5ad_path,
            measurement_name="RNA",
            obs_id_name="obs_id",
            var_id_name="feature_name",
            ingest_mode="schema_only",
        )

        # Add the experiment to the feature type collection
        ft_collection[organism] = tiledbsoma.open(experiment_uri, "w")
        print(f"Created Experiment at {experiment_uri}")
    except tiledbsoma._exception.SOMAError:
        print(f"Experiment at {experiment_uri} already exists")

    return experiment_uri


def main() -> None:
    """
    Main function to parse arguments, create the TileDB database, and register h5ad files.
    """
    # Parse command-line arguments
    args = parse_arguments()

    # Create or retrieve the database and experiment
    experiment_uri = create_db(args.db_uri, args.feature_type, args.organism, args.h5ad_files[0])

    # Register the h5ad files with the experiment
    registration_plan = tiledbsoma.io.register_h5ads(
        experiment_uri,
        args.h5ad_files,
        measurement_name="RNA",
        obs_field_name="obs_id",
        var_field_name="feature_name",
    )

    # Resize the experiment to accommodate the registered data
    tiledbsoma.io.resize_experiment(
        experiment_uri,
        nobs=registration_plan.get_obs_shape(),
        nvars=registration_plan.get_var_shapes(),
    )

    # Pickle the registration plan for later use
    with open("registration-plan.pkl", "wb") as outF:
        pickle.dump(registration_plan, outF)
    print("Pickled registration plan: registration-plan.pkl")


if __name__ == "__main__":
    main()
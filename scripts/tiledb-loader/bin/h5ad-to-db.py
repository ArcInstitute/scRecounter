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
    The database collection structure is {organism}/{experiment}.
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
        '--matrix-type', type=str, help='Matrix type', required=True
    )
    return parser.parse_args()

def main():
    """Main function to run the TileDB loader workflow."""
    args = parse_arguments()

    # unpickle the registration plan
    with open(args.registration_plan, "rb") as inF:
        registration_plan = pickle.load(inF)

    # Print the registration plan
    print("-- Registration plan --")
    print(registration_plan)

    # load h5ad file to db
    experiment_uri = os.path.join(args.db_uri, args.organism)
    print(f"Experiment URI: {experiment_uri}")
    tiledbsoma.io.from_h5ad(
        experiment_uri,
        args.h5ad_path,
        measurement_name=args.matrix_type,
        obs_id_name="obs_id",
        var_id_name="feature_name",
        registration_mapping=registration_plan,
    )

    # check results
    with tiledbsoma.open(args.db_uri) as db:
        experiment = db[args.organism]
        print(f"n_obs={experiment.obs.count}")
        print(f"n_vars={experiment.ms[args.matrix_type].var.count}")

if __name__ == "__main__":
    main()


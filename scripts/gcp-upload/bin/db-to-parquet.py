#!/usr/bin/env python3
# import
## batteries
import os
import logging
import argparse 
from uuid import uuid4
from pathlib import Path
from typing import List, Set, Tuple, Optional
## 3rd party
import pandas as pd
from pypika import Query, Table, Criterion
## package
from db_utils import db_connect

# format logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

# classes
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass

# functions
def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
    desc = 'Publish database results as parquet files.'
    epi = """DESCRIPTION:
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=CustomFormatter)
    parser.add_argument(
        '--feature-type', default='GeneFull_Ex50pAS', 
        choices=['Gene', 'GeneFull', 'GeneFull_Ex50pAS', 'GeneFull_ExonOverIntron', 'Velocyto'], 
        help='Feature type to process'
    )
    return parser.parse_args()


def load_scbasecamp_metadata(feature_type: str) -> pd.DataFrame:
    """
    Load metadata from scBasecamp database.
    Args:
        feature_type: Feature type to filter on.
    Returns:
        DataFrame with metadata.
    """
    logging.info("Obtaining scbasecamp metadata...")

    # get metadata from scRecounter postgresql database
    srx_metadata = Table("scbasecamp_metadata_tmp")
    stmt = (
        Query
        .from_(srx_metadata)
        .select("*")
        .where(srx_metadata.feature_type == feature_type)
    )
    with db_connect() as conn:
        metadata = pd.read_sql(str(stmt), conn)
    return metadata.drop(columns=['created_at', 'updated_at'])
    
def main():
    """Main function to run the TileDB loader workflow."""
    args = parse_arguments()

    # Load metadata
    metadata = load_scbasecamp_metadata(feature_type=args.feature_type)

    ## split by organism and save to parquet
    for organism, df in metadata.groupby('organism'):
        logging.info(f"Processing metadata for {organism}...")
        organism_str = organism.replace(" ", "_")
        # create directory
        out_dir = Path("metadata") / Path(args.feature_type) / Path(organism_str)
        out_dir.mkdir(parents=True, exist_ok=True)
        # write to parquet
        outfile = out_dir / 'sample_metadata.parquet'
        df.to_parquet(outfile, index=False, compression='gzip')
        logging.info(f"Saved metadata for {organism} to {outfile}")

if __name__ == "__main__":
    main()
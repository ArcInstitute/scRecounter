#!/usr/bin/env python
# import
## batteries
from __future__ import print_function
import os
import re
import sys
import argparse
import logging
import pandas as pd
from db_utils import db_connect, db_upsert

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

# argparse
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'Upload final STAR parameters to the scRecounter database'
epi = """DESCRIPTION:
Upload final STAR parameters to the scRecounter database.
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
parser.add_argument('--sample', type=str, default=None,
                    help='Sample name')
parser.add_argument('--barcodes', type=str, default=None,
                    help='Barcodes file path')
parser.add_argument('--star-index', type=str, default=None,
                    help='STAR index path')
parser.add_argument('--cell-barcode-length', type=int, default=None,
                    help='Cell barcode length')
parser.add_argument('--umi-length', type=int, default=None,
                    help='UMI length')
parser.add_argument('--strand', type=str, default=None,
                    help='Strandness')
parser.add_argument('--outfile', type=str, default="star_params.csv",
                    help='Output file path')
parser.add_argument('--use-database', action='store_true',
                    help='Use database to store STAR parameters')

# functions
def main(args):
    # set pandas display options
    pd.set_option('display.max_columns', 30)
    pd.set_option('display.width', 300)

    # create dataframe
    df = pd.DataFrame({
        'sample': [args.sample],
        'barcodes': [os.path.basename(args.barcodes)],
        'star_index': [os.path.basename(args.star_index.rstrip("/"))],
        'cell_barcode_length': [args.cell_barcode_length],
        'umi_length': [args.umi_length],
        'strand': [args.strand]
    })

    # write to file
    if os.path.exists(args.outfile):
        os.remove(args.outfile)
    df.to_csv(args.outfile, index=False)

    # upload to the scRecounter database
    if args.use_database:
        with db_connect() as conn:
            db_upsert(df, "screcounter_star_params", conn)

    # update screcounter log
    if args.use_database:
        log_df = pd.DataFrame({
            "sample": [args.sample],
            "accession": [""],
            "process": ["STAR save params"],
            "step": ["Final"],
            "status": ["Success"],
            "message": ["STAR final parameters saved to database"],
        })
        with db_connect() as conn:
            db_upsert(log_df, "screcounter_log", conn)
   

## script main
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
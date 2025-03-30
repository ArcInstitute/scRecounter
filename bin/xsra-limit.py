#!/usr/bin/env python
# import
from __future__ import print_function
import os
import sys
import json
import argparse
import logging
from glob import glob
from shutil import which, rmtree
from typing import Dict, Optional, List, Tuple
from subprocess import Popen, PIPE
import concurrent.futures
import pandas as pd
from db_utils import db_connect, db_upsert, add_to_log
from xsra_utils import rename_read_files, run_cmd, xsra_describe, xsra_dump, check_output

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

def parse_args():
    # argparse
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                          argparse.RawDescriptionHelpFormatter):    
        pass

    desc = 'Run xsra dump on an accession'
    epi = """DESCRIPTION:
    Run xsra dump on an accession.
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                    formatter_class=CustomFormatter)
    parser.add_argument('accession', type=str, help='Accession')
    parser.add_argument('--sample', type=str, default="",
                        help='Sample name')
    parser.add_argument('--threads', type=int, default=4,
                        help='Number of threads')
    parser.add_argument('--limit', type=int, default=None,
                        help='Maximum reads to write')                     
    parser.add_argument('--output-dir', type=str, default='prefetch_out',
                        help='Output directory')
    parser.add_argument('--min-read-length', type=int, default=26,
                        help='Minimum read length') 
    parser.add_argument('--use-database', action='store_true',
                        help='Use database to store STAR parameters')
    return parser.parse_args()

# functions
def main(args: argparse.Namespace, log_df: pd.DataFrame) -> Optional[None]:
    # check for executables
    for exe in ['xsra']:
        if not which(exe):
            raise OSError(f'{exe} not found in PATH')
            
    # run xsra describe
    read_names, msg = xsra_describe(
        args.accession, 
        min_read_length=args.min_read_length, 
        output_format="fasta"
    )
    add_to_log(log_df, args.sample, args.accession, "xsra", "describe", "Success", msg)

    # run xsra dump
    status,msg = xsra_dump(
        args.accession, 
        args.output_dir, 
        output_format="fasta", 
        limit=args.limit,
        threads=args.threads, 
    )
    add_to_log(log_df, args.sample, args.accession, "xsra", "dump", status, msg)

    # Check the xsra output and rename the files appropriately for each successful accession
    status, msg = check_output(read_names, args.accession, output_dir=args.output_dir)
    add_to_log(log_df, args.sample, args.accession, "xsra", "dump-check", status, msg)


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

    # write log to file
    log_df.to_csv(os.path.join(args.output_dir, "xsra.log"), index=False)
    
    # upsert log to database
    if args.use_database:
        with db_connect() as conn:
            db_upsert(log_df, "screcounter_log", conn)   
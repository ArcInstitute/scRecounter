#!/usr/bin/env python
# import
## batteries
from __future__ import print_function
import os
import re
import sys
import argparse
import logging
from typing import List, Dict, Any, Tuple
import pandas as pd
from db_utils import db_connect, db_upsert

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

# argparse
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'Summarize STAR summary files'
epi = """DESCRIPTION:
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
parser.add_argument('summary_csv', type=str, nargs='+',
                    help='STAR summary csv file(s)')
parser.add_argument('--sample', type=str, default="",
                    help='Sample name')
parser.add_argument('--accession', type=str, default="",
                    help='Accession number')
                    
# functions
def main(args):
    # set pandas display optionqs
    pd.set_option('display.max_columns', 50)
    pd.set_option('display.max_rows', 100)
    pd.set_option('display.width', 300)

    # read in all files and concatenate
    df = []
    regex = re.compile(r"_summary.csv$")
    for infile in args.summary_csv:
        x = pd.read_csv(infile, header=None)
        x.columns = ["category", "value"]
        x["feature"] = regex.sub("", os.path.basename(infile))
        df.append(x)
    df = pd.concat(df)

    # format category
    for x in ["Gene", "GeneFull", "GeneFull_Ex50pAS", "GeneFull_ExonOverIntron", "Velocyto"]:
        regex = re.compile(f" {x} ")
        df["category"] = df["category"].apply(lambda x: regex.sub(" feature ", x))

    # pivot table
    df = df.pivot(index='feature', columns='category', values='value').reset_index()
    
    # format columns: no spaces and lowercase
    df.columns = [x.lower().replace(" ", "_") for x in df.columns]
    
    # add sample and accession
    df["sample"] = args.sample
    df["accession"] = args.accession

    # coerce columns to numeric
    cols_to_convert = [
        "fraction_of_unique_reads_in_cells", "mean_gene_per_cell", "mean_umi_per_cell", 
        "mean_feature_per_cell", "median_gene_per_cell", "median_umi_per_cell",
        "median_feature_per_cell", "sequencing_saturation", "estimated_number_of_cells", 
        "number_of_reads", "reads_with_valid_barcodes", "umis_in_cells"
    ]
    for col in cols_to_convert:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    # float columns to integer
    cols_to_convert = ["estimated_number_of_cells", "number_of_reads", "reads_with_valid_barcodes", "umis_in_cells"]
    for col in cols_to_convert:
        if col in df.columns:
            df[col] = df[col].fillna(0).replace([float('inf'), -float('inf')], 0).astype(int)

    # upsert results to database
    with db_connect() as conn:
        db_upsert(df, "screcounter_star", conn)

## script main
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)

    
#!/usr/bin/env python
# import
## batteries
from __future__ import print_function
import os
import sys
import argparse
import logging
from typing import List, Dict, Any, Tuple
import pandas as pd
from db_utils import db_connect, db_upsert, add_to_log

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

# argparse
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'Select the best STAR parameters for each sample.'
epi = """DESCRIPTION:
The script reads the STAR summary CSV file and determines
the STAR parameters, based on the number of valid barcodes
for each parameter set among the test STAR runs.
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
parser.add_argument('read_stats_tsv', type=str,
                    help='Path to the seqkit stats file')
parser.add_argument('sra_stats_csv', type=str,
                    help='Path to the sra-stat file')
parser.add_argument('star_params_csv', type=str, nargs='+',
                    help='Path to the STAR parameters CSV file')
parser.add_argument('--outdir', type=str, default="results",
                    help='Path to the output JSON file')
parser.add_argument('--sample', type=str, default="",
                    help='Sample name')
parser.add_argument('--accession', type=str, default="",
                    help='SRA accession')
parser.add_argument('--reads_with_barcodes_cutoff', type=float, default=0.3,
                    help='Minimum fraction of reads with valid barcodes')

# functions
def read_seqkit_stats(stats_file: str, sample: str, accession: str) -> pd.DataFrame:
    """
    Read seqkit stats table and return as pandas dataframe.
    Args:
        stats_file: Path to seqkit stats file
        sample: Sample name
        accession: SRA accession
    Returns:
        pandas dataframe of seqkit stats
    """
    # check if stats_file is None
    if stats_file is None:
        return None
    # read in as pandas dataframe
    DF = pd.read_csv(stats_file, sep='\t')
    DF["accession"] = accession
    DF["read"] =  "read" + DF["file"].str.extract(r'read_([12]).fa.zst$') + "_length"
    DF = DF[["accession", "read", "avg_len"]]
    # convert avg_len to int
    DF["avg_len"] = DF["avg_len"].astype(int)
    # rename "avg_len" to "read_length"
    DF = DF.rename(columns={"avg_len": "read_length"})
    # pivot wider
    DF = DF.pivot(index='accession', columns='read', values='read_length').reset_index()
    DF["sample"] = sample
    return DF

def load_info(
    sra_stats_csv: str, star_params_csv: str, read_stats_tsv: str, sample: str, accession: str
    ) -> pd.DataFrame:
    """
    Load the information from the sra_stats_csv, star_params_csv, and read_stats_tsv
    and return the best parameters.
    Args:
        sra_stats_csv: Path to the sra_stats_csv file
        star_params_csv: Path to the star_params_csv file
        read_stats_tsv: Path to the read_stats_tsv file
        sample: Sample name
        accession: SRA accession
    Returns:
        pandas dataframe of parameter combinations
    """
    # read in sra stats file
    sra_stats = pd.read_csv(sra_stats_csv)

    # read in seqkit stats file
    seqkit_stats = read_seqkit_stats(read_stats_tsv, sample, accession)

    # read in param files
    params = pd.concat([pd.read_csv(x) for x in star_params_csv], axis=0)

    # merge on sample and accession, the add sra-stats to all records
    df = pd.merge(params, seqkit_stats, on=["sample", "accession"])
    return pd.merge(df, sra_stats, on=["accession"]) 


def get_strand_label(group: pd.DataFrame) -> str:
    """
    Get the strand label based on the number of reads mapped to the gene.
    Args:
        group: pandas dataframe group
    Returns:
        strand label
    """
    target_col = "Reads Mapped to GeneFull: Unique+Multiple GeneFull"
    # get the target column values for the strand
    try:
        fwd = group.loc[group["strand"]=="Forward",target_col].values[0]
    except IndexError:
        fwd = 0
    try:
        rev = group.loc[group["strand"]=="Reverse",target_col].values[0]
    except IndexError:
        rev = 0
    # return the strand label
    if fwd >= 2 * rev:
        return "Forward"
    elif rev >= 2 * fwd:
        return "Reverse"
    else:
        return "Ambiguous"

def get_best_params(
    data: pd.DataFrame, 
    reads_with_barcodes_cutoff: float=0.3
    ) -> pd.DataFrame:
    """
    Filter the data based on various criteria to select the best parameters.
    Args:
        data: pandas dataframe of all parameters
        reads_with_barcodes_cutoff: Minimum fraction of reads with valid barcodes
    Returns:
        pandas dataframe of best parameters
    """
    # group by
    group_by = ["sample", "accession", "barcodes_name", "star_index", "cell_barcode_length", "umi_length", "organism"]
    proper_strand = data.groupby(group_by).apply(get_strand_label, include_groups=False).reset_index(name="proper_strand")

    # join proper_strand to data
    data = pd.merge(data, proper_strand, on=group_by)

    # filter to proper strand
    data = data[data["strand"] == data["proper_strand"]].drop(columns="proper_strand")
    if data.shape[0] == 0:
        logging.info(f"No parameters passed the proper strand filter")

    # Filter on fraction of reads with valid barcode
    target_col = "Reads With Valid Barcodes"
    data = data[data[target_col] >= reads_with_barcodes_cutoff]
    if data.shape[0] == 0:
        logging.info(f"No parameters passed the reads_with_barcodes_cutoff filter: {reads_with_barcodes_cutoff}")

    # Filter to the max `Fraction of Unique Reads in Cells` => best parameters
    target_col = "Fraction of Unique Reads in Cells"
    data = data[data[target_col] != float("inf")]
    data = data[data[target_col] == data[target_col].max()]
    if data.shape[0] == 0:
        logging.info(f"No parameters passed the `max(Fraction of Unique Reads in Cells)` filter")

    return data

def write_all_data(data_all: pd.DataFrame, outfile_merged: str) -> None:
    """
    Write the merged data as CSV.
    Args:
        data_all: pandas dataframe of all parameters
        outfile_merged: Path to the merged parameters CSV
    """
    #-- table of all parameters --#
    # Estimate the number of cells
    data_all["saturation"] = data_all["Number of Reads"] / data_all["Sequencing Saturation"]
    data_all["num_spots"] = data_all["saturation"].where(data_all["spot_count"] > data_all["saturation"], data_all["spot_count"])
    data_all["num_cells"] = data_all["num_spots"] / data_all["Number of Reads"] * data_all["Estimated Number of Cells"] #* data_all["Reads Mapped to GeneFull: Unique+Multiple GeneFull"]
    data_all["Total Estimated Number of Cells"] = data_all["num_cells"].round().astype(int)
    data_all.drop(columns=["saturation", "num_spots", "num_cells"], inplace=True)

    # Write parameters as CSV
    data_all.to_csv(outfile_merged, index=False)
    logging.info(f"Output written to: {outfile_merged}")

def write_data(data: pd.DataFrame, data_all: pd.DataFrame, outfile_selected: str, outfile_merged: str) -> None:
    """
    Write the data as JSON and the merged data as CSV.
    Args:
        data: pandas dataframe of best parameters
        data_all: pandas dataframe of all parameters
        outfile_selected: Path to the selected parameters JSON file
        outfile_merged: Path to the merged parameters CSV
    """
    # write data as JSON
    if data is None:
         # write empty json file
        with open(outfile_selected, "w") as outF:
            outF.write("{}")
    else:
        # write data as JSON
        with open(outfile_selected, "w") as outF:
            outF.write(data.to_json(indent=4))
        # set best parameters for data_all
        target_cols = ["sample", "barcodes_file", "star_index", "cell_barcode_length", "umi_length", "strand"]
        df = data[target_cols].copy().to_frame().T
        df["Best parameters"] = True
        data_all = pd.merge(data_all, df, on=target_cols, how="left")
        data_all["Best parameters"] = data_all["Best parameters"].astype('boolean').fillna(False)
    logging.info(f"Output written to: {outfile_selected}")

    # write merged data as CSV
    write_all_data(data_all, outfile_merged)

def main(args, log_df):
    # set pandas display optionqs
    pd.set_option('display.max_columns', 40)
    pd.set_option('display.width', 300)
    process = "Select STAR params"

    # set output file paths
    outfile_selected = os.path.join(args.outdir, "selected_star_params.json")
    outfile_merged = os.path.join(args.outdir, "merged_star_params.csv")

    # load the data
    data_all = load_info(
        args.sra_stats_csv, args.star_params_csv, args.read_stats_tsv, 
        args.sample, args.accession
    )

    # Filter data for various criteria
    data_filt = get_best_params(
        data_all.copy(),  reads_with_barcodes_cutoff=args.reads_with_barcodes_cutoff
    )

    # check if data is empty
    if data_filt.shape[0] == 0:
        msg = "Best parameters not found"
        logging.warning(msg)
        add_to_log(log_df, args.sample, args.accession, process, "Get best params", "Failure", msg)
        # write empty json file
        write_data(None, data_all, outfile_selected, outfile_merged)
        return None
    else:
        add_to_log(log_df, args.sample, args.accession, process, "Get best params", "Success", "Valid barcodes found")

    # Convert dtypes
    for x in ["cell_barcode_length", "umi_length", "read1_length", "read2_length"]:
        data_filt[x] = data_filt[x].astype(int)

    # Check that read lengths are >= CELL_BARCODE_LENGTH + UMI_LENGTH
    data_filt["CHECK"] = data_filt["cell_barcode_length"] + data_filt["umi_length"] - data_filt["read1_length"]
    data_filt = data_filt[data_filt["CHECK"] <= 0]
    if data_filt.shape[0] == 0:
        msg = "No valid barcodes found in the STAR summary table after accounting for read lengths"
        logging.error(msg)
        add_to_log(log_df, args.sample, args.accession, process, "Read length filter", "Failure", msg)
        write_data(None, data_all, outfile_selected, outfile_merged)
        return None
    else:
        add_to_log(log_df, args.sample, args.accession, process, "Read length filter", "Success", "Read lengths > cell barcode + UMI")
    data_filt.drop(columns="CHECK", inplace=True)

    # If multiple rows, take the first after sorting by "READS_WITH_VALID_BARCODES"
    data_filt = data_filt.sort_values("Reads With Valid Barcodes", ascending=False).iloc[0]

    # Write output
    write_data(data_filt, data_all, outfile_selected, outfile_merged)

    # Add to log table
    add_to_log(log_df, args.sample, args.accession, process, "Final", "Success", "Best parameters selected")

## script main
if __name__ == '__main__':
    # arg parse
    args = parser.parse_args()

    # setup
    os.makedirs(args.outdir, exist_ok=True)
    log_df = pd.DataFrame(
        columns=["sample", "accession", "process", "step", "status", "message"]
    )

    # run main
    main(args, log_df)

    # write log
    #log_file = os.path.join(args.outdir, "select-star-params_log.csv")
    #log_df.to_csv(log_file, index=False)
    #logging.info(f'Log written to: {log_file}')
    
    # upsert log to database
    with db_connect() as conn:
        db_upsert(log_df, "screcounter_log", conn)
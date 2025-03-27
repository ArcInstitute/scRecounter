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
import pandas as pd
from db_utils import db_connect, db_upsert, add_to_log
from prefetch import prefetch_workflow

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

def parse_args():
    # argparse
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                          argparse.RawDescriptionHelpFormatter):    
        pass

    desc = 'Run sra-tools prefetch and xsra'
    epi = """DESCRIPTION:
    Run xsra on an accession.
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                    formatter_class=CustomFormatter)
    parser.add_argument('accession', type=str, help='Accession')
    parser.add_argument('--sample', type=str, default="",
                        help='Sample name')
    parser.add_argument('--threads', type=int, default=4,
                        help='Number of threads')
    parser.add_argument('--temp', type=str, default='TMP_FILES',
                        help='Temporary directory')
    parser.add_argument('--max-spot-id', type=int, default=None,
                        help='Maximum reads to write')
    parser.add_argument('--output-dir', type=str, default='prefetch_out',
                        help='Output directory')
    parser.add_argument('--min-read-length', type=int, default=28,
                        help='Minimum read length')  
    parser.add_argument('--output-format', type=str, choices=['fasta', 'fastq'], default='fastq',
                        help='Output format')  
    # prefetch parser
    # prefetch_parser = parser.add_argument_group('prefetch')
    # prefetch_parser.add_argument('--max-size-gb', type=int, default=300,
    #                     help='Max file size in Gb')
    # prefetch_parser.add_argument('--tries', type=int, default=3,
    #                     help='Number of tries to download')
    # prefetch_parser.add_argument('--gcp-download', action='store_true', default=False,
    #                     help='Obtain sequence data from SRA GCP mirror')
    return parser.parse_args()

# functions
def run_cmd(cmd: str) -> Tuple[int, bytes, bytes]:
    """
    Run sub-command and return returncode, output, and error.
    Args:
        cmd: Command to run
    Returns:
        tuple: (returncode, output, error)
    """
    cmd = [str(i) for i in cmd]
    logging.info(f'Running: {" ".join(cmd)}')
    try:
        p = Popen(cmd, stdout=PIPE, stderr=PIPE)
        output, err = p.communicate()
        return p.returncode, output, err
    except Exception as e:
        logging.error(f"Error running command: {str(e)}")
        return 1, b"", str(e).encode()

def rename_read_files(read_lens_filt: Dict[str, int], output_dir: str) -> Dict[str, str]:
    """
    Rename reads in `read_lens_filt` to:
    - 'read_1.fa.zstd' if there's only one file.
    - If two or more files exist:
      * Compare the top 2 by read length.
      * If lengths differ, rename the largest to 'read_2.fa.zstd' and second largest to 'read_1.fa.zstd'.
      * If lengths are the same, rename them by alphabetical order to 'read_1.fa.zstd' and 'read_2.fa.zstd'.
    Returns a dictionary { "R1": <path>, "R2": <path> } with the renamed files.
    """
    read_files_filt = {}
    num_files = len(read_lens_filt)

    # Handle only 1 file
    if num_files == 1:
        logging.info("Only one read file found; renaming to read_1.fa.zstd")
        old_name = list(read_lens_filt.keys())[0]
        new_name = os.path.join(output_dir, "read_1.fa.zstd")

        if old_name == new_name:
            raise ValueError(f"New fasta name is the same as old name: {new_name}")

        os.rename(old_name, new_name)
        read_files_filt["R1"] = new_name
        logging.info(f"Renamed {old_name} to {new_name}")
        return read_files_filt

    # Handle 2 or more files
    logging.info(">=2 read files found; picking two largest to rename")
    ## Sort by descending read length and then take the top 2
    sorted_by_length_desc = sorted(read_lens_filt.items(), key=lambda x: x[1], reverse=True)
    top_two = sorted_by_length_desc[:2]

    # If the two longest have the same read length, re-sort them by filename ascending
    if top_two[0][1] == top_two[1][1]:
        logging.info("Top 2 files have the same read length; renaming by alphabetical order")
        # Sort by filename (x[0]) ascending
        top_two.sort(key=lambda x: x[0])

        # Rename in ascending order:
        #  - first becomes read_1.fa.zstd
        #  - second becomes read_2.fa.zstd
        for i, (old_name, _) in enumerate(top_two, start=1):
            new_name = os.path.join(output_dir, f"read_{i}.fa.zstd")
            if old_name == new_name:
                raise ValueError(f"New fasta name is the same as old name: {new_name}")
            os.rename(old_name, new_name)
            read_files_filt[f"R{i}"] = new_name
            logging.info(f"Renamed {old_name} to {new_name}")
    else:
        logging.info("Top 2 files differ in read length; largest file => read_2, second => read_1")
        # Assign read_2 to the largest read, read_1 to the second largest
        for i, (old_name, _) in enumerate(top_two, start=1):
            read_num = 2 if i == 1 else 1
            new_name = os.path.join(output_dir, f"read_{read_num}.fa.zstd")

            if old_name == new_name:
                raise ValueError(f"New fasta name is the same as old name: {new_name}")

            os.rename(old_name, new_name)
            read_files_filt[f"R{read_num}"] = new_name
            logging.info(f"Renamed {old_name} to {new_name}")

    return read_files_filt

def write_log(logF, sample: str, accession: str, step: str, success: bool, msg: str) -> None:
    """
    Write skip reason to file.
    Args:
        logF: Log file handle
        sample: Sample name
        accession: SRA accession
        step: Step name
        success: Success status
        msg: Message
    """
    if len(msg) > 100:
        msg = msg[:100] + '...'
    logF.write(','.join([sample, accession, step, str(success), msg]) + '\n')

def xsra_describe(sra_file: str, min_read_length: int) -> Tuple[Optional[List[List[str]]], str]:
    """
    Run`xsra describe` to determine which, if any, of the reads are correct R1 and R2.
    Params:
        sra_file: SRA file
        min_read_length: Minimum read length
    Returns:
        A list of lists with the R1 and R2 read names, or None if no reads are correct
    """
    # xsra describe
    cmd = ["xsra", "describe", "--limit", "10000", sra_file]

    ## run command
    returncode, output, err = run_cmd(cmd)
    if returncode != 0:
        logging.warning(err)
        return None, "xsra describe command failed"

    # parse json output to access the stats field
    try:
        json_output = json.loads(output.decode())
    except json.JSONDecodeError:  
        logging.warning(err)
        return None, "Failed to parse JSON output"
    
    # determine which, if any are the reads correct R1 and R2
    read_lens_filt = {}
    for seg in json_output['stats']:
        if seg['mean_length'] >= min_read_length:
            read_lens_filt[seg['sid']] = seg['mean_length']
    # if not at least 2 reads pass the filter, return None
    if len(read_lens_filt) < 2:
        msg = f"Less than 2 reads pass the filter: {len(read_lens_filt)}"
        logging.warning(msg)
        return None,msg
    # R2 should be the largest read, while R1 should be the second largest
    sorted_by_length_desc = sorted(read_lens_filt.items(), key=lambda x: x[1], reverse=True)
    top_two = sorted_by_length_desc[:2]
    if top_two[0][1] == top_two[1][1]:
        logging.warning("Top 2 files have the same read length; renaming by alphabetical order")
        top_two.sort(key=lambda x: x[0])
    # return the read names
    return [
        [f"seg_{top_two[1][0]}.fa.zst", "read_1.fa.zst"],
        [f"seg_{top_two[0][0]}.fa.zst", "read_2.fa.zst"],
    ],"Successfully found paired-end reads via: xsra describe"

def xsra_dump(sra_file: str, output_dir: str, output_format: str, threads: int=1, max_spot_id: Optional[int]=None) -> Tuple[str, str]:
    """
    Run `xsra dump` to dump the reads.
    Params:
        sra_file: SRA file or accession
        output_dir: Output directory
        output_format: Output format
        threads: Number of threads
        max_spot_id: Maximum spot ID
    Returns:
        Tuple of (status, message)
    """
    # xsra dump
    if output_format == "fasta":
        output_format = "a"
    elif output_format == "fastq":
        output_format = "q"
    else:
        logging.warning(f"Invalid output format: {output_format}")
        return "Failure", f"Invalid output format: {output_format}"

    cmd = [
            "xsra", "dump",
            "--split", 
            "--compression", "z",
            "--format", output_format,
            "--threads", threads,
            "--outdir", output_dir,
        ]
    if max_spot_id and max_spot_id > 0:
        cmd += ["--limit", str(max_spot_id)]
    cmd.append(sra_file)

    ## run command
    returncode, output, err = run_cmd(cmd)
    if returncode == 0:
        msg = output.decode().split('\n')
    else:
        logging.warning(err)
        msg = err.decode().split('\n')
    msg = "; ".join([x for x in msg if x])        
    if msg == "":
        msg = "No output from xsra dump"
    ## add to log        
    status = "Success" if returncode == 0 else "Failure"
    return status, msg

def check_output(read_names: List[str], accession: str, output_dir: str) -> Tuple[str, str]:
    """
    Check the output of xsra dump.
    Args:
        read_names: List of read names
        accession: SRA accession
        output_dir: Output directory
    Returns:
        Tuple of (status, message)
    """
    logging.info(f"Checking output for {accession}")

    # list all files in output_dir
    out_files_str = ", ".join(glob(os.path.join(output_dir, "*")))
    logging.info(f"Files in output_dir: {out_files_str}")

    # rename the output files
    for old_name, new_name in read_names:
        old_path = os.path.join(output_dir, old_name)
        new_path = os.path.join(output_dir, new_name)
        if not os.path.exists(old_path):
            msg = f"Read file not found: {old_path}"
            logging.warning(msg)
            return "Failure", msg
        try:
            logging.info(f"Renaming {old_path} to {new_path}")  
            os.rename(old_path, new_path)
        except OSError as e:
            msg = f"Error renaming {old_path} to {new_path}: {str(e)}"
            logging.error(msg)
            return "Failure", msg

    # list output files
    read_files = glob(os.path.join(output_dir, f"read_*.fa.zst"))
    if not read_files:
        msg = f"No target read files found; files present: {out_files_str}"
        logging.warning(msg)
        return "Failure",msg

    # check that the files are not empty
    for read_file in read_files:
        if os.path.getsize(read_file) == 0:
            msg = f"Read file is empty: {read_file}"
            logging.warning(msg)
            return "Failure",msg

    # return success if all files are present and not empty
    return "Success","xsra dump successful"

def xsra_prefetch(accession: str, output_dir: str, threads: int) -> Tuple[str, str]:
    """
    Run `xsra prefetch` to prefetch the reads.
    Args:
        accession: SRA accession
        output_dir: Output directory
        threads: Number of threads
    Returns:
        Tuple of (status, message)
    """
    logging.info(f"Prefetching {accession}")
    cmd = f"xsra prefetch {accession} -o {output_dir} -p -t {threads}"
    return run_cmd(cmd)

def xsra_all(sra_file: str, output_dir: str, output_format: str, threads: int) -> Tuple[str, str]:
    """
    Run `xsra dump` to dump the reads.
    Args:
        sra_file: SRA file
        output_dir: Output directory
        output_format: Output format
        threads: Number of threads
    Returns:
        Tuple of (status, message)
    """
    # run fast(er)q-dump
    # sra_file = prefetch_workflow(
    #     sample=args.sample, 
    #     accession=args.accession, 
    #     log_df=log_df,
    #     max_size_gb=args.max_size_gb,
    #     gcp_download=args.gcp_download,
    #     tries=args.tries,
    #     outdir=os.path.join(args.temp, "prefetch")
    # )
    # if sra_file is None:
    #     return None

    # prefetch via `xsra prefetch`
    status,msg = xsra_prefetch(sra_file, output_dir, threads=threads)
    add_to_log(log_df, args.sample, args.accession, "xsra", "prefetch", status, msg)

    # run `xsra dump` to dump the reads
    status,msg = xsra_dump(sra_file, output_dir, output_format=output_format, threads=threads)
    add_to_log(log_df, args.sample, args.accession, "xsra", "dump", status, msg)

def xsra_limit(sra_file: str, output_dir: str, output_format: str, threads: int, max_spot_id: int) -> Tuple[str, str]:
    """
    Run `xsra dump` to dump the reads with a limit on the number of spots.
    Args:
        sra_file: SRA file
        output_dir: Output directory
        threads: Number of threads
        max_spot_id: Maximum spot ID
    Returns:
        Tuple of (status, message)
    """
    accession = os.path.splitext(os.path.basename(sra_file))[0]
    status,msg = xsra_dump(accession, output_dir, output_format=output_format, threads=threads, max_spot_id=max_spot_id)

def main(args: argparse.Namespace, log_df: pd.DataFrame) -> Optional[None]:
    # check for executables
    for exe in ['xsra', 'prefetch', 'vdb-dump']:
        if not which(exe):
            logging.error(f'{exe} not found in PATH')
            sys.exit(1)

    # run `xsra describe` to determine which, if any, of the reads are correct R1 and R2
    read_names,msg = xsra_describe(args.accession, args.min_read_length)
    if read_names is None:
        add_to_log(log_df, args.sample, args.accession, "xsra", "describe", "Failure", msg)
        return None

    # if args.max_spot_id, just dump reads
    if args.max_spot_id:
        xsra_limit(args.accession, args.output_dir, args.output_format, args.threads, args.max_spot_id)
    else:
        xsra_all(args.accession, args.output_dir, args.output_format, args.threads)

    # Check the xsra output and rename the files appropriately
    status,msg = check_output(read_names, args.accession, output_dir=args.output_dir)
    add_to_log(log_df, args.sample, args.accession, "xsra", "dump-check", status, msg)

    # remove temp files
    rmtree(args.temp, ignore_errors=True) 

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
    with db_connect() as conn:
       db_upsert(log_df, "screcounter_log", conn)   
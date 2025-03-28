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
from prefetch import prefetch_workflow

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

def parse_args():
    # argparse
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                          argparse.RawDescriptionHelpFormatter):    
        pass

    desc = 'Run xsra on a list of accessions'
    epi = """DESCRIPTION:
    Run xsra on a list of accessions.
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                    formatter_class=CustomFormatter)
    parser.add_argument('accessions', type=str, nargs='+', help='Accessions')
    parser.add_argument('--sample', type=str, default="",
                        help='Sample name')
    parser.add_argument('--threads', type=int, default=4,
                        help='Number of threads')
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

def xsra_describe(
    accession: str, min_read_length: int, output_format: str
) -> Tuple[Optional[List[List[str]]], str]:
    """
    Run`xsra describe` to determine which, if any, of the reads are correct R1 and R2.
    Params:
        accession: SRA accession
        min_read_length: Minimum read length
    Returns:
        A list of lists with the R1 and R2 read names, or None if no reads are correct
    """
    # xsra describe
    cmd = ["xsra", "describe", "--limit", "10000", accession]

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

    # determine file extension
    if output_format == "fasta":
        file_ext = ".fa.zst"
    elif output_format == "fastq":
        file_ext = ".fq.zst"
    else:
        raise ValueError(f"Invalid output format: {output_format}")
        
    # return the read names
    return {accession : [
        [f"seg_{top_two[1][0]}{file_ext}", f"read_1{file_ext}"],
        [f"seg_{top_two[0][0]}{file_ext}", f"read_2{file_ext}"],
    ]},"Successfully found paired-end reads via: xsra describe"

def xsra_dump(
    accession: str, output_dir: str, output_format: str, threads: int=1, max_spot_id: Optional[int]=None
    ) -> Tuple[str, str]:
    """
    Run `xsra dump` to dump the reads.
    Params:
        accession: SRA accession
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
    cmd += [accession]

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

def xsra_prefetch(accessions: List[str], output_dir: str, threads: int) -> Tuple[str, str]:
    """
    Run `xsra prefetch` to prefetch the reads.
    Args:
        accessions: List of SRA accessions
        output_dir: Output directory
        threads: Number of threads
    Returns:
        Tuple of (status, message)
    """
    logging.info(f"Prefetching {accessions}")
    project_id = os.getenv("GCP_PROJECT_ID")
    if not project_id:
        return "Failure", "GCP_PROJECT_ID environment variable not set"
    cmd = [
        "xsra", "prefetch", "--gcp-project-id", project_id, "--provider", "gcp", "--full-quality", str(threads)
    ] + accessions
    returncode, output, err = run_cmd(cmd)
    if returncode != 0:
        return "Failure", f"xsra prefetch failed: {err}"
    return "Success", f"xsra prefetch successful: {output}"
    

def xsra_all(accessions: List[str], output_dir: str, output_format: str, threads: int) -> Tuple[str, str]:
    """
    Run `xsra dump` to dump the reads.
    Args:
        accessions: List of SRA accessions
        output_dir: Output directory
        output_format: Output format
        threads: Number of threads
    Returns:
        Tuple of (status, message)
    """
    # prefetch via `xsra prefetch`
    status,msg = xsra_prefetch(accessions, output_dir, threads=threads)
    exit();
    add_to_log(log_df, args.sample, args.accession, "xsra", "prefetch", status, msg)

    # run `xsra dump` to dump the reads
    status,msg = xsra_dump(accessions, output_dir, output_format=output_format, threads=threads)
    add_to_log(log_df, args.sample, args.accession, "xsra", "dump", status, msg)

def xsra_limit(
    accessions: List[str], output_dir: str, output_format: str, threads: int, max_spot_id: int
    ) -> Tuple[str, str]:
    """
    Run `xsra dump` to dump the reads with a limit on the number of spots.
    Args:
        accessions: List of accessions
        output_dir: Output directory
        output_format: Output format
        threads: Number of threads
        max_spot_id: Maximum spot ID
    Returns:
        Tuple of (status, message)
    """
    for accession in accessions:
        status,msg = xsra_dump(
            accession, output_dir, 
            output_format=output_format, 
            threads=threads, 
            max_spot_id=max_spot_id
        )

def describe_accessions(accessions: List[str], min_read_length: int, output_format: str, 
                           threads: int, sample: str, log_df: pd.DataFrame) -> Optional[List[Tuple]]:
    """
    Process multiple accessions in parallel using ThreadPoolExecutor.
    
    Args:
        accessions: List of SRA accessions to process
        min_read_length: Minimum read length
        output_format: Output format (fastq or fasta)
        threads: Number of threads to use
        sample: Sample name for logging
        log_df: DataFrame to log results
        
    Returns:
        List of tuples (accession, read_names, msg) or None if all failed
    """
    def describe_accession(accession):
        read_names, msg = xsra_describe(accession, min_read_length, output_format)
        if read_names is None:
            add_to_log(log_df, sample, accession, "xsra", "describe", "Failure", msg)
            return None
        return (accession, read_names, msg)
    
    # Use ThreadPoolExecutor to parallelize processing
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        future_to_accession = {
            executor.submit(describe_accession, accession): accession for accession in accessions
        }
        for future in concurrent.futures.as_completed(future_to_accession):
            accession = future_to_accession[future]
            try:
                result = future.result()
                if result is not None:
                    results.append(result)
                else:
                    logging.warning(f"Failed to process accession {accession}")
            except Exception as exc:
                logging.error(f"Accession {accession} generated an exception: {exc}")
                add_to_log(log_df, sample, accession, "xsra", "describe", "Failure", str(exc))
    
    if not results:
        logging.error("All accessions failed to process")
        return None
    
    return results

def main(args: argparse.Namespace, log_df: pd.DataFrame) -> Optional[None]:
    # check for executables
    for exe in ['xsra', 'prefetch', 'vdb-dump']:
        if not which(exe):
            raise OSError(f'{exe} not found in PATH')
            
    # Run `xsra describe` on accessions in parallel
    results = describe_accessions(
        args.accessions,
        args.min_read_length,
        args.output_format,
        args.threads,
        args.sample,
        log_df
    )

    # dump reads
    if args.max_spot_id:
        # if args.max_spot_id, just dump reads
        xsra_limit(
            args.accessions, 
            output_dir=args.output_dir, 
            output_format=args.output_format, 
            threads=args.threads, 
            max_spot_id=args.max_spot_id
        )
    else:
        xsra_all(
            args.accessions, 
            output_dir=args.output_dir, 
            output_format=args.output_format, 
            threads=args.threads
        )

    # Check the xsra output and rename the files appropriately for each successful accession
    for accession, read_names, _ in results:
        status, msg = check_output(read_names, accession, output_dir=args.output_dir)
        add_to_log(log_df, args.sample, accession, "xsra", "dump-check", status, msg)

    # remove temp files
    #rmtree(args.temp, ignore_errors=True) 

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
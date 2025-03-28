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
    return parser.parse_args()

# functions
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
    logging.info(f"Prefetching {', '.join(accessions)}")

    # project_id = os.getenv("GCP_PROJECT_ID")
    # if not project_id:
    #     return "Failure", "GCP_PROJECT_ID environment variable not set"  
    # "xsra", "prefetch", "--gcp-project-id", project_id, "--provider", "gcp", "--full-quality", str(threads)
    
    cmd = ["xsra", "prefetch", "--full-quality"] + accessions
    returncode, output, err = run_cmd(cmd)
    if returncode != 0:
        return "Failure", f"xsra prefetch failed: {err}"
    return "Success", f"xsra prefetch successful: {output}"

def describe_accessions(
    accessions: List[str], sample: str, min_read_length: int, output_format: str, threads: int, log_df: pd.DataFrame
    ) -> Optional[List[Tuple]]:
    """
    Process multiple accessions in parallel using ThreadPoolExecutor.
    
    Args:
        accessions: List of SRA accessions to process
        sample: Sample name for logging
        min_read_length: Minimum read length
        output_format: Output format (fastq or fasta)
        threads: Number of threads to use
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
    for exe in ['xsra']:
        if not which(exe):
            raise OSError(f'{exe} not found in PATH')
            
    # Run `xsra describe` on accessions
    read_idx = describe_accessions(
        args.accessions,
        args.sample,
        min_read_length=args.min_read_length,
        output_format="fastq",
        threads=args.threads,
        log_df=log_df
    )

    # run `xsra prefetch` on accessions
    status,msg = xsra_prefetch(args.accessions, args.output_dir, threads=args.threads)
    for accession in args.accessions:
        add_to_log(log_df, args.sample, accession, "xsra", "prefetch", status, msg)

    # run `xsra dump` to dump the reads
    for accession in args.accessions:
        status,msg = xsra_dump(f"{accession}.sra", args.output_dir, output_format="fastq", threads=args.threads)
        add_to_log(log_df, args.sample, accession, "xsra", "dump", status, msg)

        # get read names
        read_names = [read_names for accession, read_names, _ in read_idx if accession == accession][0]
        check_output(read_names, accession, output_dir=args.output_dir, append=True)
        add_to_log(log_df, args.sample, accession, "xsra", "dump-check", status, msg)

        # delete {output_dir}/seg*.fq.zst files
        logging.info(f"Deleting temp read files")
        for file in glob(os.path.join(args.output_dir, f"seg*.fq.zst")):
            os.remove(file)

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
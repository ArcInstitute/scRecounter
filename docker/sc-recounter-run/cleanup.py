#!/usr/bin/env python3
import os
import argparse
from typing import Tuple, List
from google.cloud import storage

# argparse
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'Clean up after a scRecounter production run'
epi = """DESCRIPTION:
Examples:
python cleanup.py gs://arc-ctc-nextflow/scRecounter/prod/work/SCRECOUNTER_2025-01-06_15-46-04/ gs://arc-ctc-screcounter/prod/SCRECOUNTER_2025-01-06_15-46-04/ 
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
parser.add_argument(
    'work_dir', type=str, 
    help='GCP bucket path to work directory (e.g., gs://bucket-name/path/to/folder)'
)
parser.add_argument(
    'output_dir', type=str, 
    help='GCP bucket path to output directory (e.g., gs://bucket-name/path/to/folder)'
)


# functions
def list_directories(bucket_name: str, prefix: str) -> List[str]:
    """
    List directories in a GCP bucket path
    Args:
        bucket_name: GCP bucket name
        prefix: GCP bucket prefix    
    Returns:
        List of directories in the bucket path
    """
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blobs = bucket.list_blobs(prefix=prefix, delimiter='/')
    
    directories = []
    for page in blobs.pages:
        directories.extend(page.prefixes)
    return [d.rstrip('/') for d in directories]

def delete_bucket_path(bucket_name: str, path: str) -> None:
    """
    Delete all objects in a GCP bucket path
    Args:
        bucket_name: GCP bucket name
    """

    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blobs = bucket.list_blobs(prefix=path)
    
    for blob in blobs:
        blob.delete()

def parse_gs_path(gs_path: str) -> Tuple[str, str]:
    """
    Parse a GCP bucket path
    Args:
       gs_path: GCP bucket path
    Returns:
        Tuple of bucket name and prefix
    """
    if not gs_path.startswith("gs://"):
        raise ValueError("Path must start with 'gs://'")
    parts = gs_path[5:].split("/", 1)
    bucket_name = parts[0]
    prefix = parts[1] if len(parts) > 1 else ""
    return bucket_name, prefix.rstrip("/") + "/"

def clean_output_dir(output_dir: str) -> None:
    """
    Delete the contents of the output directory, 
    if it only contains 'nf-report', 'nf-trace', and optionally 'logs' directories
    Args:
       output_dir: GCP bucket path to output directory
    """
    # parse the bucket path  
    bucket_name, path_prefix = parse_gs_path(output_dir)
    
    # list directories in the bucket path
    directories = list_directories(bucket_name, path_prefix)
    directories = [os.path.basename(d) for d in directories]
    print(f"Directories found: {', '.join(directories)}")
    
    allowed_dirs = {"nf-report", "nf-trace", "logs"}
    if set(directories).issubset(allowed_dirs):
        print("No STAR results found. Deleting the bucket path...")
        delete_bucket_path(bucket_name, path_prefix)
        print(f"Deleted path: {output_dir}")
    else:
        print("Bucket path contains other directories. No deletion performed.")

def clean_work_dir(work_dir:  str) -> None:
    """
    Delete the contents of the work directory
    Args:
       work_dir: GCP bucket path to work directory
    """
    # parse the bucket path  
    bucket_name, path_prefix = parse_gs_path(work_dir)
    
    print("Deleting the contents of the working directory...")
    delete_bucket_path(bucket_name, path_prefix)
    print(f"Deleted path: {work_dir}")

def main(args): 
    clean_work_dir(args.work_dir)
    clean_output_dir(args.output_dir)
    

if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
#!/usr/bin/env python3
import os
import sys
import argparse
from typing import Tuple, List, Dict
import pandas as pd
from google.cloud import storage
from psycopg2.extensions import connection
from db_utils import db_connect, db_update


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter): pass

def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.
    Returns:
        argparse.Namespace containing arguments.
    """
    desc = 'Purge SRX accessions from the scRecounter system.'
    epi = """DESCRIPTION:
    Purging:
     - Removes SRX records from scRecounter SQL database.
     - Removes the SRX directories from the GCP output folder of the scRecounter pipeline.

    Note: only scRecounter is purged, not SRAgent.

    Examples:
    purge-srx.py ERX10024831 ERX10086874
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=CustomFormatter)
    parser.add_argument('srx_accession', type=str, nargs='+',
                        help='>=1 SRX accession to purge from the scRecounter system.')
    parser.add_argument('--dry-run', action='store_true', default=False,
                        help='Print actions without executing.')
    parser.add_argument('--gcs-dir', type=str, default='gs://arc-ctc-screcounter/prod3/',
                        help='Base directory in GCP bucket where SCRECOUNTER directories are stored.') 
    parser.add_argument('--tenant', type=str, default='prod',
                        choices=['prod', 'test'],
                        help='SQL database tenant')                
    return parser.parse_args()

def parse_gs_path(gs_path: str) -> Tuple[str, str]:
    """
    Parse a GCP bucket path.
    Args:
        gs_path: GCP bucket path starting with gs://
    Returns:
        A tuple of (bucket_name, prefix).
    """
    if not gs_path.startswith("gs://"):
        raise ValueError("Path must start with 'gs://'")
    parts = gs_path[5:].split("/", 1)
    bucket_name = parts[0]
    prefix = parts[1] if len(parts) > 1 else ""
    return bucket_name, prefix.rstrip("/") + "/"

def list_screcounter_directories(
    bucket: storage.bucket.Bucket,
    prefix: str,
    srx_accesions: List[str],
) -> Dict[str,str]:
    """
    List directories named 'SCRECOUNTER_YYYY-MM-DD_hh-mm-ss' in the bucket under the given prefix.
    Args:
        bucket: The GCS bucket object.
        prefix: The prefix (subfolder) in which to look for SCRECOUNTER directories.
    Returns:
        A dictionary of {srx_accession: directory_path} for the target SRX accessions.
    """
    print(f"Searching for SRX directories...", file=sys.stderr)
    srx_dirs = {}
    for blob in bucket.list_blobs(prefix=prefix):
        blob_dir = os.path.dirname(blob.name)
        blob_dir_base = os.path.basename(blob_dir)
        blob_dir_parent = os.path.basename(os.path.dirname(blob_dir))
        if blob_dir_parent == "STAR" and blob_dir_base in srx_accesions:
            srx_dirs[blob_dir_base] = blob_dir
    print(f"  Found {len(srx_dirs)} SRX directories", file=sys.stderr)
    return srx_dirs

def purge_accession_tables(
    srx_dirs: Dict[str,str], bucket: storage.bucket.Bucket, dry_run: bool=False
    ) -> None:
    """
    Purge SRX accessions from the accession tables in the GCP bucket.
    Args:
        srx_dirs: Dictionary of {srx_accession: directory_path} for the target SRX accessions.
        bucket: The GCS bucket object.
        dry_run: If True, only print actions without executing.
    """   
    if len(srx_dirs) == 0:
        return None
    print(f"Purging accession tables...", file=sys.stderr)
    target_parent_dirs = set()
    for srx, srx_dir in srx_dirs.items():
        target_parent_dirs.add(os.path.dirname(os.path.dirname(srx_dir)))
    
    for parent_dir in target_parent_dirs:
        for blob in bucket.list_blobs(prefix=parent_dir):
            if os.path.basename(blob.name) == "accessions.csv":
                # read in accessions file
                if not dry_run:
                    df = pd.read_csv(pd.io.common.StringIO(blob.download_as_text()))
                    # filter out the SRX accessions
                    df = df[~df["sample"].isin(srx_dirs.keys())]
                    # write back to GCP
                    blob.upload_from_string(df.to_csv(index=False))
                print(f"  Purged {blob.name}", file=sys.stderr)

def delete_srx(srx_accessions: List[str], conn: connection, dry_run: bool=False):
    """
    Delete SRX accessions from scRecounter tables
    Args:
        srx_accessions: list of SRX accessions to delete
        conn: database connection
        dry_run: if True, only print actions without executing
    """
    if len(srx_accessions) == 0:
        return None
    print("Purging SRX accessions from scRecounter DB tables...", file=sys.stderr)
    target_tables = ["screcounter_log", "screcounter_star_params", "screcounter_star_results"]
    with db_connect() as conn:
        for srx in srx_accessions:
            if not dry_run:
                for tbl_name in target_tables:
                    with conn.cursor() as cur:
                        cur.execute(f"DELETE FROM {tbl_name} WHERE sample = '{srx}'")
                        conn.commit()
            print(f"  Deleted: {srx}", file=sys.stderr)

def delete_srx_star_dirs(srx_dirs: Dict[str,str], bucket: storage.bucket.Bucket, dry_run: bool=False):
    """
    Delete SRX directories from the GCP bucket
    Args:
        srx_dirs: Dictionary of {srx_accession: directory_path} for the target SRX accessions.
        bucket: The GCS bucket object.
        dry_run: If True, only print actions without executing.
    """
    if len(srx_dirs) == 0:
        return None
    print(f"Deleting SRX STAR directories...", file=sys.stderr)
    for srx_dir in srx_dirs.values():
        print(f"  Deleting: {srx_dir}", file=sys.stderr)
        if not dry_run:
            for blob in bucket.list_blobs(prefix=srx_dir):
                blob.delete()

def main(args: argparse.Namespace) -> None:
    """
     - Input: 
       - >=1 NCBI SRX accession
       - Path to GCP directory
     - Method
       - Recursively search in the GCP directory for folders named the same as the SRX accession
         - Delete each target folder
       - For each target folder, find the "accessions.csv" file 2 levels up from the target folder
       - Also delete the SRX from all scRecounter tables in the SQL database
    """
    os.environ["GCP_SQL_DB_NAME"] = f"sragent-{args.tenant}"
    print(f"GCP_SQL_DB_NAME: {os.getenv('GCP_SQL_DB_NAME')}", file=sys.stderr)

    # Parse the GCP bucket path
    bucket_name, path_prefix = parse_gs_path(args.gcs_dir)
    
    # Initialize GCP client and bucket
    client = storage.Client()
    bucket = client.bucket(bucket_name)

    # Dind target SRX directories in GCP bucket
    srx_dirs = list_screcounter_directories(bucket, path_prefix, args.srx_accession)

    # Selete SRX accessions from scRecounter tables
    purge_accession_tables(srx_dirs, bucket, dry_run=args.dry_run)

    # Selete SRX directories from GCP bucket
    delete_srx_star_dirs(srx_dirs, bucket, dry_run=args.dry_run)

    # Selete SRX accessions from scRecounter tables
    with db_connect() as conn:
        delete_srx(args.srx_accession, conn, dry_run=args.dry_run)


if __name__ == "__main__":
    from dotenv import load_dotenv
    load_dotenv()
    args = parse_args()
    main(args)
#!/usr/bin/env python3
import argparse
import sys
import os
import logging
import gzip
from typing import TextIO, Union

# Set up logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

# Custom formatter for argparse
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass

def parse_arguments():
    desc = 'Filter cells in a matrix to match target barcodes'
    epi = """DESCRIPTION:
    This script filters a matrix to keep only cells that match barcodes in the target file.
    The script supports gzipped input files (.gz extension).
    
    Example:
    # Filter matrix to match target barcodes
    ./filter-mtx.py barcodes.txt matrix.mtx target_barcodes.tsv --output-matrix matrix_filt.mtx --output-barcodes barcodes_filt.tsv

    ../../scripts/gcp-upload/bin/filter-mtx.py ./raw/barcodes.tsv ./raw/UniqueAndMult-EM.mtx ./filtered/barcodes.tsv
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=CustomFormatter)
    # Required arguments
    parser.add_argument('barcode_file', help='Input barcode file (can be gzipped)')
    parser.add_argument('matrix_file', help='Input matrix file (can be gzipped)')
    parser.add_argument('target_barcodes', help='Target barcodes file to match with (can be gzipped)')
    # Optional output file paths
    parser.add_argument('--output-matrix', default='matrix_filt.mtx', help='Output matrix file path')
    parser.add_argument('--output-barcodes', default='barcodes_filt.tsv', help='Output barcodes file path')
    # Parse arguments
    return parser.parse_args()

def open_file(filename: str, mode: str = 'r') -> Union[TextIO, gzip.GzipFile]:
    """Open a file, handling gzip if the filename ends with .gz"""
    if filename.endswith('.gz'):
        return gzip.open(filename, mode + 't')  # Text mode for gzip
    else:
        return open(filename, mode)

def match_barcodes(barcode_file: str, matrix_file: str, target_barcode_file: str, 
                  out_cb_file: str, out_mat_file: str) -> None:
    """
    Filter matrix to keep only cells that match target barcodes
    
    Args:
        barcode_file: Path to input barcode file
        matrix_file: Path to input matrix file
        target_barcode_file: Path to target barcodes file
        out_cb_file: Path to output barcodes file
        out_mat_file: Path to output matrix file
    """
    
    logger.info(f"Starting barcode matching")
    logger.info(f"  Input barcode file: {barcode_file}")
    logger.info(f"  Input matrix file: {matrix_file}")
    logger.info(f"  Target barcodes file: {target_barcode_file}")
    
    # Read original barcodes and create mapping
    barcode_to_index = {}
    with open_file(barcode_file) as f:
        for i, line in enumerate(f, 1):
            barcode = line.strip()
            barcode_to_index[barcode] = i
    
    logger.info(f"  Found {len(barcode_to_index)} original barcodes")
    
    # Read target barcodes
    target_barcodes = set()
    with open_file(target_barcode_file) as f:
        for line in f:
            target_barcodes.add(line.strip())
    
    logger.info(f"  Found {len(target_barcodes)} target barcodes (used to filter)")
    
    # Find matching barcodes and create new index mapping
    matched_indices = {}  # Maps original index to new index
    matched_barcodes = []  # List of matched barcodes
    
    for barcode, orig_idx in sorted(barcode_to_index.items(), key=lambda x: x[1]):
        if barcode in target_barcodes:
            new_idx = len(matched_indices) + 1
            matched_indices[orig_idx] = new_idx
            matched_barcodes.append(barcode)
    
    logger.info(f"  Found {len(matched_barcodes)} matching barcodes")
    
    # Write matched barcodes
    with open(out_cb_file, 'w') as out_f:
        for barcode in matched_barcodes:
            out_f.write(f"{barcode}\n")
    
    # Process matrix file
    header_lines = []
    filtered_entries = []
    n_features = 0
    
    with open_file(matrix_file) as f:
        # Process header lines (comments)
        line = f.readline().strip()
        while line.startswith('%'):
            header_lines.append(line)
            line = f.readline().strip()
        
        # Process dimensions line
        dimensions = line.split()
        n_features = int(dimensions[0])
        
        # Process data lines
        for line in f:
            cols = line.strip().split()
            gene_idx = cols[0]
            cell_idx = int(cols[1])
            umi_count = cols[2]
            
            if cell_idx in matched_indices:
                filtered_entries.append((gene_idx, matched_indices[cell_idx], umi_count))
    
    # Write filtered matrix
    with open(out_mat_file, 'w') as out_f:
        # Write header comments
        for line in header_lines:
            out_f.write(f"{line}\n")
        
        # Write new dimensions
        out_f.write(f"{n_features} {len(matched_barcodes)} {len(filtered_entries)}\n")
        
        # Write filtered entries
        for gene_idx, cell_idx, umi_count in filtered_entries:
            out_f.write(f"{gene_idx} {cell_idx} {umi_count}\n")
    
    logger.info(f"  Filtered barcodes saved to {out_cb_file}")
    logger.info(f"  Filtered matrix saved to {out_mat_file}")
    logger.info(f"  Barcode matching complete: {len(matched_barcodes)} cells, {len(filtered_entries)} matrix entries")

def main():
    args = parse_arguments()
    
    # Set default output file paths if not specified
    out_mat_file = args.output_matrix if args.output_matrix else f"{args.matrix_file}.filtered"
    out_cb_file = args.output_barcodes if args.output_barcodes else f"{args.barcode_file}.filtered"
    
    # Run barcode matching
    match_barcodes(
        barcode_file=args.barcode_file,
        matrix_file=args.matrix_file,
        target_barcode_file=args.target_barcodes,
        out_cb_file=out_cb_file,
        out_mat_file=out_mat_file
    )

if __name__ == "__main__":
    main()
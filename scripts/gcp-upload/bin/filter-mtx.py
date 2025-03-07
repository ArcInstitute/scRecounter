#!/usr/bin/env python3
import argparse
import gzip
import os

class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Filter cell barcodes and matrix based on UMI counts'
    )
    parser.add_argument('barcode_file', help='Input barcode file (can be gzipped)')
    parser.add_argument('matrix_file', help='Input matrix file (can be gzipped)')
    parser.add_argument('--exactCells', type=int, default=0, 
                        help='Exact number of cells to keep')
    parser.add_argument('--maxCells', type=int, default=3000, 
                        help='Maximum number of cells to consider')
    parser.add_argument('--maxPercentile', type=float, default=0.99, 
                        help='Percentile threshold for filtering')
    parser.add_argument('--maxMinRatio', type=float, default=10, 
                        help='Ratio between max and min UMI counts')
    return parser.parse_args()

def is_gzipped(filename: str) -> bool:
    """Check if a file is gzipped based on its extension."""
    return filename.endswith('.gz')

def open_file(filename: str, mode: str='r') -> str:
    """Open a file, gzipped or not, in the appropriate mode."""
    if is_gzipped(filename):
        return gzip.open(filename, mode + 't')  # 't' for text mode
    else:
        return open(filename, mode)

def main():
    args = parse_arguments()
    
    # Set parameters
    exactCells = args.exactCells
    maxCells = args.maxCells
    maxPercentile = args.maxPercentile
    maxMinRatio = args.maxMinRatio
    
    # Print parameters
    print(f"Parameters: exactCells={exactCells}, maxCells={maxCells}, "
          f"maxPercentile={maxPercentile}, maxMinRatio={maxMinRatio}")
    
    # Define constants
    n_header_lines = 3
    
    # Output file paths - remove .gz extension if present
    base_barcode = args.barcode_file[:-3] if is_gzipped(args.barcode_file) else args.barcode_file
    base_matrix = args.matrix_file[:-3] if is_gzipped(args.matrix_file) else args.matrix_file
    out_cb_file = f"{base_barcode}.filtered"
    out_mat_file = f"{base_matrix}.filtered"
    
    # Read barcodes
    cb = {}
    with open_file(args.barcode_file) as f:
        for line_num, line in enumerate(f, 1):
            cb[line_num] = line.strip()
    
    # Initialize data structures
    a = {}  # Header lines
    cell_g = {}  # Gene indices
    cell_i = {}  # Cell indices
    cell_n = {}  # UMI counts
    cell_tot = {}  # Total UMIs per cell
    
    # Process matrix file
    with open_file(args.matrix_file) as f:
        for line_num, line in enumerate(f, 1):
            cols = line.strip().split()
            
            if line_num <= n_header_lines:
                a[line_num] = line.strip()
                if line_num == n_header_lines:
                    n_genes = int(cols[0])
            else:
                gene_idx = cols[0]
                cell_idx = int(cols[1])
                umi_count = int(cols[2])
                
                cell_g[line_num] = gene_idx
                cell_i[line_num] = cell_idx
                cell_n[line_num] = umi_count
                
                cell_tot[cell_idx] = cell_tot.get(cell_idx, 0) + umi_count
        
        n_lines = line_num  # Total number of lines read
    
    # Sort cell totals (equivalent to asort in AWK)
    cell_tot_sorted = sorted(cell_tot.values())
    
    # Determine thresholds
    if len(cell_tot_sorted) > 0:
        if exactCells > 0:
            if len(cell_tot) < exactCells:
                n_min = cell_tot_sorted[0]  # Minimum value
            else:
                # Find the threshold that keeps exactly 'exactCells' number of cells
                n_min = cell_tot_sorted[len(cell_tot_sorted) - exactCells]
        else:
            # Calculate based on percentile and ratio
            idx = -int((1 - maxPercentile) * maxCells) + len(cell_tot_sorted)
            idx = min(max(0, idx), len(cell_tot_sorted) - 1)  # Ensure index is within bounds
            n_max = cell_tot_sorted[idx]
            n_min = n_max / maxMinRatio
    else:
        n_min = 0
        n_max = 0
    
    # Filter cells and write to output files
    n_cell = 0
    cell_i_new = {}
    
    with open(out_cb_file, 'w') as out_cb, open(f"{out_cb_file}.counts", 'w') as out_counts:
        for ii in range(1, len(cb) + 1):
            if ii in cell_tot and cell_tot[ii] >= n_min:
                out_cb.write(f"{cb[ii]}\n")
                n_cell += 1
                cell_i_new[ii] = n_cell
                out_counts.write(f"{n_cell} {cell_tot[ii]}\n")
    
    # Print statistics
    if exactCells == 0:
        max_umi = cell_tot_sorted[-1] if cell_tot_sorted else 0
        print(f"maxUMIperCell={max_umi} Robust maxUMIperCel={n_max} minUMIperCell={n_min} Filtered N cells={n_cell}")
    else:
        print(f"total N cells={len(cell_tot)} exactCells={exactCells} minUMIperCell={n_min}")
    
    # Count filtered matrix entries
    n_mat = 0
    for ii in range(n_header_lines + 1, n_lines + 1):
        cell_idx = cell_i[ii]
        if cell_idx in cell_tot and cell_tot[cell_idx] >= n_min:
            n_mat += 1
    
    # Write filtered matrix
    with open(out_mat_file, 'w') as out_mat:
        # Write header lines
        for ii in range(1, n_header_lines):
            out_mat.write(f"{a[ii]}\n")
        
        # Write dimensions line
        out_mat.write(f"{n_genes} {n_cell} {n_mat}\n")
        
        # Write filtered entries
        for ii in range(n_header_lines + 1, n_lines + 1):
            cell_idx = cell_i[ii]
            if cell_idx in cell_tot and cell_tot[cell_idx] >= n_min:
                out_mat.write(f"{cell_g[ii]} {cell_i_new[cell_idx]} {cell_n[ii]}\n")

if __name__ == "__main__":
    main()
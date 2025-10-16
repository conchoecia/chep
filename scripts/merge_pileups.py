#!/usr/bin/env python3
"""
Merge multiple pileup files by position.

Takes a directory of pileup files and combines them by genomic position,
summing depths and concatenating base calls and quality scores.

Loads all data into memory for guaranteed correct sorting.

Handles both .pileup and .pileup.gz files automatically.

Usage:
    python merge_pileups.py <pileup_directory> [--ref-fai <reference.fa.fai>]
    
    --ref-fai: Optional. If provided, output will be sorted in the same order
               as chromosomes appear in the reference .fai file.

Output is written to stdout in pileup format (6 columns):
    chrom, pos, ref, total_depth, bases, quals
"""

import sys
import glob
import os
import gzip
import time
import argparse

def get_memory_usage_mb():
    """Get current memory usage in MB"""
    try:
        import psutil
        process = psutil.Process(os.getpid())
        return process.memory_info().rss / 1024 / 1024
    except ImportError:
        return -1

def read_chrom_order_from_fai(fai_file):
    """
    Read chromosome order from reference .fai file.
    Returns dict mapping chromosome name to its order (0-indexed).
    """
    chrom_order = {}
    with open(fai_file, 'r') as f:
        for idx, line in enumerate(f):
            chrom = line.split('\t')[0]
            chrom_order[chrom] = idx
    print("Read " + str(len(chrom_order)) + " chromosomes from " + fai_file, file=sys.stderr)
    return chrom_order

def get_sort_key(chrom, pos, chrom_order):
    """
    Get sort key for a position.
    If chrom_order is provided, use reference order. Otherwise use lexicographic.
    """
    if chrom_order is not None:
        order = chrom_order.get(chrom, len(chrom_order) + 1000000)
        return (order, pos)
    else:
        return (chrom, pos)

def get_chromosomes_from_files(pileup_files):
    """
    Scan pileup files to get list of all chromosomes.
    Returns set of chromosome names.
    """
    chroms = set()
    for pfile in pileup_files:
        if pfile.endswith('.gz'):
            fh = gzip.open(pfile, 'rt')
        else:
            fh = open(pfile, 'r')
        
        for line in fh:
            chrom = line.split('\t', 1)[0]
            chroms.add(chrom)
        fh.close()
    
    return chroms

def merge_pileups_by_chromosome(pileup_dir, chrom_order=None):
    """
    Merge multiple pileup files by processing one chromosome at a time.
    Memory usage: O(n) where n is positions in the largest chromosome.
    Much more memory-efficient than loading entire genome at once.
    
    Handles both .pileup and .pileup.gz files automatically.
    
    Args:
        pileup_dir: Directory containing .pileup or .pileup.gz files
        chrom_order: Optional dict mapping chromosome names to their order
    
    Yields:
        Merged pileup lines in sorted order
    """
    pileup_files = sorted(glob.glob(pileup_dir + "/*.pileup.gz") +
                         glob.glob(pileup_dir + "/*.pileup"))
    
    pileup_files = sorted(list(set(pileup_files)))
    
    if not pileup_files:
        print("Warning: No pileup files found in " + pileup_dir, file=sys.stderr)
        return
    
    print("Merging " + str(len(pileup_files)) + " pileup files chromosome-by-chromosome...", file=sys.stderr)
    start_time = time.time()
    
    # Get list of all chromosomes from files
    print("Scanning files to identify chromosomes...", file=sys.stderr)
    all_chroms = get_chromosomes_from_files(pileup_files)
    print("Found " + str(len(all_chroms)) + " chromosomes", file=sys.stderr)
    
    # Sort chromosomes by reference order or lexicographically
    if chrom_order is not None:
        sorted_chroms = sorted(all_chroms, key=lambda c: chrom_order.get(c, len(chrom_order) + 1000000))
    else:
        sorted_chroms = sorted(all_chroms)
    
    total_positions = 0
    
    # Process one chromosome at a time
    for chrom_idx, current_chrom in enumerate(sorted_chroms):
        print("Processing chromosome " + str(chrom_idx+1) + "/" + str(len(sorted_chroms)) + ": " + current_chrom + "...", file=sys.stderr)
        
        # Dictionary for this chromosome only: {(pos, ref): [depth, bases, quals]}
        positions = {}
        
        # Read all files and accumulate data for this chromosome
        for pfile in pileup_files:
            if pfile.endswith('.gz'):
                fh = gzip.open(pfile, 'rt')
            else:
                fh = open(pfile, 'r')
            
            for line in fh:
                parts = line.split('\t', 5)
                chrom = parts[0]
                
                # Skip if not the chromosome we're processing
                if chrom != current_chrom:
                    continue
                
                pos = int(parts[1])
                ref = parts[2]
                depth = int(parts[3])
                bases = parts[4] if len(parts) > 4 and parts[4] != '*' else ''
                quals = parts[5].rstrip('\n') if len(parts) > 5 and parts[5] != '*' else ''
                
                key = (pos, ref)
                
                if key in positions:
                    positions[key][0] += depth
                    if bases:
                        positions[key][1] += bases
                    if quals:
                        positions[key][2] += quals
                else:
                    positions[key] = [depth, bases, quals]
            
            fh.close()
        
        # Sort and output this chromosome
        num_positions = len(positions)
        total_positions += num_positions
        
        for key in sorted(positions.keys()):
            pos, ref = key
            depth, bases, quals = positions[key]
            
            bases = bases if bases else '*'
            quals = quals if quals else '*'
            
            yield current_chrom + "\t" + str(pos) + "\t" + ref + "\t" + str(depth) + "\t" + bases + "\t" + quals
        
        mem_mb = get_memory_usage_mb()
        if mem_mb > 0:
            mem_str = str(round(mem_mb, 1)) + " MB"
        else:
            mem_str = "N/A"
        print("  Completed " + current_chrom + ": " + "{:,}".format(num_positions) + " positions, using " + mem_str + " RAM", file=sys.stderr)
        
        # Clear dictionary to free memory
        positions.clear()
    
    total_time = time.time() - start_time
    print("Merge complete: " + "{:,}".format(total_positions) + " total positions in " + str(round(total_time, 1)) + "s", file=sys.stderr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Merge multiple pileup files by genomic position',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    parser.add_argument('pileup_dir', help='Directory containing pileup files')
    parser.add_argument('--ref-fai', help='Reference genome .fai file for chromosome order (optional)')
    
    args = parser.parse_args()
    
    chrom_order = None
    if args.ref_fai:
        chrom_order = read_chrom_order_from_fai(args.ref_fai)

    output_buffer = []
    BUFFER_SIZE = 10000
    
    for line in merge_pileups_by_chromosome(args.pileup_dir, chrom_order):
        output_buffer.append(line)
        if len(output_buffer) >= BUFFER_SIZE:
            sys.stdout.write('\n'.join(output_buffer) + '\n')
            output_buffer.clear()
    
    if output_buffer:
        sys.stdout.write('\n'.join(output_buffer) + '\n')

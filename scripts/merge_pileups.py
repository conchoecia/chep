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

def merge_pileups_in_memory(pileup_dir, chrom_order=None):
    """
    Merge multiple pileup files by loading all data into memory.
    Uses a dictionary keyed by (chrom, pos, ref) to accumulate data.
    Memory usage: O(n) where n is total unique positions across all files.
    
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
    
    print("Merging " + str(len(pileup_files)) + " pileup files by loading into memory...", file=sys.stderr)
    start_time = time.time()
    
    positions = {}
    
    for idx, pfile in enumerate(pileup_files):
        print("Reading file " + str(idx+1) + "/" + str(len(pileup_files)) + ": " + os.path.basename(pfile) + "...", file=sys.stderr)
        
        if pfile.endswith('.gz'):
            fh = gzip.open(pfile, 'rt')
        else:
            fh = open(pfile, 'r')
        
        lines_read = 0
        for line in fh:
            parts = line.split('\t', 5)
            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[2]
            depth = int(parts[3])
            bases = parts[4] if len(parts) > 4 and parts[4] != '*' else ''
            quals = parts[5].rstrip('\n') if len(parts) > 5 and parts[5] != '*' else ''
            
            key = (chrom, pos, ref)
            
            if key in positions:
                positions[key][0] += depth
                if bases:
                    positions[key][1].append(bases)
                if quals:
                    positions[key][2].append(quals)
            else:
                bases_list = [bases] if bases else []
                quals_list = [quals] if quals else []
                positions[key] = [depth, bases_list, quals_list]
            
            lines_read += 1
            if lines_read % 1_000_000 == 0:
                mem_mb = get_memory_usage_mb()
                if mem_mb > 0:
                    mem_str = str(round(mem_mb, 1)) + " MB"
                else:
                    mem_str = "N/A"
                print("  Read " + "{:,}".format(lines_read) + " lines, " + "{:,}".format(len(positions)) + " unique positions, using " + mem_str + " RAM", file=sys.stderr)
        
        fh.close()
        mem_mb = get_memory_usage_mb()
        if mem_mb > 0:
            mem_str = str(round(mem_mb, 1)) + " MB"
        else:
            mem_str = "N/A"
        print("  Completed file " + str(idx+1) + ": " + "{:,}".format(lines_read) + " lines read, " + "{:,}".format(len(positions)) + " total unique positions, using " + mem_str + " RAM", file=sys.stderr)
    
    elapsed = time.time() - start_time
    mem_mb = get_memory_usage_mb()
    if mem_mb > 0:
        mem_str = str(round(mem_mb, 1)) + " MB"
    else:
        mem_str = "N/A"
    print("All files read in " + str(round(elapsed, 1)) + "s. Total unique positions: " + "{:,}".format(len(positions)) + ", using " + mem_str + " RAM", file=sys.stderr)
    
    print("Sorting and outputting " + "{:,}".format(len(positions)) + " positions...", file=sys.stderr)
    output_start = time.time()
    
    def sort_key_func(key):
        chrom, pos, ref = key
        return (get_sort_key(chrom, pos, chrom_order), ref)
    
    last_chrom = None
    positions_in_chrom = 0
    
    for key in sorted(positions.keys(), key=sort_key_func):
        chrom, pos, ref = key
        depth, bases_list, quals_list = positions[key]
        
        bases = ''.join(bases_list) if bases_list else '*'
        quals = ''.join(quals_list) if quals_list else '*'
        
        yield chrom + "\t" + str(pos) + "\t" + ref + "\t" + str(depth) + "\t" + bases + "\t" + quals
        
        if chrom != last_chrom:
            if last_chrom is not None:
                print("Done with " + last_chrom + ": " + "{:,}".format(positions_in_chrom) + " positions", file=sys.stderr)
            last_chrom = chrom
            positions_in_chrom = 0
        positions_in_chrom += 1
    
    if last_chrom is not None:
        print("Done with " + last_chrom + ": " + "{:,}".format(positions_in_chrom) + " positions", file=sys.stderr)
    
    total_time = time.time() - start_time
    output_time = time.time() - output_start
    print("Merge complete: " + "{:,}".format(len(positions)) + " positions in " + str(round(total_time, 1)) + "s (" + str(round(output_time, 1)) + "s for sorting/output)", file=sys.stderr)

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
    
    for line in merge_pileups_in_memory(args.pileup_dir, chrom_order):
        output_buffer.append(line)
        if len(output_buffer) >= BUFFER_SIZE:
            sys.stdout.write('\n'.join(output_buffer) + '\n')
            output_buffer.clear()
    
    if output_buffer:
        sys.stdout.write('\n'.join(output_buffer) + '\n')

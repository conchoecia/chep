#!/usr/bin/env python3
"""
Merge multiple pileup files by position using streaming merge.

Takes a directory of pileup files and combines them by genomic position,
summing depths and concatenating base calls and quality scores.

Uses a memory-efficient streaming merge algorithm that assumes pileup files
are sorted by position (which they are from samtools mpileup).

Usage:
    python merge_pileups.py <pileup_directory>

Output is written to stdout in pileup format (6 columns):
    chrom, pos, ref, total_depth, bases, quals
"""

import sys
import glob
import heapq
import os

def get_memory_usage_mb():
    """Get current memory usage in MB"""
    try:
        import psutil
        process = psutil.Process(os.getpid())
        return process.memory_info().rss / 1024 / 1024
    except ImportError:
        # Fallback if psutil not available
        return -1

class PileupLine:
    """Wrapper for pileup line with comparison for heapq"""
    def __init__(self, line, file_handle, file_idx):
        parts = line.strip().split('\t')
        self.chrom = parts[0]
        self.pos = int(parts[1])
        self.ref = parts[2]
        self.depth = int(parts[3])
        self.bases = parts[4] if len(parts) > 4 else '*'
        self.quals = parts[5] if len(parts) > 5 else '*'
        self.file_handle = file_handle
        self.file_idx = file_idx
    
    def __lt__(self, other):
        # Sort by chromosome name, then position
        return (self.chrom, self.pos) < (other.chrom, other.pos)
    
    def get_key(self):
        return (self.chrom, self.pos, self.ref)

def merge_pileups_streaming(pileup_dir):
    """
    Merge multiple pileup files using a streaming algorithm.
    Memory usage: O(k) where k is number of files, not O(n) where n is positions.
    
    Args:
        pileup_dir: Directory containing .pileup files
    
    Yields:
        Merged pileup lines
    """
    pileup_files = sorted(glob.glob(f"{pileup_dir}/*.pileup"))
    
    if not pileup_files:
        print(f"Warning: No pileup files found in {pileup_dir}", file=sys.stderr)
        return
    
    print(f"Merging {len(pileup_files)} pileup files using streaming merge...", file=sys.stderr)
    
    # Open all files and create initial heap
    file_handles = []
    heap = []
    
    for idx, pfile in enumerate(pileup_files):
        fh = open(pfile, 'r')
        file_handles.append(fh)
        line = fh.readline()
        if line:
            heapq.heappush(heap, PileupLine(line, fh, idx))
    
    # Merge loop
    current_key = None
    current_depth = 0
    current_bases = []
    current_quals = []
    positions_merged = 0
    last_chrom = None
    positions_in_chrom = 0
    
    while heap:
        # Get next line from heap
        pline = heapq.heappop(heap)
        
        # Check if we've moved to a new position
        if current_key is None:
            current_key = pline.get_key()
            last_chrom = pline.chrom
        
        # Check if we've moved to a new chromosome
        if pline.chrom != last_chrom:
            mem_mb = get_memory_usage_mb()
            mem_str = f"{mem_mb:.1f} MB" if mem_mb > 0 else "N/A"
            print(f"Done with {last_chrom}: {positions_in_chrom:,} positions merged, currently using {mem_str} RAM", 
                  file=sys.stderr)
            last_chrom = pline.chrom
            positions_in_chrom = 0
        
        if pline.get_key() == current_key:
            # Same position, accumulate data
            current_depth += pline.depth
            if pline.bases != '*':
                current_bases.append(pline.bases)
            if pline.quals != '*':
                current_quals.append(pline.quals)
        else:
            # New position, output accumulated data
            chrom, pos, ref = current_key
            bases = ''.join(current_bases) if current_bases else '*'
            quals = ''.join(current_quals) if current_quals else '*'
            yield f"{chrom}\t{pos}\t{ref}\t{current_depth}\t{bases}\t{quals}"
            positions_merged += 1
            positions_in_chrom += 1
            
            # Progress reporting every 1M positions
            if positions_merged % 1_000_000 == 0:
                mem_mb = get_memory_usage_mb()
                mem_str = f"{mem_mb:.1f} MB" if mem_mb > 0 else "N/A"
                print(f"Progress: {positions_merged:,} positions merged on {chrom}, currently using {mem_str} RAM", 
                      file=sys.stderr)
            
            # Reset for new position
            current_key = pline.get_key()
            current_depth = pline.depth
            current_bases = [pline.bases] if pline.bases != '*' else []
            current_quals = [pline.quals] if pline.quals != '*' else []
        
        # Read next line from this file
        line = pline.file_handle.readline()
        if line:
            heapq.heappush(heap, PileupLine(line, pline.file_handle, pline.file_idx))
    
    # Output last position
    if current_key is not None:
        chrom, pos, ref = current_key
        bases = ''.join(current_bases) if current_bases else '*'
        quals = ''.join(current_quals) if current_quals else '*'
        yield f"{chrom}\t{pos}\t{ref}\t{current_depth}\t{bases}\t{quals}"
        positions_merged += 1
        positions_in_chrom += 1
        
        # Final chromosome report
        mem_mb = get_memory_usage_mb()
        mem_str = f"{mem_mb:.1f} MB" if mem_mb > 0 else "N/A"
        print(f"Done with {last_chrom}: {positions_in_chrom:,} positions merged, currently using {mem_str} RAM", 
              file=sys.stderr)

    # Close all files
    for fh in file_handles:
        fh.close()

    print(f"Merged {positions_merged} positions", file=sys.stderr)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__, file=sys.stderr)
        sys.exit(1)

    pileup_dir = sys.argv[1]

    for line in merge_pileups_streaming(pileup_dir):
        print(line)

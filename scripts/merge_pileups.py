#!/usr/bin/env python3
"""
Merge multiple pileup files by position.

Takes a directory of pileup files and combines them by genomic position,
summing depths and concatenating base calls and quality scores.

Usage:
    python merge_pileups.py <pileup_directory>

Output is written to stdout in pileup format (6 columns):
    chrom, pos, ref, total_depth, bases, quals
"""

import sys
import glob
from collections import defaultdict

def merge_pileups(pileup_dir):
    """
    Merge multiple pileup files from a directory.
    
    Args:
        pileup_dir: Directory containing .pileup files
    
    Yields:
        Tuples of (chrom, pos, ref, depth, bases, quals) sorted by position
    """
    # Read all pileup files
    pileups = defaultdict(lambda: {'depth': 0, 'bases': '', 'quals': ''})
    
    pileup_files = glob.glob(f"{pileup_dir}/*.pileup")
    
    if not pileup_files:
        print(f"Warning: No pileup files found in {pileup_dir}", file=sys.stderr)
        return
    
    print(f"Merging {len(pileup_files)} pileup files...", file=sys.stderr)
    
    for pfile in pileup_files:
        with open(pfile) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    chrom, pos, ref, depth, bases, quals = parts[:6]
                    key = (chrom, int(pos), ref)
                    pileups[key]['depth'] += int(depth)
                    if bases != '*':
                        pileups[key]['bases'] += bases
                    if quals != '*':
                        pileups[key]['quals'] += quals
    
    print(f"Merged {len(pileups)} unique positions", file=sys.stderr)
    
    # Output merged pileups sorted by (chrom, position)
    for (chrom, pos, ref), data in sorted(pileups.items()):
        depth = data['depth']
        bases = data['bases'] if data['bases'] else '*'
        quals = data['quals'] if data['quals'] else '*'
        yield f"{chrom}\t{pos}\t{ref}\t{depth}\t{bases}\t{quals}"

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__, file=sys.stderr)
        sys.exit(1)
    
    pileup_dir = sys.argv[1]
    
    for line in merge_pileups(pileup_dir):
        print(line)

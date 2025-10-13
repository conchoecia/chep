#!/usr/bin/env python3
"""
Merge multiple chep_3D.txt files by summing frequencies for matching (depth, count) pairs.
Usage: merge_chep3d.py file1.txt file2.txt ... > merged.txt
       Or with a file list: merge_chep3d.py $(cat file_list.txt) > merged.txt
"""
import sys
from collections import defaultdict

def merge_chep3d_files(file_paths):
    """Merge multiple chep_3D.txt files by summing frequencies."""
    # Use nested defaultdict to store depth -> count -> frequency
    merged = defaultdict(lambda: defaultdict(int))
    
    for file_path in file_paths:
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split('\t')
                    if len(parts) != 3:
                        continue
                    depth = int(parts[0])
                    count = int(parts[1])
                    freq = int(parts[2])
                    merged[depth][count] += freq
        except Exception as e:
            print(f"Error processing {file_path}: {e}", file=sys.stderr)
            continue
    
    # Output in sorted order (by depth, then by count)
    for depth in sorted(merged.keys()):
        for count in sorted(merged[depth].keys()):
            print(f"{depth}\t{count}\t{merged[depth][count]}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: merge_chep3d.py file1.txt file2.txt ... > merged.txt", file=sys.stderr)
        sys.exit(1)
    
    merge_chep3d_files(sys.argv[1:])

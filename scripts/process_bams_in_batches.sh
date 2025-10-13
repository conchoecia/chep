#!/bin/bash
# Process multiple BAM files in batches using GNU parallel and merge results
# Usage: process_bams_in_batches.sh <reference.fa> <bam_list.txt> <output.txt> [batch_size] [threads]

set -e

REF=$1
BAM_LIST=$2
OUTPUT=$3
BATCH_SIZE=${4:-50}    # Default to 50 BAMs per batch
THREADS=${5:-8}        # Default to 8 parallel jobs

if [ $# -lt 3 ]; then
    echo "Usage: $0 <reference.fa> <bam_list.txt> <output.txt> [batch_size] [threads]"
    echo "  reference.fa  - Reference genome FASTA file"
    echo "  bam_list.txt  - File with one BAM path per line"
    echo "  output.txt    - Output merged chep_3D.txt file"
    echo "  batch_size    - Number of BAMs per batch (default: 50)"
    echo "  threads       - Number of parallel jobs (default: 8)"
    exit 1
fi

# Check if GNU parallel is available
if ! command -v parallel &> /dev/null; then
    echo "Error: GNU parallel is not installed or not in PATH"
    echo "Install with: conda install -c conda-forge parallel"
    echo "Or: brew install parallel (on macOS)"
    exit 1
fi

# Get the directory where this script lives (to find merge_chep3d.py)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Create temporary directory for batch outputs
TMPDIR=$(mktemp -d)
echo "Processing BAMs in batches of ${BATCH_SIZE} using ${THREADS} parallel jobs..."
echo "Temporary directory: ${TMPDIR}"

# Split BAM list into batches
split -l ${BATCH_SIZE} ${BAM_LIST} ${TMPDIR}/batch_

# Count batches
BATCH_COUNT=$(ls ${TMPDIR}/batch_* | wc -l)
echo "Created ${BATCH_COUNT} batches"

# Function to process a single batch (will be called by parallel)
process_batch() {
    local batch=$1
    local ref=$2
    local tmpdir=$3
    local batch_name=$(basename ${batch})
    
    samtools mpileup -B -Q 0 -q 0 -A -R -f ${ref} -b ${batch} | \
        chep_pileup_to_array > ${tmpdir}/output_${batch_name}.txt
    
    echo "Completed: ${batch_name}"
}

# Export the function and variables so parallel can use them
export -f process_batch
export REF TMPDIR

# Process all batches in parallel
echo "Processing batches in parallel..."
ls ${TMPDIR}/batch_* | parallel -j ${THREADS} process_batch {} ${REF} ${TMPDIR}

# Merge all batch outputs
echo "Merging ${BATCH_COUNT} batch outputs..."
python3 ${SCRIPT_DIR}/merge_chep3d.py ${TMPDIR}/output_batch_*.txt > ${OUTPUT}

# Cleanup
echo "Cleaning up temporary files..."
rm -rf ${TMPDIR}

echo "Done! Output written to: ${OUTPUT}"

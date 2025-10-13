#!/bin/bash
#SBATCH --job-name=chep_snakemake
#SBATCH --output=logs/chep_snakemake_%j.log
#SBATCH --error=logs/chep_snakemake_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# SLURM submission script for CHEP Snakemake workflow
# This is the main controller job that submits individual BAM processing jobs
#
# Usage:
#   sbatch submit_chep_snakemake.sh <reference.fa> <bam_list.txt> <output.txt> [max_jobs]
#
# Example:
#   sbatch submit_chep_snakemake.sh ref.fa bams.txt output.txt 100

set -e

# Load required modules (adjust for your cluster)
# module load samtools
# module load python3
# module load snakemake

# Check arguments
if [ $# -lt 3 ]; then
    echo "Usage: sbatch $0 <reference.fa> <bam_list.txt> <output.txt> [max_jobs]"
    exit 1
fi

REF=$1
BAM_LIST=$2
OUTPUT=$3
MAX_JOBS=${4:-100}

# Create logs directory if it doesn't exist
mkdir -p logs

# Get the directory where this script lives
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SNAKEFILE="${SCRIPT_DIR}/../snakemake_scripts/Snakefile_chep_multi_bam.snakemake"

echo "Starting CHEP Snakemake workflow"
echo "Reference: ${REF}"
echo "BAM list: ${BAM_LIST}"
echo "Output: ${OUTPUT}"
echo "Max parallel jobs: ${MAX_JOBS}"
echo "Snakefile: ${SNAKEFILE}"

# Run Snakemake with SLURM executor
snakemake -s ${SNAKEFILE} \
    --config ref=${REF} bam_list=${BAM_LIST} output=${OUTPUT} \
    --executor slurm \
    --jobs ${MAX_JOBS} \
    --default-resources slurm_account=default slurm_partition=default \
    --latency-wait 60 \
    --rerun-incomplete \
    --printshellcmds

echo "Workflow complete! Output: ${OUTPUT}"

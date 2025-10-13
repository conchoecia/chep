#!/bin/bash
# Wrapper script to run CHEP multi-BAM Snakemake workflow
# Usage: run_chep_snakemake.sh <reference.fa> <bam_list.txt> <output.txt> [max_jobs]

###### SET SLURM OPTIONS IF NEEDED ######
#SBATCH --job-name=chp
#SBATCH --output=HC2_%j.out
#SBATCH --error=HC2_%j.err
#SBATCH --partition=tttt,owners
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00

module load biology
module load samtools
### LOCAL MODULE LOADING FOR SAMTOOLS ###

set -e

REF=$1
BAM_LIST=$2
OUTPUT=$3
MAX_JOBS=${4:-100}

if [ $# -lt 3 ]; then
    echo "Usage: $0 <reference.fa> <bam_list.txt> <output.txt> [max_jobs]"
    echo ""
    echo "Arguments:"
    echo "  reference.fa  - Reference genome FASTA file"
    echo "  bam_list.txt  - File with one BAM path per line"
    echo "  output.txt    - Output merged chep_3D.txt file"
    echo "  max_jobs      - Maximum parallel jobs (default: 100)"
    echo ""
    echo "Example:"
    echo "  $0 ref.fa bams.txt output.txt 50"
    echo ""
    echo "For SLURM cluster, this will submit jobs via the Snakemake SLURM executor."
    echo "For local execution, use: $0 ref.fa bams.txt output.txt 8"
    exit 1
fi

# Get the directory where this script lives
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SNAKEFILE="${SCRIPT_DIR}/../snakemake_scripts/Snakefile_chep_multi_bam.snakemake"

# Create logs directory if it doesn't exist
mkdir -p logs

# Check if we're on a SLURM cluster
if command -v sbatch &> /dev/null; then
    echo "SLURM detected - using Snakemake SLURM executor"
    snakemake -s ${SNAKEFILE} \
        --config ref=${REF} bam_list=${BAM_LIST} output=${OUTPUT} \
        --executor slurm \
        --jobs ${MAX_JOBS} \
        --default-resources slurm_account=default slurm_partition=default \
        --latency-wait 60 \
        --rerun-incomplete
else
    echo "Running locally with ${MAX_JOBS} cores"
    snakemake -s ${SNAKEFILE} \
        --config ref=${REF} bam_list=${BAM_LIST} output=${OUTPUT} \
        --cores ${MAX_JOBS} \
        --rerun-incomplete
fi

echo "Workflow complete!"

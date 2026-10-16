#!/bin/bash
#SBATCH --job-name=find.sft.power
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=70G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=72:50:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err

set -euo pipefail

# ============================================================
# Identify WGCNA soft-thresholding power within sex
# Reference: Age group 1
# ============================================================

# Configure HPC directories
output_path="..."
scratch="..."
RSCRIPT="..."

# Create temporary working directory
mkdir -p "${scratch}/workdir" "$output_path"

workdir=$(mktemp -d "${scratch}/workdir/tmpXXXXXX")
trap 'rm -rf "$workdir"' EXIT

# Load R
module load r/4.4.2-gfbf-2024a

# Run soft-threshold analysis
cd "$workdir"

Rscript --slave "${RSCRIPT}/01_find_soft_power_by_sex.R"

# Copy generated figures, if any
shopt -s nullglob
png_files=( *.png )

if (( ${#png_files[@]} > 0 )); then
    rsync -av "${png_files[@]}" "$output_path/"
fi

echo "Soft-threshold analysis completed."

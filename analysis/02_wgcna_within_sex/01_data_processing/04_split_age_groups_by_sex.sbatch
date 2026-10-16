#!/bin/bash
#SBATCH --job-name=split.chr
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=30G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err

set -euo pipefail

# ============================================================
# Split pre-corrected DNAm data into eight age groups
# within sex, processing one chromosome per SLURM job.
#
# Example submission:
# for chr in {1..22}; do
#     sbatch --export=ALL,chr=$chr 04_split_age_groups_by_sex.sbatch
# done
# ============================================================

# Configure HPC directories
dnam_inpath="..."
scratch="..."
RSCRIPT="..."

# Check chromosome identifier
: "${chr:?Please supply chr using sbatch --export=ALL,chr=N}"

echo "Present working directory: $(pwd)"

# Create temporary working directory
workdir=$(mktemp -d "${scratch}/tmpXXXXXX")
trap 'rm -rf "$workdir"' EXIT

echo "Temporary directory: $workdir"

# Substitute chromosome placeholder in R script
sed "s/@chr@/${chr}/g" \
    "${RSCRIPT}/04_split_age_groups_by_sex.R" \
    > "${workdir}/run_script.R"

# Load R
module load r/4.4.2-gfbf-2024a

# Run analysis
cd "$workdir"

Rscript --slave run_script.R

echo "Completed chromosome: $chr"

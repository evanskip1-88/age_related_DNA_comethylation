#!/bin/bash
#SBATCH --job-name=comb.age.grps
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=70G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=00:50:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err
#SBATCH --array=1-8

set -euo pipefail

# ============================================================
# Combine chromosome-level DNAm data within age groups
# Each SLURM array task processes one age group (Age1–Age8).
#
# Submission:
# sbatch 05_combine_age_group_chromosomes.sbatch
# ============================================================

# Configure HPC directories
output_path="..."
scratch="..."
RSCRIPT="..."

# Age group assigned by SLURM array
AgeGrp="${SLURM_ARRAY_TASK_ID}"

echo "Processing age group: $AgeGrp"
echo "Present working directory: $(pwd)"

# Create temporary working directory
workdir=$(mktemp -d "${scratch}/workdir/tmpXXXXXX")
trap 'rm -rf "$workdir"' EXIT

echo "Temporary directory: $workdir"

mkdir -p "$output_path"

# Substitute age-group placeholder in R script
sed "s/@AgeGrp@/${AgeGrp}/g" \
    "${RSCRIPT}/05_combine_age_group_chromosomes.R" \
    > "${workdir}/run_script.R"

# Load R
module load r/4.4.2-gfbf-2024a

# Run analysis
cd "$workdir"

Rscript --slave run_script.R

echo "Completed age group: $AgeGrp"

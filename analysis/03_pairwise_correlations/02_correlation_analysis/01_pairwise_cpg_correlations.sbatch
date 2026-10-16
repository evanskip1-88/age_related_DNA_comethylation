#!/bin/bash
#SBATCH --job-name=GS.pairwise.cor
#SBATCH --account=<HPC_ACCOUNT>
#SBATCH --partition=<HPC_PARTITION>
#SBATCH --mem=400G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --array=2-8
#SBATCH --output=log.out/%x.%A_%a.out
#SBATCH --error=log.out/%x.%A_%a.err

set -euo pipefail

# ============================================================
# Pairwise CpG correlations across age groups
#
# Each array task processes:
#   One chromosome
#   One reference age group (Age2-Age8)
#
# CpG pairs with |r| > 0.8 in the reference age group
# are evaluated across all eight age groups.
# ============================================================

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

dnam_path="..."
cov_inpath="..."
output_path="..."
script_path="/analysis/03_pairwise_correlations/02_correlation_analysis"

chr="${chr:?ERROR: Supply chr=1..22}"

ageGrp="${SLURM_ARRAY_TASK_ID}"

# ------------------------------------------------------------
# Validate inputs
# ------------------------------------------------------------

if ! [[ "$chr" =~ ^([1-9]|1[0-9]|2[0-2])$ ]]; then
    echo "ERROR: Invalid chromosome: $chr"
    exit 1
fi

mkdir -p "$output_path"

# ------------------------------------------------------------
# Load R
# ------------------------------------------------------------

module load r/4.4.2-gfbf-2024a

echo "Analysis started: $(date)"
echo "Chromosome: $chr"
echo "Reference age group: Age${ageGrp}"

# ------------------------------------------------------------
# Run correlation analysis
# ------------------------------------------------------------

Rscript \
    "${script_path}/01_pairwise_cpg_correlations.R" \
    "$chr" \
    "$ageGrp" \
    "$dnam_path" \
    "$cov_inpath" \
    "$output_path"

echo "Analysis completed: $(date)"

#!/bin/bash
#SBATCH --job-name=covar
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=30G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err

set -euo pipefail

# =============================================================================
# Pairwise CpG covariance analysis
#
# Runs the Age1-reference covariance ANOVA analysis across chromosomes 1–22.
#
# Usage:
#   mkdir -p logs
#   sbatch 07_calculate_pairwise_covariance.sbatch
#
# Edit the paths below before submission.
# =============================================================================

SCRIPT_DIR="/analysis/03_pairwise_correlations/02_correlation_analysis"
SCRATCH_DIR="..."

# Load R
module load r/4.4.2-gfbf-2024a

# Create temporary working directory
mkdir -p "$SCRATCH_DIR"

WORKDIR=$(mktemp -d "${SCRATCH_DIR}/covar_XXXXXX")

# Remove temporary files on exit
trap 'rm -rf "$WORKDIR"' EXIT

echo "Job ID: ${SLURM_JOB_ID:-NA}"
echo "Working directory: $WORKDIR"
echo "Started: $(date)"

# Run the R analysis
cd "$WORKDIR"

Rscript "${SCRIPT_DIR}/07_calculate_pairwise_covariance.R"

echo "Completed: $(date)"

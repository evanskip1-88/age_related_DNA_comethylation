#!/bin/bash

#SBATCH --job-name=find_soft_power
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=70G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=96:50:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err


# ==============================================================================
# Evaluate WGCNA soft-thresholding powers
#
# Runs 01_find_soft_power.R using the youngest age group (Age group 1) as the
# reference dataset for soft-threshold selection.
# ==============================================================================


# ------------------------------------------------------------------------------
# User-defined paths
# ------------------------------------------------------------------------------

output_path="..."
scratch="..."
R_SCRIPT="..."


# Exit immediately if a command fails
set -e


# ------------------------------------------------------------------------------
# Create temporary working directory
# ------------------------------------------------------------------------------

echo "Present working directory: $(pwd)"

workdir=$(mktemp -d "${scratch}/workdir/tmpXXXXXX")

echo "Temporary working directory: ${workdir}"

mkdir -p "${output_path}"

trap 'rm -rf "${workdir}"' EXIT

cd "${workdir}" || exit 1


# ------------------------------------------------------------------------------
# Load R
# ------------------------------------------------------------------------------

module load r


# ------------------------------------------------------------------------------
# Run soft-threshold analysis
# ------------------------------------------------------------------------------

echo "Starting soft-threshold analysis..."

Rscript --slave "${R_SCRIPT}/01_find_soft_power.R"


# ------------------------------------------------------------------------------
# Copy output files
# ------------------------------------------------------------------------------

rsync -av ./*.png "${output_path}/"

echo "Soft-threshold analysis complete."

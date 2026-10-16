#!/bin/bash

#SBATCH --job-name=construct_reference_network
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=150G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err


# ==============================================================================
# Construct the reference DNA co-methylation network
#
# Runs WGCNA on residualised DNAm data from the youngest age group
# (Age group 1) using 10 CPU threads.
#
# Network parameters are defined in:
#   02_construct_reference_network.R
# ==============================================================================


# ------------------------------------------------------------------------------
# User-defined paths
# ------------------------------------------------------------------------------

scratch="..."
R_SCRIPT="analysis/01_wgcna_combined_sex/02_network_analysis/"
output_path="..."


# Exit immediately if a command fails
set -e


# ------------------------------------------------------------------------------
# Start timer
# ------------------------------------------------------------------------------

start_time=$(date +%s)

echo "Start time: $(date)"


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
# Construct reference network
# ------------------------------------------------------------------------------

echo "Starting reference network construction..."

Rscript "${R_SCRIPT}/02_construct_reference_network.R"


# ------------------------------------------------------------------------------
# Finish
# ------------------------------------------------------------------------------

end_time=$(date +%s)

execution_time=$((end_time - start_time))
execution_time_min=$((execution_time / 60))

echo "End time: $(date)"
echo "Total execution time: ${execution_time_min} minutes"
echo "Reference network construction complete."

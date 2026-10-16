#!/bin/bash
#SBATCH --job-name=wgcna.preservation
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err

# -------------------------------------------------------------------------
# Run whole-genome WGCNA module preservation analysis
#
# Reference:
#   Age group 1 (<36 years)
#
# Test datasets:
#   Age groups 2-8
#
# The corresponding R script evaluates preservation of modules defined
# in the reference age group across the remaining age groups.
# -------------------------------------------------------------------------

set -euo pipefail

# Paths
scratch="..."
RSCRIPT="analysis/01_wgcna_combined_sex/02_network_analysis"
output_path="results/01_wgcna_combined_sex"

mkdir -p "${output_path}"
mkdir -p log.out

# Record start time
start_time=$(date +%s)
echo "Analysis started: $(date)"

# Create temporary working directory
workdir=$(mktemp -d "${scratch}/workdir/tmpXXXXXX")
echo "Temporary working directory: ${workdir}"

# Remove temporary directory on exit
trap 'rm -rf "${workdir}"' EXIT

cd "${workdir}"

# Load R
module load r/4.4.2-gfbf-2024a

# Run module-preservation analysis
Rscript "${RSCRIPT}/03_module_preservation.R"

# Record runtime
end_time=$(date +%s)
execution_time=$((end_time - start_time))
execution_time_min=$((execution_time / 60))

echo "Analysis ended: $(date)"
echo "Total execution time: ${execution_time_min} minutes"

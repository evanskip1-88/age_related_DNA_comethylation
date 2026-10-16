#!/bin/bash
#SBATCH --job-name=wgcna.eigengenes
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=05:00:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err

# -------------------------------------------------------------------------
# Calculate WGCNA module eigengenes and average module profiles
#
# Modules identified in the reference age group (Age1) are used to calculate
# module eigengenes and average module profiles across age groups.
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

echo "Present working directory: $(pwd)"

# Create temporary working directory
workdir=$(mktemp -d "${scratch}/workdir/tmpXXXXXX")
echo "Temporary working directory: ${workdir}"

cd "${workdir}"

# Load R
module load r

# Run module eigengene analysis
Rscript "${RSCRIPT}/04_module_eigengenes.R"

# Clean up temporary files
cd ..
rm -rf "${workdir}"

# Record runtime
end_time=$(date +%s)
execution_time=$((end_time - start_time))
execution_time_min=$((execution_time / 60))

echo "Analysis ended: $(date)"
echo "Total execution time: ${execution_time_min} minutes"

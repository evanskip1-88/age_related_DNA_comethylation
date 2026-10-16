#!/bin/bash
#SBATCH --job-name=pres.wgcna.by.sex
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err

set -euo pipefail

# ============================================================
# Assess WGCNA module preservation across age groups
# within and between sexes
#
# Reference: Age1
# Test groups: Age2–Age8
# Observed preservation statistics (nPermutations = 0)
# ============================================================

# Record start time
start_time=$(date +%s)
echo "Analysis started: $(date)"

# Configure HPC directories
scratch="..."
RSCRIPT="/analysis/02_wgcna_within_sex/02_network_analysis"
output_path="/results/01_wgcna_combined_sex"

# Create temporary working directory
mkdir -p "${scratch}/workdir" "$output_path"

workdir=$(mktemp -d "${scratch}/workdir/tmpXXXXXX")
trap 'rm -rf "$workdir"' EXIT

echo "Temporary directory: $workdir"

# Load R
module load r/4.4.2-gfbf-2024a

# Run module preservation analysis
cd "$workdir"

Rscript --slave "${RSCRIPT}/03_module_preservation_by_sex.R"

# Record execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))

echo "Analysis completed: $(date)"
echo "Total execution time: $((execution_time / 60)) minutes"

#!/bin/bash
#SBATCH --job-name=wgcna.cluster.by.sex
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err

set -euo pipefail

# ============================================================
# Construct sex-specific WGCNA reference network (Age1)
# Soft-threshold power = 6, empirically selected from females
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

# Load R
module load r/4.4.2-gfbf-2024a

# Run WGCNA reference network construction
cd "$workdir"

Rscript --slave "${RSCRIPT}/02_construct_reference_network_by_sex.R"

# Record execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))

echo "Analysis completed: $(date)"
echo "Total execution time: $((execution_time / 60)) minutes"

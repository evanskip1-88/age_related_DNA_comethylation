#!/bin/bash
#SBATCH --job-name=rGREAT.wg.by.sex
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=30G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=5:00:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err

set -euo pipefail

# ============================================================
# GREAT enrichment of sex-biased co-methylation modules
#
# Reference network: Male Age1
# ============================================================

# Record start time
start_time=$(date +%s)
echo "Analysis started: $(date)"

# Configure HPC directories
scratch="..."
RSCRIPT="..."

# Create temporary working directory
mkdir -p "${scratch}/workdir"

workdir=$(mktemp -d "${scratch}/workdir/tmpXXXXXX")
trap 'rm -rf "$workdir"' EXIT

echo "Temporary directory: $workdir"

# Load R
module load r/4.4.2-gfbf-2024a

# Run GREAT enrichment
cd "$workdir"

Rscript --slave "${RSCRIPT}/05_rGREAT_sex_biased_modules.R"

# Record execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))

echo "Analysis completed: $(date)"
echo "Total execution time: $((execution_time / 60)) minutes"

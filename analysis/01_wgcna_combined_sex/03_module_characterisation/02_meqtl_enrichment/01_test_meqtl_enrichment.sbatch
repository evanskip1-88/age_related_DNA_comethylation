#!/bin/bash
#SBATCH --job-name=meqtl.enrichment
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err

# -------------------------------------------------------------------------
# Test enrichment of age-related co-methylation CpG sets among
# cis- and trans-meQTL-associated CpGs from GoDMC
#
# The analysis is restricted to CpGs represented in both the Generation
# Scotland EPIC dataset and the Illumina HumanMethylation450 array.
# -------------------------------------------------------------------------

set -euo pipefail

# Paths
scratch="..."
RSCRIPT="analysis/01_wgcna_combined_sex/03_module_characterisation/02_meqtl_enrichment"

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
module load r/4.4.2-gfbf-2024a

# Run meQTL enrichment analysis
Rscript "${RSCRIPT}/01_test_meqtl_enrichment.R"

# Clean temporary directory
cd ..
rm -rf "${workdir}"

# Record execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))
execution_time_min=$((execution_time / 60))

echo "Analysis ended: $(date)"
echo "Total execution time: ${execution_time_min} minutes"

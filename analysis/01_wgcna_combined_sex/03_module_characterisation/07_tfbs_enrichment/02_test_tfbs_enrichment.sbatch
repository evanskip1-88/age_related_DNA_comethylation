#!/bin/bash
#SBATCH --job-name=tfbsNoPerm.enrich
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err

#==============================================================================
# ENCODE TFBS enrichment analysis
#
# Runs one-sided hypergeometric enrichment tests for age-stable and
# age-variable CpG sets without permutation.
#==============================================================================

set -euo pipefail

#------------------------------------------------------------------------------
# Paths
#------------------------------------------------------------------------------

scratch="..."
RSCRIPT="/analysis/01_wgcna_combined_sex/03_module_characterisation/07_tfbs_enrichment"

#------------------------------------------------------------------------------
# Set up temporary working directory
#------------------------------------------------------------------------------

start_time=$(date +%s)

echo "Analysis started: $(date)"
echo "Submission directory: $(pwd)"

mkdir -p "${scratch}/workdir"

workdir=$(mktemp -d "${scratch}/workdir/tmpXXXXXX")

echo "Temporary working directory: ${workdir}"

# Automatically remove temporary files when the job exits
trap 'rm -rf -- "${workdir}"' EXIT

cd "${workdir}"

#------------------------------------------------------------------------------
# Load R and run enrichment
#------------------------------------------------------------------------------

module load r/4.4.2-gfbf-2024a

Rscript "${RSCRIPT}/02_test_tfbs_enrichment.R"

#------------------------------------------------------------------------------
# Report execution time
#------------------------------------------------------------------------------

end_time=$(date +%s)
execution_time=$((end_time - start_time))

echo "Analysis completed: $(date)"
echo "Total execution time: $((execution_time / 60)) minutes"

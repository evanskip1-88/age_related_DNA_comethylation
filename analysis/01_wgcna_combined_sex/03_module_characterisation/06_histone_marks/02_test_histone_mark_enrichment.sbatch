#!/bin/bash
#SBATCH --job-name=histone.enrich
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err

#==============================================================================
# Histone-mark enrichment analysis
#
# Runs enrichment of stable and changing WGCNA and pairwise CpG sets
# across published histone-mark regions.
#==============================================================================

set -euo pipefail

#------------------------------------------------------------------------------
# Paths
#------------------------------------------------------------------------------

scratch="..."
RSCRIPT="..."

#------------------------------------------------------------------------------
# Set up temporary working directory
#------------------------------------------------------------------------------

start_time=$(date +%s)

echo "Analysis started: $(date)"
echo "Submission directory: $(pwd)"

mkdir -p "${scratch}/workdir"

workdir=$(mktemp -d "${scratch}/workdir/tmpXXXXXX")

echo "Temporary directory: ${workdir}"

# Clean up temporary files when the job exits
trap 'rm -rf -- "${workdir}"' EXIT

cd "${workdir}"

#------------------------------------------------------------------------------
# Run enrichment analysis
#------------------------------------------------------------------------------

module load r

Rscript "${RSCRIPT}/01_test_histone_mark_enrichment.R"

#------------------------------------------------------------------------------
# Report execution time
#------------------------------------------------------------------------------

end_time=$(date +%s)
execution_time=$((end_time - start_time))

echo "Analysis completed: $(date)"
echo "Total execution time: $((execution_time / 60)) minutes"

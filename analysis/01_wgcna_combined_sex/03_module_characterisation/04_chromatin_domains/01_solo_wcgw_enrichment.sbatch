#!/bin/bash
#SBATCH --job-name=PMD.enrich
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
# Enrichment of co-methylated CpG sets in PMDs and solo-WCGW regions
#
# Runs the human hg19 PMD/solo-WCGW enrichment analysis.
#==============================================================================

set -euo pipefail

#------------------------------------------------------------------------------
# Paths
#------------------------------------------------------------------------------

scratch="..."
RSCRIPT="/analysis/01_wgcna_combined_sex/03_module_characterisation/04_chromatin_domains"

#------------------------------------------------------------------------------
# Set up temporary working directory
#------------------------------------------------------------------------------

start_time=$(date +%s)

echo "Analysis started: $(date)"
echo "Present working directory: $(pwd)"

mkdir -p "${scratch}/workdir"

workdir=$(mktemp -d "${scratch}/workdir/tmpXXXXXX")

echo "Temporary directory: ${workdir}"

# Remove temporary files on exit
trap 'rm -rf -- "${workdir}"' EXIT

cd "${workdir}"

#------------------------------------------------------------------------------
# Run enrichment analysis
#------------------------------------------------------------------------------

module load r

Rscript "${RSCRIPT}/01_solo_wcgw_enrichment.R"

#------------------------------------------------------------------------------
# Report execution time
#------------------------------------------------------------------------------

end_time=$(date +%s)
execution_time=$((end_time - start_time))
execution_time_min=$((execution_time / 60))

echo "Analysis completed: $(date)"
echo "Total execution time: ${execution_time_min} minutes"

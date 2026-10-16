#!/bin/bash
#SBATCH --job-name=sub.raw.dnam
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err
#SBATCH --array=1-4

#==============================================================================
# Subset raw DNAm values for selected co-methylation CpG sets
#
# Each array task selects one CpG set from the combined pairwise/WGCNA
# CpG-set list and extracts the corresponding DNAm values.
#==============================================================================

set -euo pipefail

cpg_hit_num="${SLURM_ARRAY_TASK_ID}"

echo "Processing CpG-set index: ${cpg_hit_num}"

#------------------------------------------------------------------------------
# Paths
#------------------------------------------------------------------------------

output_path="results/01_wgcna_combined_sex"
scratch="..."
RSCRIPT="analysis/01_wgcna_combined_sex/03_module_characterisation/03_smoking_variance"

mkdir -p "${output_path}"
mkdir -p log.out

#------------------------------------------------------------------------------
# Temporary working directory
#------------------------------------------------------------------------------

start_time=$(date +%s)

echo "Analysis started: $(date)"

workdir=$(mktemp -d "${scratch}/workdir/tmpXXXXXX")

echo "Temporary working directory: ${workdir}"

cd "${workdir}"

# Insert CpG-set index into R script
sed \
  "s/@cpgHitSet@/${cpg_hit_num}/g" \
  "${RSCRIPT}/01_subset_raw_beta_cpg_sets.R" \
  > run_script.R

#------------------------------------------------------------------------------
# Run analysis
#------------------------------------------------------------------------------

module load r

Rscript run_script.R

#------------------------------------------------------------------------------
# Clean up
#------------------------------------------------------------------------------

cd /
rm -rf "${workdir}"

end_time=$(date +%s)

execution_time=$((end_time - start_time))
execution_time_min=$((execution_time / 60))

echo "Analysis ended: $(date)"
echo "Total execution time: ${execution_time_min} minutes"

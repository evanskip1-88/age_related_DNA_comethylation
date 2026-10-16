#!/bin/bash
#SBATCH --job-name=smoking.oreml
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
# Estimate smoking-related variance explained by co-methylated CpG sets
#
# For each of the four CpG sets (aSCPs, aVCPs, aSCMs and aVCMs):
#   1. Convert beta-value data to OSCA binary format
#   2. Update probe information
#   3. Construct an omics relationship matrix (ORM)
#   4. Estimate variance explained by the CpG set using OREML
#
# Smoking phenotypes:
#   - Ever smoking
#   - Smoking pack-years
#
# OREML models are adjusted for age, sex and estimated blood-cell composition.
#==============================================================================

set -euo pipefail

#------------------------------------------------------------------------------
# Paths
#------------------------------------------------------------------------------

gs_cov_path="..." # Indivdial level data
gs_opi="data/gsDNAm"
output="results/01_wgcna_combined_sex"
gs_dnam="..." # Indivdial level data
scratch="..."

mkdir -p "${output}"
mkdir -p log.out

#------------------------------------------------------------------------------
# Set up temporary working directory
#------------------------------------------------------------------------------

start_time=$(date +%s)

echo "Start time: $(date)"
echo "Present working directory: $(pwd)"

workdir=$(mktemp -d "${scratch}/tmpXXXXXX")

echo "Temporary directory: ${workdir}"

cd "${workdir}"

#------------------------------------------------------------------------------
# Select DNAm dataset for this array task
#------------------------------------------------------------------------------

ls "${gs_dnam}"/res_raw_beta_* > file_list.txt

dnam_dat=$(sed -n "${SLURM_ARRAY_TASK_ID}p" file_list.txt)
dnam_fname=$(basename "${dnam_dat}" .txt)

echo "Processing file: ${dnam_fname}"

#------------------------------------------------------------------------------
# Smoking phenotypes
#------------------------------------------------------------------------------

gs_cov_dat="${gs_cov_path}/GS_18869_all.cov"

# Ever-smoking status
awk '{print $2, $2, $6}' \
  "${gs_cov_dat}" > everSmoke.cov

# Smoking pack-years
awk '{print $2, $2, $7}' \
  "${gs_cov_dat}" > smokePackYrs.cov

#------------------------------------------------------------------------------
# 1. Convert DNAm data to OSCA binary format
#------------------------------------------------------------------------------

osca \
  --efile "${dnam_dat}" \
  --methylation \
  --make-bod \
  --out "${dnam_fname}"

#------------------------------------------------------------------------------
# 2. Update probe information
#------------------------------------------------------------------------------

osca \
  --befile "${dnam_fname}" \
  --update-opi "${gs_opi}/GS_sub_qc1.opi"

#------------------------------------------------------------------------------
# 3. Construct omics relationship matrix
#------------------------------------------------------------------------------

osca \
  --befile "${dnam_fname}" \
  --make-orm \
  --out "${dnam_fname}"

#------------------------------------------------------------------------------
# 4. OREML: ever smoking
#------------------------------------------------------------------------------

echo "Running OREML for ever smoking"

osca \
  --reml \
  --orm "${dnam_fname}" \
  --pheno everSmoke.cov \
  --qcovar "${gs_cov_path}/GS_18869_all_gs_environ_sub.cov" \
  --out "${dnam_fname}_everSmoke_osca"

#------------------------------------------------------------------------------
# 5. OREML: smoking pack-years
#------------------------------------------------------------------------------

echo "Running OREML for smoking pack-years"

osca \
  --reml \
  --orm "${dnam_fname}" \
  --pheno smokePackYrs.cov \
  --qcovar "${gs_cov_path}/GS_18869_all_gs_environ_sub.cov" \
  --out "${dnam_fname}_smokePackYrs_osca"

#------------------------------------------------------------------------------
# Save results
#------------------------------------------------------------------------------

cp ./*.log ./*.rsq "${output}/"

#------------------------------------------------------------------------------
# Clean up
#------------------------------------------------------------------------------

cd /
rm -rf "${workdir}"

end_time=$(date +%s)
execution_time=$((end_time - start_time))
execution_time_min=$((execution_time / 60))

echo "End time: $(date)"
echo "Total execution time: ${execution_time_min} minutes"

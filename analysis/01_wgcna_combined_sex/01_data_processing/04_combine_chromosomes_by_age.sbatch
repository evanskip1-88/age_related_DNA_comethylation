#!/bin/bash

#SBATCH --job-name=combine_age_groups
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=70G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=00:50:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err
#SBATCH --array=1-8


# ==============================================================================
# Combine chromosome-specific DNAm matrices within each age group
#
# One SLURM array task is run for each of the eight age groups.
# The age-group number is inserted into 04_combine_chromosomes_by_age.R by
# replacing the @AgeGrp@ placeholder.
#
# Inputs:
#   precorr_gs_qc_chr1_age_grp{AGE}.txt
#   ...
#   precorr_gs_qc_chr22_age_grp{AGE}.txt
#
# Output:
#   One whole-genome residualised DNAm matrix for each age group.
# ==============================================================================


# ------------------------------------------------------------------------------
# User-defined paths
# ------------------------------------------------------------------------------

output_path="..."
scratch="..."
R_SCRIPT="..."


# ------------------------------------------------------------------------------
# Age group supplied by SLURM array
# ------------------------------------------------------------------------------

AgeGrp=${SLURM_ARRAY_TASK_ID}

echo "Processing age group ${AgeGrp}"


# Exit immediately if a command fails
set -e


# ------------------------------------------------------------------------------
# Create temporary working directory
# ------------------------------------------------------------------------------

echo "Present working directory: $(pwd)"

workdir=$(mktemp -d "${scratch}/workdir/tmpXXXXXX")

echo "Temporary working directory: ${workdir}"

mkdir -p "${output_path}"

trap 'rm -rf "${workdir}"' EXIT

cd "${workdir}" || exit 1


# ------------------------------------------------------------------------------
# Prepare age-group-specific R script
# ------------------------------------------------------------------------------

sed \
    "s/@AgeGrp@/${AgeGrp}/g" \
    "${R_SCRIPT}/04_combine_chromosomes_by_age.R" \
    > run_script.R


# ------------------------------------------------------------------------------
# Run analysis
# ------------------------------------------------------------------------------

module load r/4.2.1-foss-2021a

Rscript --slave run_script.R


echo "Age group ${AgeGrp} complete."

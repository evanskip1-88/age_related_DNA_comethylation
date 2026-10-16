#!/bin/bash

#SBATCH --job-name=split_age_groups
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=30G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err


# ==============================================================================
# Split chromosome-specific residualised DNAm data into age groups
#
# Each SLURM job processes one chromosome. The chromosome number is inserted
# into 03_split_age_groups.R by replacing the @chr@ placeholder.
#
# Example submission:
#
# for chr in {1..22}; do
#     sbatch --export=chr=${chr} 03_split_age_groups.sbatch
# done
# ==============================================================================


# ------------------------------------------------------------------------------
# User-defined paths
# ------------------------------------------------------------------------------

scratch="..."
R_SCRIPT="..."


# Exit immediately if a command fails
set -e


# ------------------------------------------------------------------------------
# Create temporary working directory
# ------------------------------------------------------------------------------

echo "Present working directory: $(pwd)"

workdir=$(mktemp -d "${scratch}/tmpXXXXXX")

echo "Temporary working directory: ${workdir}"


# Ensure temporary files are removed when the job exits
trap 'rm -rf "${workdir}"' EXIT


# ------------------------------------------------------------------------------
# Prepare chromosome-specific R script
# ------------------------------------------------------------------------------

cd "${workdir}" || exit 1

echo "Processing chromosome ${chr}"

sed \
    "s/@chr@/${chr}/g" \
    "${R_SCRIPT}/03_split_age_groups.R" \
    > run_script.R


# ------------------------------------------------------------------------------
# Run analysis
# ------------------------------------------------------------------------------

module load r/4.2.1-foss-2021a

Rscript --slave run_script.R


echo "Chromosome ${chr} complete."

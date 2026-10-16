#!/bin/bash

#SBATCH --job-name=combine_residuals
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=00:50:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err


# ==============================================================================
# Combine block-wise residualised DNAm data by chromosome
#
# CpG residualisation is performed in separate blocks by
# 01_residualise_dnam.sbatch. This script combines those block-level residual
# files to produce one residualised DNAm matrix for each autosome.
#
# Inputs:
#   gs_qc_tr_lm_chr{CHR}.IDs
#   pre_corr_gs_*_chr{CHR}.res
#
# Outputs:
#   precorr_gs_qc_chr{CHR}.txt
#
# Each output matrix contains sample IDs in the first column followed by
# residualised CpG methylation values.
# ==============================================================================


# ------------------------------------------------------------------------------
# Defined path
# ------------------------------------------------------------------------------

output_path="..."

set -e

cd "${output_path}" || exit 1


# ------------------------------------------------------------------------------
# Combine residualised CpG blocks within each chromosome
# ------------------------------------------------------------------------------

echo "Combining residualised DNAm blocks by chromosome..."

for chr in {1..22}; do

    echo "Chromosome ${chr}"

    paste \
        "gs_qc_tr_lm_chr${chr}.IDs" \
        <(paste *chr${chr}.res) \
        > "precorr_gs_qc_chr${chr}.txt"

done


echo "All chromosome-specific residualised DNAm files created."

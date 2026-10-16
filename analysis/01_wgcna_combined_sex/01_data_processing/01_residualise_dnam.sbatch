#!/bin/bash

#SBATCH --job-name=residualise_dnam
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=70G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=00:40:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err


# ==============================================================================
# Residualise DNAm data in CpG blocks
#
# Each SLURM job processes one CpG block across chromosomes 1-22.
# Block definitions are obtained from gs_75_1kblocks.idx.
#
# Example submission:
#
# for block in {1..N}; do
#     sbatch --export=block=${block} 01_residualise_dnam.sbatch
# done
#
# Replace N with the number of blocks defined in the block index file.
# ==============================================================================


# ------------------------------------------------------------------------------
# User-defined paths
# ------------------------------------------------------------------------------

output_path="..."
dnam_inpath="..."
scratch="..."
R_SCRIPT="..."


# ------------------------------------------------------------------------------
# Create temporary working directory
# ------------------------------------------------------------------------------

echo "Present working directory: $(pwd)"

workdir=$(mktemp -d "${scratch}/tmpXXXXXX")

echo "Temporary working directory: ${workdir}"

mkdir -p "${output_path}"


# ------------------------------------------------------------------------------
# Copy chromosome-specific DNAm files to scratch
# ------------------------------------------------------------------------------

for chr in {1..22}; do

    echo "Copying chromosome ${chr}"

    cp \
        "${dnam_inpath}/chr${chr}_gs_qc_tr.txt" \
        "${workdir}/"

done


# ------------------------------------------------------------------------------
# Obtain CpG row range corresponding to the requested block
# ------------------------------------------------------------------------------

cd "${workdir}" || exit 1

block_range=$(awk -v block="${block}" \
    '$1 == block {print $2}' \
    "${dnam_inpath}/gs_75_1kblocks.idx")

echo "Block ${block}: ${block_range}"


# Insert block range into R script
sed \
    "s/@block@/${block_range}/g" \
    "${R_SCRIPT}/01_residualise_dnam.R" \
    > run_script.R


# ------------------------------------------------------------------------------
# Run residualisation
# ------------------------------------------------------------------------------

module load r/4.2.1-foss-2021a

Rscript --slave run_script.R 2>&1 | \
    tee "${output_path}/gs_qc_tr_lm_block${block}_chr_all.log"


# ------------------------------------------------------------------------------
# Copy results from scratch
# ------------------------------------------------------------------------------

rsync -av "${workdir}"/*.res "${output_path}/"
rsync -av "${workdir}"/*.IDs "${output_path}/"


# ------------------------------------------------------------------------------
# Clean temporary directory
# ------------------------------------------------------------------------------

cd ..
rm -rf "${workdir}"

echo "Block ${block} complete."

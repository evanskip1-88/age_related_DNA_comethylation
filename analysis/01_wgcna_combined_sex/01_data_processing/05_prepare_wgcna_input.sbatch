#!/bin/bash

#SBATCH --job-name=prepare_wgcna
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err
#SBATCH --array=1-8


# ==============================================================================
# Prepare age-stratified DNAm data for WGCNA
#
# Residualised whole-genome DNAm data were rounded to four decimal places to
# reduce file size and memory requirements. Each age-group matrix was then
# divided into two column-wise subsets to enable processing of the large
# methylation matrices on the HPC system.
#
# One SLURM array task is run for each of the eight age groups.
# ==============================================================================


# Age group supplied by SLURM array
AgeNum=${SLURM_ARRAY_TASK_ID}


# ------------------------------------------------------------------------------
# User-defined paths
# ------------------------------------------------------------------------------

dnam_inpath="..."
scratch="..."


# Exit immediately if a command fails
set -e


# ------------------------------------------------------------------------------
# Start
# ------------------------------------------------------------------------------

start_time=$(date +%s)

echo "Start time: $(date)"
echo "Processing age group ${AgeNum}"


# Create temporary working directory
workdir=$(mktemp -d "${scratch}/workdir/tmpXXXXXX")

echo "Temporary working directory: ${workdir}"

trap 'rm -rf "${workdir}"' EXIT

cd "${workdir}" || exit 1


# ------------------------------------------------------------------------------
# Step 1: Round residualised DNAm values
#
# Values are rounded to four decimal places to reduce storage and memory
# requirements when working with the full methylation matrix.
# ------------------------------------------------------------------------------

input_file="${dnam_inpath}/precorr_gs_qc_chr_all_age${AgeNum}.txt"
rounded_file="${dnam_inpath}/precorr_gs_qc_chr_all_age${AgeNum}_roundoff.txt"

awk '
    NR == 1 {
        print
        next
    }
    {
        printf "%s %.4f %s %.4f ", $1, $2, $3, $4

        for (i = 5; i <= NF; i++) {
            printf "%.4f%s", $i, (i < NF ? OFS : ORS)
        }
    }
' "${input_file}" > "${rounded_file}"


# ------------------------------------------------------------------------------
# Step 2: Split the rounded matrix into two column-wise subsets
#
# The first four columns contain sample/covariate information and are retained
# in both subsets.
# ------------------------------------------------------------------------------

echo "Splitting age group ${AgeNum} matrix..."


# Subset 1
awk '{
    print $1, $2, $3, $4, substr($0, index($0, $5))
}' "${rounded_file}" |
    cut -d' ' -f1-400000 \
    > "${dnam_inpath}/sub1_precorr_gs_qc_chr_all_age${AgeNum}_roundoff.txt"


# Subset 2
awk '{
    print $1, $2, $3, $4, substr($0, index($0, $5))
}' "${rounded_file}" |
    cut -d' ' -f1-4,400001- \
    > "${dnam_inpath}/sub2_precorr_gs_qc_chr_all_age${AgeNum}_roundoff.txt"


# ------------------------------------------------------------------------------
# Finish
# ------------------------------------------------------------------------------

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

echo "Age group ${AgeNum} complete."
echo "End time: $(date)"
echo "Elapsed time: ${elapsed_time} seconds"

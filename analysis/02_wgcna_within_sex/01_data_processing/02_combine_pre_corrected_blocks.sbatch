#!/bin/bash
#SBATCH --job-name=comb.pre.corr.dat
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=00:50:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err

set -euo pipefail

# ============================================================
# Combine pre-corrected DNAm blocks within each sex
# Combines residuals across CpG blocks for each chromosome.
# ============================================================

# Configure HPC directories
output_path="..."

mkdir -p "$output_path"
cd "$output_path"

# Combine residuals separately for males and females
for sex in M F; do

    echo "Combining residuals for sex: $sex"

    for chr in {1..22}; do

        echo "Processing chromosome: $chr"

        # Find residual files for this chromosome and sex
        res_files=( *chr${chr}_${sex}.res )

        # Check that residual files exist
        if (( ${#res_files[@]} == 0 )); then
            echo "ERROR: No residual files for chr${chr}_${sex}"
            exit 1
        fi

        # Combine sample IDs and residuals
        paste \
            "gs_qc_tr_lm_chr${chr}_${sex}.IDs" \
            <(paste "${res_files[@]}") \
            > "precorr_gs_qc_chr${chr}_${sex}.txt"

    done
done

echo "All pre-corrected blocks combined successfully!"

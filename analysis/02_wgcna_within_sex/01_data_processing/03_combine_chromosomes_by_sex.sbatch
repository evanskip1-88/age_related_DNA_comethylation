#!/bin/bash
#SBATCH --job-name=comb.chr.by.sex
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=70G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=00:50:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err

set -euo pipefail

# ============================================================
# Combine pre-corrected DNAm data across chromosomes by sex
# Verify sample ID consistency before combining.
# ============================================================

# Configure HPC directory
output_path="..."

cd "$output_path"

# Process males and females separately
for sex in M F; do

    echo "Combining chromosome files for sex: $sex"

    # Extract sample IDs from chromosome 1
    cut -f1 "precorr_gs_qc_chr1_${sex}.txt" > "IID_${sex}.txt"

    # Check sample ID consistency and extract DNAm residuals
    chr_files=()

    for chr in {1..22}; do

        echo "Processing chromosome: $chr"

        input_file="precorr_gs_qc_chr${chr}_${sex}.txt"
        cut_file="chr${chr}_${sex}_cut.txt"

        # Verify identical sample IDs and ordering
        if ! diff -q \
            "IID_${sex}.txt" \
            <(cut -f1 "$input_file") > /dev/null; then

            echo "ERROR: Sample ID mismatch in chromosome $chr ($sex)"
            exit 1
        fi

        # Remove sample ID column
        cut -f2- "$input_file" > "$cut_file"

        chr_files+=("$cut_file")

    done

    # Combine chromosome data in numerical order
    paste "IID_${sex}.txt" "${chr_files[@]}" \
        > "precorr_gs_qc_chr_all_${sex}.txt"

    # Remove temporary files
    rm -f "IID_${sex}.txt" "${chr_files[@]}"

    echo "Completed: precorr_gs_qc_chr_all_${sex}.txt"

done

echo "All chromosome files combined successfully."

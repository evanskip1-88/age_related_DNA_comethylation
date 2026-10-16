#!/bin/bash

set -euo pipefail

# ============================================================
# Combine heteroskedasticity results across chromosomes
# and reference age groups (Age2-Age8).
#
# Retain the header from the first file only.
# ============================================================

input_path="..."
output_path="/results/03_pairwise_correlations"

mkdir -p "$output_path"

output_file="${output_path}/chr_all_AgeBaseAll_age_beta_est_resiq_sq.res"

files=()

for chr in {1..22}; do
    for age in {2..8}; do
        file="${input_path}/chr${chr}_AgeBase${age}_age_beta_est_resiq_sq.txt"

        if [[ ! -f "$file" ]]; then
            echo "ERROR: Missing result file: $file" >&2
            exit 1
        fi

        files+=("$file")
    done
done

awk '
    FNR == 1 && NR != 1 {next}
    {print}
' "${files[@]}" > "$output_file"

echo "Combined ${#files[@]} files."
echo "Output: $output_file"

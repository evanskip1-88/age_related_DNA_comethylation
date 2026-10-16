#!/bin/bash

set -euo pipefail

input_path="/results/03_pairwise_correlations"
output_path="/results/03_pairwise_correlations"

mkdir -p "$output_path"

output_file="${output_path}/precorr_gs_qc_chr_all_decay_cor.res"

shopt -s nullglob

files=("${input_path}"/precorr_gs_qc_chr*_Age*_decay.cor)

if [[ ${#files[@]} -eq 0 ]]; then
    echo "ERROR: No correlation-decay results found."
    exit 1
fi

awk '
    FNR == 1 && NR != 1 {next}
    {print}
' "${files[@]}" > "$output_file"

echo "Combined ${#files[@]} correlation-decay files."
echo "Output: $output_file"

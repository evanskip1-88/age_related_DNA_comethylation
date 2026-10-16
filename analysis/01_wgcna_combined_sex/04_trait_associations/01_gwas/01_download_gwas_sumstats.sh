#!/usr/bin/env bash
set -euo pipefail

trait_file="/data/gwas/02_gwas_catalog_traits.txt"
output_dir="/data/gwas"

mkdir -p "$output_dir"

for row in $(seq 2 91); do
    trait_url=$(awk -F '\t' -v nr="$row" 'NR==nr {print $5}' "$trait_file")

    if [[ -z "$trait_url" ]]; then
        echo "Missing URL at row $row" >&2
        exit 1
    fi

    wget -P "$output_dir" "$trait_url"
done

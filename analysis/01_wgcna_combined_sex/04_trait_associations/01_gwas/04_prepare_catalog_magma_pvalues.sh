#!/usr/bin/env bash
set -euo pipefail

# Prepare curated GWAS Catalog summary statistics for MAGMA.
# Input files use GRCh37 coordinates.
#
# Expected source columns:
#   1 = CHR, 2 = SNP, 3 = BP, 6 = SS, 13 = P
#
# Retain variants with rs identifiers, as in the original analysis.

input_dir="data/gwas/sumstats"
output_dir="data/gwas/magma"

mkdir -p "$output_dir"

trait_files=(
    "GCST90041821_buildGRCh37.tsv.gz"
    "GCST90041858_buildGRCh37.tsv.gz"
)

for file in "${trait_files[@]}"; do

    trait_id="${file%%_*}"

    echo "Preparing ${trait_id}"

    gzip -dc "${input_dir}/${file}" |
        awk -v trait="$trait_id" '
            BEGIN {
                OFS = "\t"
                print "SNP", "CHR", "BP", "LP", "P", "SS", "traitID"
            }
            NR > 1 && $2 ~ /^rs/ {
                print $2, $1, $3, -log($13)/log(10), $13, $6, trait
            }
        ' > "${output_dir}/${trait_id}.pval"

done


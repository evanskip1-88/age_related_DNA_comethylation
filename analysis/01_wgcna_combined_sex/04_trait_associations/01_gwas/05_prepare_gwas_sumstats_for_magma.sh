#!/usr/bin/env bash
set -euo pipefail

# Prepare OpenGWAS VCF summary statistics for MAGMA.
# Original dataset: ukb-b-8909
# Output: SNP CHR BP LP P SS traitID

input_dir="/data/gwas/sumstats"
output_dir="/data/gwas/magma"
traits_file="/data/gwas/03_open_gwas_traits.txt"

mkdir -p "$output_dir"

gwas_id="ukb-b-8909"
input_file="${input_dir}/${gwas_id}.vcf.gz"
output_file="${output_dir}/${gwas_id}.pval"

# Sample size from the original trait manifest (column 3)
sample_size=$(awk -v id="$gwas_id" '$1 == id {print $3; exit}' "$traits_file")

if [[ -z "$sample_size" ]]; then
    echo "Sample size not found for ${gwas_id}" >&2
    exit 1
fi

echo "Preparing ${gwas_id}"

gzip -dc "$input_file" |
    awk -v sample_size="$sample_size" -v trait_id="$gwas_id" '
        BEGIN {
            OFS = "\t"
            print "SNP", "CHR", "BP", "LP", "P", "SS", "traitID"
        }
        !/^#/ {
            split($10, fields, ":")
            lp = fields[3]
            if ($3 ~ /^rs/ && lp != "" && lp != ".") {
                print $3, $1, $2, lp, 10^(-lp), sample_size, trait_id
            }
        }
    ' > "$output_file"

echo "Saved: $output_file"


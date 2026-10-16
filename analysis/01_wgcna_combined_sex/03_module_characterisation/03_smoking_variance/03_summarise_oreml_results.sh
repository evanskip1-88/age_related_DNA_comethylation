#!/bin/bash

#==============================================================================
# Summarise OSCA OREML results
#
# Extract variance components, variance explained (R2), P-values and
# standard errors from the OSCA .rsq output files.
# 
# NB: Run this script within the OSCA-OREML out put folder
#==============================================================================

set -euo pipefail

output_file="results/01_wgcna_combined_sex/01a_oreml_variance_summary.tsv"

# Write header
echo -e "analysis\tV_o\tV_e\tV_p\tR2\tPval\tSE" > "${output_file}"

#------------------------------------------------------------------------------
# Extract results from each OREML output
#------------------------------------------------------------------------------

for file in raw_beta_*.rsq res_*.rsq; do

    fname=$(basename "${file}" .rsq)

    # Remove DNAm input prefixes
    fname=${fname#raw_beta_}
    fname=${fname#res_raw_beta_}

    # Extract variance components and model statistics
    Vo=$(awk '$1=="V(O)" {print $2}' "${file}")
    Ve=$(awk '$1=="V(e)" {print $2}' "${file}")
    Vp=$(awk '$1=="Vp" {print $2}' "${file}")
    R2=$(awk '$1=="V(O)/Vp" {print $2}' "${file}")
    Pval=$(awk '$1=="Pval" {print $2}' "${file}")
    SE=$(awk '$1=="V(O)/Vp" {print $3}' "${file}")

    echo -e \
      "${fname}\t${Vo}\t${Ve}\t${Vp}\t${R2}\t${Pval}\t${SE}" \
      >> "${output_file}"

done

echo "OREML summary written to: ${output_file}"

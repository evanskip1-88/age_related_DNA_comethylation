#!/usr/bin/env bash
#SBATCH --job-name=magma_gene_analysis
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=05:00:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err
#SBATCH --array=3-91

# =================================================================================
# Purpose: Run MAGMA SNP-to-gene annotation and gene-level GWAS association
#          analysis using GRCh37 coordinates and a European LD reference.
#
# Inputs:  GWAS p-value files, trait manifest, GRCh37 gene locations,
#          and 1000 Genomes European reference panel.
# Outputs: MAGMA gene annotation (.genes.annot) and association (.genes.out) files.
# =================================================================================

set -euo pipefail

start_time=$(date +%s)
echo "Started: $(date)"

# Run from repository root
cd "${SLURM_SUBMIT_DIR}"

# Configure MAGMA executable and reference data
MAGMA="<path_to_magma_executable>"

gene_loc="/data/gwas/NCBI37.3.gene.loc"
reference="/data/gwas/reference/g1000_eur"

# Input and output directories
pval_dir="/data/gwas/magma"
trait_list="${pval_dir}/all_pval_file_names.txt"

annot_dir="results/trait_associations/gwas/gene_annotation"
gene_dir="results/trait_associations/gwas/magma_gene_analysis"

mkdir -p "$annot_dir" "$gene_dir"

# Select GWAS trait from the original manifest
trait_n="${SLURM_ARRAY_TASK_ID}"

gwas_file=$(awk -v n="$trait_n" 'NR == n {print $2}' "$trait_list")

if [[ -z "$gwas_file" ]]; then
    echo "No GWAS filename found at row ${trait_n}" >&2
    exit 1
fi

gwas_id="${gwas_file%.pval}"
pval_file="${pval_dir}/${gwas_id}.pval"

echo "Analysing GWAS: ${gwas_id}"

# Step 1: SNP-to-gene annotation
"$MAGMA" \
    --annotate window=50 \
    --snp-loc "$pval_file" \
    --gene-loc "$gene_loc" \
    --out "${annot_dir}/${gwas_id}"

# Step 2: MAGMA gene-level association analysis
"$MAGMA" \
    --bfile "$reference" \
    --pval "$pval_file" ncol=SS \
    --gene-annot "${annot_dir}/${gwas_id}.genes.annot" \
    --genes-only \
    --out "${gene_dir}/${gwas_id}"

echo "Finished: $(date)"
echo "Execution time: $(( ($(date +%s) - start_time) / 60 )) minutes"


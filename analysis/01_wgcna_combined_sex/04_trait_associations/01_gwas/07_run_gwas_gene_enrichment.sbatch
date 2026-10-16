#!/usr/bin/env bash
#SBATCH --job-name=gwas_gene_enrichment
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=30G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err

# =============================================================================
# Purpose: Test enrichment of methylation-associated genes among top MAGMA
#          GWAS genes using hypergeometric tests.
#
# Inputs:  WGCNA gene lists, CpG-to-gene mappings, MAGMA gene-level results.
# Outputs: GWAS enrichment statistics and overlapping gene lists.
# =============================================================================

set -euo pipefail

start_time=$(date +%s)
echo "Started: $(date)"

# Run from the repository root
cd "${SLURM_SUBMIT_DIR}"

# Load R
module load r

# Run GWAS gene enrichment analysis
Rscript analysis/04_trait_associations/01_gwas/07_run_gwas_gene_enrichment.R

echo "Finished: $(date)"
echo "Execution time: $(( ($(date +%s) - start_time) / 60 )) minutes"


#!/bin/bash
#SBATCH --job-name=GS.heteroskedasticity
#SBATCH --account=<HPC_ACCOUNT>
#SBATCH --partition=<HPC_PARTITION>
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00:00
#SBATCH --array=2-8
#SBATCH --output=log.out/%x.%A_%a.out
#SBATCH --error=log.out/%x.%A_%a.err

set -euo pipefail

# ============================================================
# Age-related heteroskedasticity of CpG-pair relationships
#
# For each strongly correlated CpG pair:
#
# Model 1:
#   CpG1 ~ CpG2 + age_group
#
# Normalised squared residuals:
#   residual^2 / mean(residual^2)
#
# Model 2:
#   normalised_squared_residuals ~ age_group
#
# Each SLURM array task processes one reference age group
# for a specified chromosome.
# ============================================================

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

dnam_path="..."
cov_inpath="..."
probe_inpath="..."
output_path="..."
script_path="/analysis/03_pairwise_correlations/02_correlation_analysis"

chr="${chr:?ERROR: Supply chr=1..22}"
AgeGrp="${SLURM_ARRAY_TASK_ID}"

# ------------------------------------------------------------
# Validate inputs
# ------------------------------------------------------------

if ! [[ "$chr" =~ ^([1-9]|1[0-9]|2[0-2])$ ]]; then
    echo "ERROR: Invalid chromosome: $chr"
    exit 1
fi

pair_file="${probe_inpath}/precorr_gs_qc_chr${chr}_Age${AgeGrp}.cor"
dnam_file="${dnam_path}/precorr_gs_qc_chr${chr}.txt"
cov_file="${cov_inpath}/GS_18869_all.cov"

for file in "$pair_file" "$dnam_file" "$cov_file"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: Missing input file: $file"
        exit 1
    fi
done

mkdir -p "$output_path"

# ------------------------------------------------------------
# Load R
# ------------------------------------------------------------

module load r/4.4.2-gfbf-2024a

echo "Analysis started: $(date)"
echo "Chromosome: $chr"
echo "Reference age group: Age${AgeGrp}"

# ------------------------------------------------------------
# Run heteroskedasticity analysis
# ------------------------------------------------------------

Rscript \
    "${script_path}/02_test_pairwise_heteroskedasticity.R" \
    "$chr" \
    "$AgeGrp" \
    "$dnam_path" \
    "$cov_inpath" \
    "$probe_inpath" \
    "$output_path"

echo "Analysis completed: $(date)"

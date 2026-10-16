#!/bin/bash
#SBATCH --job-name=GS.corr_decay
#SBATCH --account=<HPC_ACCOUNT>
#SBATCH --partition=<HPC_PARTITION>
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err

# ============================================================
# Correlation decay analysis
#
# One job per chromosome-by-age DNAm file.
#
# Usage:
#   mkdir -p log.out
#
#   for num in {1..N}; do
#       sbatch --export=ALL,fnum="$num" \
#           04_calculate_correlation_decay.sbatch
#   done
#
# Replace N with the number of prepared input files.
# ============================================================

set -euo pipefail

# ------------------------------------------------------------
# 1. Configure paths
# ------------------------------------------------------------

dnam_inpath="..." # Indvidual level data
opi_file="/data/gsDNAm"
output_path="/results/03_pairwise_correlations"
scratch="..."
script_path="/analysis/03_pairwise_correlations/02_correlation_analysis/"

fnum="${fnum:?ERROR: Supply fnum when submitting the job}"

mkdir -p "$output_path"

# ------------------------------------------------------------
# 2. Identify the input file
# ------------------------------------------------------------

shopt -s nullglob

matches=(
    "${dnam_inpath}/f${fnum}_precorr_gs_qc_chr"*.txt
)

if [[ ${#matches[@]} -ne 1 ]]; then
    echo "ERROR: Expected exactly one input file for fnum=${fnum}"
    echo "Number of matches: ${#matches[@]}"
    exit 1
fi

input_file="${matches[0]}"
filename="$(basename "$input_file")"
basename_no_ext="${filename%.txt}"

echo "Input file: $input_file"

# Original naming convention:
# f<num>_precorr_gs_qc_chr<chr>_Age<age>.txt
#
# Extract chromosome and age from the filename.

if [[ "$filename" =~ _chr([0-9]+)_Age([1-8])\.txt$ ]]; then
    chr="${BASH_REMATCH[1]}"
    age_group="${BASH_REMATCH[2]}"
else
    echo "ERROR: Could not extract chromosome and age from: $filename"
    exit 1
fi

echo "Chromosome: $chr"
echo "Age group: $age_group"

# ------------------------------------------------------------
# 3. Create temporary working directory
# ------------------------------------------------------------

mkdir -p "$scratch"

workdir="$(mktemp -d "${scratch}/tmpXXXXXX")"

trap 'rm -rf -- "$workdir"' EXIT

echo "Temporary directory: $workdir"

cp "$input_file" "$workdir/"
cp "$opi_file" "$workdir/"

# ------------------------------------------------------------
# 4. Load R
# ------------------------------------------------------------

module load r/4.4.2-gfbf-2024a

# ------------------------------------------------------------
# 5. Run correlation-decay analysis
# ------------------------------------------------------------

Rscript "${script_path}/04_calculate_correlation_decay.R" \
    "${workdir}/${filename}" \
    "${workdir}/$(basename "$opi_file")" \
    "$chr" \
    "$age_group" \
    "$workdir"

# ------------------------------------------------------------
# 6. Copy results to output directory
# ------------------------------------------------------------

result_file="${workdir}/precorr_gs_qc_chr${chr}_Age${age_group}_decay.cor"

if [[ ! -f "$result_file" ]]; then
    echo "ERROR: Expected result file not found: $result_file"
    exit 1
fi

cp "$result_file" "$output_path/"

echo "Correlation-decay analysis completed successfully."
echo "Result: ${output_path}/$(basename "$result_file")"

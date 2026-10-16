#!/bin/bash
#SBATCH --job-name=GS.chr.export
#SBATCH --account=<HPC_ACCOUNT>
#SBATCH --partition=<HPC_PARTITION>
#SBATCH --mem=150G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=100:00:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err

set -euo pipefail

# ============================================================
# Export QC-filtered Generation Scotland DNA methylation data
# separately for chromosomes 1-22 using OSCA.
#
# Input:
#   Existing OSCA binary dataset: GS_sub_qc1
#
# Output:
#   chr1_gs_qc.txt ... chr22_gs_qc.txt
#
# Note:
#   The upstream QC step is not rerun by this script.
# ============================================================

# Configure paths
scratch="..."
input_path="..."
output_path="..."

input_prefix="GS_sub_qc1"

# Chromosome passed at submission
chr="${chr:?Specify chr=1..22 when submitting}"

if ! [[ "$chr" =~ ^([1-9]|1[0-9]|2[0-2])$ ]]; then
    echo "ERROR: Chromosome must be an integer from 1 to 22."
    exit 1
fi

echo "Analysis started: $(date)"
echo "Chromosome: ${chr}"

# Prepare directories
mkdir -p "$scratch" "$output_path"

workdir=$(mktemp -d "${scratch}/tmpXXXXXX")
trap 'rm -rf "$workdir"' EXIT

cd "$workdir"

# OSCA must be available in PATH
command -v osca >/dev/null 2>&1 || {
    echo "ERROR: OSCA executable not found."
    exit 1
}

# Export non-transposed methylation data
osca \
    --befile "${input_path}/${input_prefix}" \
    --chr "$chr" \
    --make-efile \
    --thread-num 2 \
    --out "${output_path}/chr${chr}_gs_qc.txt"

echo "Chromosome ${chr} completed: $(date)"

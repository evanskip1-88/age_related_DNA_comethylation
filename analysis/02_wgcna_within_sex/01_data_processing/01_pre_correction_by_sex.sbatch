#!/bin/bash
#SBATCH --job-name=lm.gs.by.sex
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=70G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=00:50:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err

# ============================================================
# DNAm pre-correction within sex
# Residualise DNAm values separately by sex and chromosome
# using predefined CpG blocks.
#
# Example submission:
# for block in {1..50}; do
#     sbatch --export=ALL,block=$block 01_pre_correction_by_sex.sbatch
# done
#
# ============================================================

# Configure HPC directories
output_path="..."
dnam_inpath="..."
scratch="..."
RSCRIPT="..."

# Check block identifier
: "${block:?Please supply block using sbatch --export=ALL,block=N}"

echo "Present working directory: $(pwd)"

# Create temporary working directory
workdir=$(mktemp -d "${scratch}/tmpXXXXXX")
echo "Temporary directory: $workdir"

mkdir -p "$output_path"

# Retrieve the genomic block range
block_range=$(awk -v block="$block" \
    '$1 == block {print $2}' \
    "${dnam_inpath}/gs_75_1kblocks.idx")

if [[ -z "$block_range" ]]; then
    echo "ERROR: No block range found for block $block"
    exit 1
fi

echo "Processing block: $block_range"

# Substitute block placeholder in R script
sed "s/@block@/${block_range}/g" \
    "${RSCRIPT}/01_pre_correction_by_sex.R" \
    > "${workdir}/run_script.R"

# Load R
module load r/4.4.2-gfbf-2024a

# Run analysis
cd "$workdir"

Rscript --slave run_script.R \
    2>&1 | tee "${output_path}/gs_qc_tr_lm_${block_range}_chr_all.log"

# Copy residuals and sample IDs
shopt -s nullglob

res_files=( *.res )
id_files=( *.IDs )

if (( ${#res_files[@]} > 0 )); then
    rsync -av "${res_files[@]}" "$output_path/"
fi

if (( ${#id_files[@]} > 0 )); then
    rsync -av "${id_files[@]}" "$output_path/"
fi

# Clean temporary directory
cd /
rm -rf "$workdir"

echo "Analysis completed for block: $block_range"

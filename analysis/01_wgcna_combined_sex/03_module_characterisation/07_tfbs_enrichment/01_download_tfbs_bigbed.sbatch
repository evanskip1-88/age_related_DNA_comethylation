#!/bin/bash
#SBATCH --job-name=download_TFBS
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err

#==============================================================================
# Download JASPAR 2024 transcription factor binding-site annotations
#
# Genome assembly: hg19
# Format: bigBed (.bb)
#==============================================================================

set -euo pipefail

# Output directory
output_path="/data/encode_tfbs"

mkdir -p "${output_path}"

# Download JASPAR TFBS annotations
wget \
  --continue \
  --show-progress \
  --directory-prefix="${output_path}" \
  "https://frigg.uio.no/JASPAR/JASPAR_TFBSs/2024/JASPAR2024_hg19.bb"

echo "JASPAR 2024 hg19 TFBS download completed."

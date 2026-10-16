#!/bin/bash
#SBATCH --job-name=rGREAT.wg
#SBATCH --account=a_mcrae
#SBATCH --partition=general
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=05:00:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err

# =============================================================================
# Run rGREAT functional enrichment on Bunya
#
# Submit from the repository root:
#   mkdir -p log.out
#   sbatch analysis/03_module_characterisation/10_functional_enrichment/02_submit_rgreat.slurm
# =============================================================================

set -euo pipefail

echo "Start time: $(date)"

# Paths
REPO_DIR="..."
SCRATCH_DIR="..."

RSCRIPT="${REPO_DIR}/analysis/03_module_characterisation/10_functional_enrichment/01_run_rgreat_enrichment.R"

# Temporary working directory
mkdir -p "${SCRATCH_DIR}/workdir"
WORKDIR=$(mktemp -d "${SCRATCH_DIR}/workdir/tmpXXXXXX")

trap 'rm -rf "${WORKDIR}"' EXIT

echo "Temporary directory: ${WORKDIR}"

module load r/4.4.2-gfbf-2024a

cd "${REPO_DIR}"

Rscript "${RSCRIPT}"

echo "End time: $(date)"


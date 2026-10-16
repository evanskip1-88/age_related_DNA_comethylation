
#!/usr/bin/env bash
#SBATCH --job-name=cpgs2genes
#SBATCH --account=<project_account>
#SBATCH --partition=<partition>
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --output=log.out/%x.%j.out
#SBATCH --error=log.out/%x.%j.err

# Map CpGs to NCBI build 37 genes using the original
# coordinate window (gene start - 5 kb to gene end + 1 kb).

set -euo pipefail

start_time=$(date +%s)
echo "Started: $(date)"

# Run from the repository root
cd "${SLURM_SUBMIT_DIR}"

module load r/4.4.2-gfbf-2024a

Rscript analysis/04_trait_associations/01_gwas/03_map_cpgs_to_magma_genes.R

echo "Finished: $(date)"
echo "Execution time: $(( ($(date +%s) - start_time) / 60 )) minutes"


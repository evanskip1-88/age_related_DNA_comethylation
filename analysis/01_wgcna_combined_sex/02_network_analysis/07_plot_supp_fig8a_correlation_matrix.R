# =============================================================================
#
# Purpose: Plot correlation matrix of WGCNA preservation statistics (Fig. S8a).
# Inputs:  WGCNA observed preservation statistics.
# Outputs: Supplementary Figure 8a.
#
# =============================================================================

library(data.table)
library(pheatmap)

# Input and output
input_file <- "results/01_wgcna_combined_sex/precorr_gs_qc_chr_all_refage1_wg.obsStat"
output_file <- "figures/Figure_S8a_correlation_matrix.pdf"

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

# Load observed preservation statistics
dt_obs_stats <- fread(input_file)

# Correlation matrix of preservation statistics
dt_obs_stats_cor <- cor(dt_obs_stats[, 5:8])

# Supplementary Figure 8a
pdf(output_file, width = 7, height = 7)

pheatmap(
  dt_obs_stats_cor,
  angle_col = 45,
  scale = "none",
  fontsize = 16,
  display_numbers = TRUE,
  number_format = "%.2f"
)

dev.off()


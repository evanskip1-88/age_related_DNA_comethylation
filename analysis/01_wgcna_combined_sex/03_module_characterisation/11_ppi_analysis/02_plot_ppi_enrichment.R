
# =============================================================================
# Purpose: Visualise STRING PPI enrichment across stable and changing modules.
# Inputs:  PPI enrichment results from 01_run_ppi_analysis.R.
# Outputs: Figure 2d – module-level PPI enrichment plot.
# =============================================================================

library(data.table)
library(readxl)
library(ggplot2)

# -------------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------------

input_file <- "/results/01_wgcna_combined_sex/02a_ppi_enrichment.xlsx"
output_file <- "/figures/Figure_2d_PPI_enrichment.pdf"

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Load and prepare PPI results
# -------------------------------------------------------------------------

dt_res <- as.data.table(read_xlsx(input_file))

# Cap -log10(P) at 30
dt_res[, logPval := -log10(pmax(ppi_pval.enrichment, 1e-300))]
dt_res[, logPval := pmin(logPval, 30)]

# Exclude pooled module groups and retain original module ordering
dt_res_sub <- dt_res[
  analysis != "combined_mods"
][order(logPval)]

dt_res_sub <- dt_res_sub[
  order(analysis, logPval)
]

dt_res_sub[
  ,
  moduleNumID := factor(moduleNumID, levels = unique(moduleNumID))
]

# -------------------------------------------------------------------------
# Plot module-level PPI enrichment
# -------------------------------------------------------------------------

p <- ggplot(
  dt_res_sub[n_input_genes >= 3],
  aes(
    x = moduleNumID,
    y = logPval,
    color = analysis,
    size = n_input_genes
  )
) +
  geom_segment(
    aes(
      xend = moduleNumID,
      y = 0,
      yend = logPval
    ),
    color = "grey60",
    linewidth = 0.5,
    alpha = 0.7
  ) +
  geom_point(alpha = 0.9) +
  coord_flip() +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.5
  ) +
  labs(
    x = "Module ID",
    y = "PPI enrichment",
    size = "Gene count",
    color = NULL
  ) +
  scale_size(range = c(2, 6)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = c(0.8, 0.25),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.5),
    plot.margin = margin(5, 10, 5, 5),
    panel.grid = element_blank()
  )

# -------------------------------------------------------------------------
# Export
# -------------------------------------------------------------------------

ggsave(
  filename = output_file,
  plot = p,
  width = 8,
  height = 7,
  units = "in"
)

print(p)


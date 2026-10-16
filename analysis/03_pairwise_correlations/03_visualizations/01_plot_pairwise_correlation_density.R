#!/usr/bin/env Rscript

# =============================================================================
# Figure 1: Pairwise CpG correlation dynamics across age groups
#
# Input:
#   precorr_gs_qc_chr_all_Ages_vsOthers_long.cor
#
# Output:
#   Fig_1_pairwise_correlations_positive.png
#   Pairwise_correlations_negative.png
#
# Description:
#   Visualises correlations of CpG pairs identified in Age1 across
#   Age2-Age8 using two-dimensional density heatmaps.
#
#   Positive and negative Age1 correlations are plotted separately.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# 1. Configure paths
# -----------------------------------------------------------------------------

input_file <- "/results/03_pairwise_correlations/precorr_gs_qc_chr_all_AgebaseAll_vsOthers.cor"
output_dir <- "/figures"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

options(bitmapType = "cairo")

# -----------------------------------------------------------------------------
# 2. Load and prepare data
# -----------------------------------------------------------------------------

dt_corr <- fread(input_file)

required_cols <- c("AgeBase", "Age.grp", "AgeBaseCorr", "corr")

stopifnot(all(required_cols %in% names(dt_corr)))

# Retain correlations identified using Age1 as the reference
dt_age1 <- dt_corr[
  AgeBase == "Age1" &
    Age.grp %in% paste0("Age", 2:8) &
    !is.na(AgeBaseCorr) &
    !is.na(corr)
]

if (nrow(dt_age1) == 0L) {
  stop("No valid Age1 reference correlations found.")
}

dt_age1[, Age.grp := factor(
  Age.grp,
  levels = paste0("Age", 2:8)
)]

# -----------------------------------------------------------------------------
# 3. Common plotting function
# -----------------------------------------------------------------------------

plot_pairwise_density <- function(
    dt,
    x_breaks,
    x_limits,
    y_limits,
    y_label
) {

  ggplot(dt, aes(x = AgeBaseCorr, y = corr)) +
    geom_bin2d(bins = 120) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      colour = "black"
    ) +
    facet_wrap(~Age.grp, ncol = 7) +
    scale_x_continuous(
      breaks = x_breaks,
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    coord_cartesian(
      xlim = x_limits,
      ylim = y_limits
    ) +
    scale_fill_viridis_c(
      trans = "log10",
      option = "plasma"
    ) +
    labs(
      x = "Correlation (Age1)",
      y = y_label,
      fill = "CpG pairs"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      legend.position = "right",
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )
}

# -----------------------------------------------------------------------------
# 4. Positive correlations: Figure 1
# -----------------------------------------------------------------------------

dt_pos <- dt_age1[AgeBaseCorr > 0]

if (nrow(dt_pos) == 0L) {
  stop("No positively correlated Age1 CpG pairs found.")
}

p_pos <- plot_pairwise_density(
  dt = dt_pos,
  x_breaks = c(0.8, 0.9, 1),
  x_limits = c(0.78, 1),
  y_limits = c(min(dt_pos$corr) - 0.1, 1),
  y_label = "Correlation (Ages 2-8)"
)

ggsave(
  filename = file.path(
    output_dir,
    "Fig_1_pairwise_correlations_positive.png"
  ),
  plot = p_pos,
  width = 16,
  height = 5,
  dpi = 300
)

# -----------------------------------------------------------------------------
# 5. Negative correlations
# -----------------------------------------------------------------------------

dt_neg <- dt_age1[AgeBaseCorr < 0]

if (nrow(dt_neg) > 0L) {

  p_neg <- plot_pairwise_density(
    dt = dt_neg,
    x_breaks = c(-1, -0.9, -0.8),
    x_limits = c(-1, -0.78),
    y_limits = c(-1, max(dt_neg$corr) + 0.1),
    y_label = "Correlation (older age groups)"
  )

  ggsave(
    filename = file.path(
      output_dir,
      "Pairwise_correlations_negative.png"
    ),
    plot = p_neg,
    width = 16,
    height = 5,
    dpi = 300
  )
}

message("Figure generation completed.")

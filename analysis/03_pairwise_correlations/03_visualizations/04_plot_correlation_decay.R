#!/usr/bin/env Rscript

# ============================================================
# Supplementary Figure S5b
# Pairwise CpG correlation decay with genomic distance
#
# Input:
#   precorr_gs_qc_chr_all_decay_cor.res
#
# Analysis:
#   Age1 reference
#   Genomic distances <= 10 kb
#   Positive and negative correlation groups
#
# Plot:
#   Mean correlation +/- SD
#   Zoomed region: 0-2 kb
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggforce)
})

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

input_file <- "/results/03_pairwise_correlations/precorr_gs_qc_chr_all_decay_cor.res"
output_dir <- "/figures"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Load correlation decay summaries
# ------------------------------------------------------------

dt_corr <- fread(input_file)

required_cols <- c(
  "age.grp",
  "dist.mean.bin",
  "corr.mean.bin",
  "corr.sd.bin",
  "corr.grp"
)

stopifnot(all(required_cols %in% names(dt_corr)))

# ------------------------------------------------------------
# 3. Select Age1 and distances <= 10 kb
# ------------------------------------------------------------

dt_age1 <- dt_corr[
  age.grp == "Age1" &
    dist.mean.bin <= 10000 &
    !is.na(dist.mean.bin) &
    !is.na(corr.mean.bin) &
    !is.na(corr.sd.bin) &
    !is.na(corr.grp)
]

if (nrow(dt_age1) == 0L) {
  stop("No valid Age1 correlation decay data found.")
}

# ------------------------------------------------------------
# 4. Generate Supplementary Figure S5b
# ------------------------------------------------------------

p_s5b <- ggplot(
  dt_age1,
  aes(
    x = dist.mean.bin,
    y = corr.mean.bin,
    colour = corr.grp
  )
) +
  geom_line(
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  geom_errorbar(
    aes(
      ymin = corr.mean.bin - corr.sd.bin,
      ymax = corr.mean.bin + corr.sd.bin
    ),
    width = 0.03,
    alpha = 0.06,
    linewidth = 0.4
  ) +
  scale_color_manual(
    values = c("#CC6666", "#9999CC")
  ) +
  labs(
    x = "Genomic distance (BP)",
    y = "Mean pairwise correlation"
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.text = element_text(size = 12),
    legend.position = "none"
  ) +
  facet_zoom(
    xlim = c(0, 2000),
    split = TRUE,
    zoom.size = 1.5
  )

# ------------------------------------------------------------
# 5. Save figure
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "Fig_S5b_correlation_decay.pdf"
  ),
  plot = p_s5b,
  width = 12,
  height = 5,
  units = "in"
)

message("Supplementary Figure S5b generated successfully.")

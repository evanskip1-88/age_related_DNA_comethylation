#!/usr/bin/env Rscript

# =============================================================================
# Supplementary Figure S5a, S5c and S5d
# Genomic distance and pairwise CpG co-methylation
#
# Input:
#   precorr_gs_qc_chr_all_Ages_vsOthers_long.cor
#   GS_sub_qc1.opi
#
# Outputs:
#   Fig_S5a_pairwise_distance_distribution.pdf
#   Fig_S5c_correlation_vs_distance.pdf
#   Fig_S5d_correlation_loss_vs_distance.pdf
#
# Notes:
#   CpG pair distances are measured in Mb.
#   Age1 reference correlations are compared with Age8.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggforce)
  library(ggpubr)
})

# -----------------------------------------------------------------------------
# 1. Paths
# -----------------------------------------------------------------------------

corr_file <- "/results/03_pairwise_correlations/precorr_gs_qc_chr_all_Ages_vsOthers_long.cor"
opi_file <- "/data/gsDNAm"
output_dir <- "figures"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 2. Load data
# -----------------------------------------------------------------------------

dt_corr <- fread(corr_file)
opi <- fread(opi_file)

required_corr <- c(
  "probe1", "probe2", "AgeBaseCorr",
  "AgeBase", "Age.grp", "corr"
)

required_opi <- c("V1", "V2", "V3")

stopifnot(
  all(required_corr %in% names(dt_corr)),
  all(required_opi %in% names(opi))
)

# -----------------------------------------------------------------------------
# 3. Prepare CpG annotations
# -----------------------------------------------------------------------------

probe_map <- unique(
  opi[
    !is.na(V1) & !is.na(V2) & !is.na(V3),
    .(
      probe = as.character(V2),
      chr = as.character(V1),
      bp = as.numeric(V3)
    )
  ],
  by = "probe"
)

map1 <- copy(probe_map)
setnames(
  map1,
  c("probe", "chr", "bp"),
  c("probe1", "probe1_chr", "probe1_bp")
)

map2 <- copy(probe_map)
setnames(
  map2,
  c("probe", "chr", "bp"),
  c("probe2", "probe2_chr", "probe2_bp")
)

# -----------------------------------------------------------------------------
# 4. Annotate both CpGs
# -----------------------------------------------------------------------------

dt <- merge(
  dt_corr,
  map1,
  by = "probe1",
  all = FALSE
)

dt <- merge(
  dt,
  map2,
  by = "probe2",
  all = FALSE
)

# Retain within-chromosome pairs for genomic distance analysis
dt <- dt[
  probe1_chr == probe2_chr &
    is.finite(probe1_bp) &
    is.finite(probe2_bp)
]

dt[, dist_diff := abs(probe1_bp - probe2_bp) / 1e6]

# Direction-independent pair identifier
dt[, probe_pair := paste(
  pmin(probe1, probe2),
  pmax(probe1, probe2),
  sep = "_"
)]

# -----------------------------------------------------------------------------
# 5. Figure S5a: Genomic distance distribution
# -----------------------------------------------------------------------------

dt_dist <- unique(
  dt[, .(probe_pair, dist_diff)],
  by = "probe_pair"
)

p_s5a <- ggplot(
  dt_dist,
  aes(x = dist_diff)
) +
  geom_histogram(
    binwidth = 2,
    colour = "black",
    fill = "lightblue"
  ) +
  labs(
    x = "Probe pairwise distance (Mb)",
    y = "Number of probe pairs"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 15)
  ) +
  facet_zoom(
    xlim = c(0, 50),
    split = TRUE,
    zoom.size = 1.5
  )

ggsave(
  filename = file.path(
    output_dir,
    "Fig_S5a_pairwise_distance_distribution.pdf"
  ),
  plot = p_s5a,
  width = 12,
  height = 6
)

# -----------------------------------------------------------------------------
# 6. Prepare Age1 versus Age8 comparison
# -----------------------------------------------------------------------------

dt_age1_age8 <- unique(
  dt[
    AgeBase == "Age1" &
      Age.grp == "Age8",
    .(
      probe_pair,
      probe1,
      probe2,
      probe1_chr,
      probe2_chr,
      probe1_bp,
      probe2_bp,
      AgeBaseCorr,
      corr,
      dist_diff
    )
  ],
  by = "probe_pair"
)

if (nrow(dt_age1_age8) == 0L) {
  stop("No Age1 versus Age8 comparison records found.")
}

# Original definition:
# Loss of absolute correlation from Age1 to Age8

dt_age1_age8[
  ,
  corrDiff := abs(AgeBaseCorr) - abs(corr)
]

dt_age1_age8[
  ,
  corrDirection := factor(
    ifelse(AgeBaseCorr > 0, "Positive", "Negative"),
    levels = c("Positive", "Negative")
  )
]

# -----------------------------------------------------------------------------
# 7. Common plotting function
# -----------------------------------------------------------------------------

plot_distance_relationship <- function(dat, y_variable, y_label) {

  ggplot(
    dat,
    aes(
      x = dist_diff,
      y = .data[[y_variable]],
      colour = corrDirection
    )
  ) +
    geom_point(
      alpha = 0.35,
      size = 0.6
    ) +
    geom_smooth(
      method = "lm",
      colour = "black",
      se = FALSE
    ) +
    stat_cor(
      aes(
        group = corrDirection,
        colour = corrDirection,
        label = after_stat(
          gsub("R", "italic(r)", r.label)
        )
      ),
      method = "pearson",
      parse = TRUE,
      label.x = 30,
      label.y = 0.05,
      size = 6
    ) +
    facet_wrap(~corrDirection) +
    theme_bw() +
    labs(
      x = "Distance between CpG sites (Mb)",
      y = y_label,
      colour = ""
    ) +
    theme(
      text = element_text(size = 18),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      )
    )
}

# -----------------------------------------------------------------------------
# 8. Figure S5c: Correlation versus genomic distance
# -----------------------------------------------------------------------------

p_s5c <- plot_distance_relationship(
  dat = dt_age1_age8,
  y_variable = "corr",
  y_label = "Pairwise correlation in Age8"
)

ggsave(
  filename = file.path(
    output_dir,
    "Fig_S5c_correlation_vs_distance.pdf"
  ),
  plot = p_s5c,
  width = 12,
  height = 6
)

# -----------------------------------------------------------------------------
# 9. Figure S5d: Correlation loss versus genomic distance
# -----------------------------------------------------------------------------

p_s5d <- plot_distance_relationship(
  dat = dt_age1_age8,
  y_variable = "corrDiff",
  y_label = "Correlation loss (Age1 vs Age8)"
) +
  theme(
    text = element_text(size = 16)
  )

ggsave(
  filename = file.path(
    output_dir,
    "Fig_S5d_correlation_loss_vs_distance.pdf"
  ),
  plot = p_s5d,
  width = 12,
  height = 6
)

message("Supplementary Figures S5a, S5c and S5d generated.")

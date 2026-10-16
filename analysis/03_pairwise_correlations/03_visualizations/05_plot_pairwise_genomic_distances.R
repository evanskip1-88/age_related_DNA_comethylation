#!/usr/bin/env Rscript

# =============================================================================
# Supplementary Figure S6
# Genomic distribution of co-methylated CpG pairs
#
# Inputs:
#   precorr_gs_qc_chr_all_Ages_vsOthers_long.cor
#   GS_sub_qc1.opi
#
# Outputs:
#   Fig_S6_pairwise_positions.pdf
#   Fig_S6_pairwise_distance_gradient.pdf
#   Fig_S6_pairwise_distance_categories.pdf
#
# Description:
#   Maps genomic positions of strongly co-methylated CpG pairs
#   identified using Age1 as the reference age group.
#
#   Three visualisations:
#     1. Genomic positions
#     2. Continuous genomic distance
#     3. Categorical genomic distance
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# 1. Configure paths
# -----------------------------------------------------------------------------

corr_file <- "/results/03_pairwise_correlations/precorr_gs_qc_chr_all_AgebaseAll_vsOthers.cor"
opi_file  <- "/data/gsDNAm"
output_dir <- "/figures"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# -----------------------------------------------------------------------------
# 2. Load data
# -----------------------------------------------------------------------------

dt_corr <- fread(corr_file)
opi <- fread(opi_file)

required_cols <- c(
  "probe1", "probe2", "AgeBase"
)

stopifnot(all(required_cols %in% names(dt_corr)))
stopifnot(all(c("V1", "V2", "V3") %in% names(opi)))

# -----------------------------------------------------------------------------
# 3. Retain Age1 reference pairs
# -----------------------------------------------------------------------------

dt_age1 <- dt_corr[
  AgeBase == "Age1",
  .(probe1, probe2)
]

# Each pair appears across multiple comparison age groups.
# Retain one record per unordered CpG pair.

dt_age1[, probe_low := pmin(probe1, probe2)]
dt_age1[, probe_high := pmax(probe1, probe2)]

dt_age1 <- unique(
  dt_age1[, .(
    probe1 = probe_low,
    probe2 = probe_high
  )]
)

# -----------------------------------------------------------------------------
# 4. Prepare probe annotation
# -----------------------------------------------------------------------------

probe_map <- unique(
  opi[
    !is.na(V1) &
      !is.na(V2) &
      !is.na(V3),
    .(
      probe = as.character(V2),
      chr = as.character(V1),
      bp = as.numeric(V3)
    )
  ],
  by = "probe"
)

# -----------------------------------------------------------------------------
# 5. Annotate both CpGs
# -----------------------------------------------------------------------------

probe1_map <- copy(probe_map)

setnames(
  probe1_map,
  c("probe", "chr", "bp"),
  c("probe1", "probe1_chr", "probe1_bp")
)

probe2_map <- copy(probe_map)

setnames(
  probe2_map,
  c("probe", "chr", "bp"),
  c("probe2", "probe2_chr", "probe2_bp")
)

plot_dt <- merge(
  dt_age1,
  probe1_map,
  by = "probe1",
  all = FALSE
)

plot_dt <- merge(
  plot_dt,
  probe2_map,
  by = "probe2",
  all = FALSE
)

# -----------------------------------------------------------------------------
# 6. Retain within-chromosome pairs
# -----------------------------------------------------------------------------

plot_dt <- plot_dt[
  probe1_chr == probe2_chr &
    is.finite(probe1_bp) &
    is.finite(probe2_bp)
]

plot_dt[, chr := probe1_chr]

# Standardise plotting orientation
plot_dt[, x_bp := pmin(probe1_bp, probe2_bp)]
plot_dt[, y_bp := pmax(probe1_bp, probe2_bp)]

# Calculate genomic distance
plot_dt[, distance_mb := (y_bp - x_bp) / 1e6]

# Order chromosomes
plot_dt[, chr := factor(
  chr,
  levels = as.character(1:22)
)]

if (nrow(plot_dt) == 0L) {
  stop("No annotated within-chromosome CpG pairs found.")
}

# -----------------------------------------------------------------------------
# 7. Common plot formatting
# -----------------------------------------------------------------------------

base_plot <- ggplot(
  plot_dt,
  aes(
    x = x_bp / 1e6,
    y = y_bp / 1e6
  )
) +
  facet_wrap(
    ~chr,
    scales = "free",
    labeller = labeller(
      chr = function(x) paste0("chr", x)
    )
  ) +
  theme_classic() +
  labs(
    x = "Lower CpG genomic position (Mb)",
    y = "Higher CpG genomic position (Mb)"
  ) +
  theme(
    text = element_text(size = 18),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

# -----------------------------------------------------------------------------
# 8. Plot 1: Genomic positions
# -----------------------------------------------------------------------------

p_positions <- base_plot +
  geom_point(
    size = 0.4,
    alpha = 0.4
  )

ggsave(
  filename = file.path(
    output_dir,
    "Fig_S6_pairwise_positions.pdf"
  ),
  plot = p_positions,
  width = 14,
  height = 12
)

# -----------------------------------------------------------------------------
# 9. Plot 2: Continuous genomic distance
# -----------------------------------------------------------------------------

# Log10 colour scale requires positive distances.
dt_gradient <- plot_dt[distance_mb > 0]

p_gradient <- ggplot(
  dt_gradient,
  aes(
    x = x_bp / 1e6,
    y = y_bp / 1e6,
    colour = distance_mb
  )
) +
  geom_point(
    size = 0.5,
    alpha = 0.5
  ) +
  facet_wrap(
    ~chr,
    scales = "free",
    labeller = labeller(
      chr = function(x) paste0("chr", x)
    )
  ) +
  scale_colour_viridis_c(
    name = "Distance (Mb)",
    trans = "log10"
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.3
  ) +
  theme_classic() +
  labs(
    x = "Lower CpG genomic position (Mb)",
    y = "Higher CpG genomic position (Mb)"
  ) +
  theme(
    text = element_text(size = 18),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(
  filename = file.path(
    output_dir,
    "Fig_S6_pairwise_distance_gradient.pdf"
  ),
  plot = p_gradient,
  width = 14,
  height = 12
)

# -----------------------------------------------------------------------------
# 10. Plot 3: Genomic distance categories
# -----------------------------------------------------------------------------

plot_dt[, distance_group := cut(
  distance_mb,
  breaks = c(-Inf, 0.1, 1, 10, Inf),
  labels = c(
    "<0.1 Mb",
    "0.1–1 Mb",
    "1–10 Mb",
    ">10 Mb"
  ),
  right = FALSE
)]

p_categories <- base_plot +
  geom_point(
    aes(colour = distance_group),
    size = 0.5,
    alpha = 0.5
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.3
  ) +
  labs(
    colour = "Genomic distance"
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(
        size = 4,
        alpha = 1
      )
    )
  )

ggsave(
  filename = file.path(
    output_dir,
    "Fig_S6_pairwise_distance_categories.pdf"
  ),
  plot = p_categories,
  width = 14,
  height = 12
)

message("Supplementary Figure S6 plots generated successfully.")

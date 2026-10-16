#!/usr/bin/env Rscript

# =============================================================================
# Figure 3: Overlap of co-methylated CpG pairs across age groups
#
# Input:
#   precorr_gs_qc_chr_all_AgebaseAll_vsOthers.cor
#
# Outputs:
#   Fig_3_positive_upset.pdf
#   Fig_3_negative_upset.pdf
#
# Description:
#   UpSet plots of CpG-pair overlap across eight age groups,
#   separated by the sign of their reference-group correlation.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(UpSetR)
})

# -----------------------------------------------------------------------------
# 1. Configure paths
# -----------------------------------------------------------------------------

input_file <- "/results/03_pairwise_correlations/precorr_gs_qc_chr_all_AgebaseAll_vsOthers.cor"
output_dir <- "figures"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 2. Load correlation data
# -----------------------------------------------------------------------------

dt_corr <- fread(input_file)

required_cols <- c(
  "probe1",
  "probe2",
  "AgeBase",
  "AgeBaseCorr",
  paste0("Age", 2:8)
)

stopifnot(all(required_cols %in% names(dt_corr)))

# -----------------------------------------------------------------------------
# 3. Reshape to long format
# -----------------------------------------------------------------------------

dt_long <- melt(
  dt_corr,
  id.vars = c(
    "probe1",
    "probe2",
    "AgeBase",
    "AgeBaseCorr"
  ),
  measure.vars = paste0("Age", 2:8),
  variable.name = "AgeGrp",
  value.name = "Corr"
)

# Exclude comparisons of an age group with itself
dt_long <- dt_long[AgeGrp != AgeBase]

# Direction-independent CpG pair identifiers
dt_long[, probe_pairs := paste(
  pmin(probe1, probe2),
  pmax(probe1, probe2),
  sep = "_v_"
)]

# -----------------------------------------------------------------------------
# 4. Function to generate an UpSet plot
# -----------------------------------------------------------------------------

plot_upset <- function(dt, correlation_sign, output_file) {

  if (correlation_sign == "positive") {

    dt_subset <- dt[AgeBaseCorr > 0]

  } else if (correlation_sign == "negative") {

    dt_subset <- dt[AgeBaseCorr < 0]

  } else {
    stop("correlation_sign must be positive or negative.")
  }

  # Unique CpG pairs identified in each reference age group
  probe_lists <- split(
    dt_subset$probe_pairs,
    dt_subset$AgeBase
  )

  probe_lists <- lapply(probe_lists, unique)

  age_groups <- paste0("Age", 1:8)

  # Ensure all age groups are represented
  for (age in age_groups) {
    if (!age %in% names(probe_lists)) {
      probe_lists[[age]] <- character(0)
    }
  }

  probe_lists <- probe_lists[age_groups]

  if (all(lengths(probe_lists) == 0L)) {
    stop("No CpG pairs available for ", correlation_sign, " correlations.")
  }

  # UpSetR draws using grid graphics
  pdf(
    file = output_file,
    width = 12,
    height = 7
  )

  upset(
    fromList(probe_lists),
    number.angles = 30,
    point.size = 3.0,
    sets = age_groups,
    mainbar.y.label = "Intersection Size",
    sets.x.label = "#Probe pairs Per Age Group",
    text.scale = rep(1.5, 6),
    order.by = "freq",
    keep.order = TRUE,
    show.numbers = "no"
  )

  dev.off()

  message("Saved: ", output_file)
}

# -----------------------------------------------------------------------------
# 5. Positive correlations
# -----------------------------------------------------------------------------

plot_upset(
  dt = dt_long,
  correlation_sign = "positive",
  output_file = file.path(
    output_dir,
    "Fig_3_positive_upset.pdf"
  )
)

# -----------------------------------------------------------------------------
# 6. Negative correlations
# -----------------------------------------------------------------------------

plot_upset(
  dt = dt_long,
  correlation_sign = "negative",
  output_file = file.path(
    output_dir,
    "Fig_3_negative_upset.pdf"
  )
)

message("Figure 3 UpSet plots generated successfully.")

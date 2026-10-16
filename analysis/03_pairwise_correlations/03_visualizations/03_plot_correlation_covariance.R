#!/usr/bin/env Rscript

# =============================================================================
# Supplementary Figure S2: Correlation and covariance trajectories
#
# Illustrative subset of Age1-reference co-methylated CpG pairs.
#
# For each correlation direction (positive and negative), select the
# 10 CpG pairs with the smallest heteroskedasticity p-values and plot
# their correlation and covariance across eight age groups.
#
# These are illustrative examples, not all co-methylated CpG pairs.
#
# Inputs:
#   precorr_gs_qc_chr_all_Age1_vsOthers.cor
#   chr_all_AgeBase1PosCorr_age_beta_est_resiq_sq.res
#   chr_all_AgeBase1NegCorr_age_beta_est_resiq_sq.res
#   chr_all_AgeBase1_covar.res
#
# Outputs:
#   Fig_S2c_positive_correlation_covariance.pdf
#   Fig_S2d_negative_correlation_covariance.pdf
#   Fig_S2_selected_pairs.tsv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# 1. Input/output paths
# -----------------------------------------------------------------------------

correlation_file <- "/results/03_pairwise_correlations/precorr_gs_qc_chr_all_AgebaseAll_vsOthers.cor"
positive_heteroskedasticity_file <- "/results/03_pairwise_correlations/chr_all_AgeBase1PosCorr_age_beta_est_resiq_sq.res"
negative_heteroskedasticity_file <- "/results/03_pairwise_correlations/chr_all_AgeBase1NegCorr_age_beta_est_resiq_sq.res"
covariance_file <- "/results/03_pairwise_correlations/chr_all_AgeBase1_covar.res"
output_dir <- "/figures"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_pairs <- 10L
age_levels <- paste0("Age", 1:8)

# -----------------------------------------------------------------------------
# 2. Load correlation data
# -----------------------------------------------------------------------------

dt_corr <- fread(correlation_file)

required_corr <- c("probe1", "probe2", age_levels)

if (!all(required_corr %in% names(dt_corr))) {
  stop("Correlation file is missing required columns.")
}

dt_corr[
  ,
  Probe_pair := paste(probe1, probe2, sep = "_")
]

# -----------------------------------------------------------------------------
# 3. Load heteroskedasticity results
# -----------------------------------------------------------------------------

hetero_columns <- c(
  "Probe_pair",
  "probe1",
  "probe2",
  "Chr",
  "AgeBase",
  "beta",
  "beta.se",
  "t_val",
  "p_val"
)

read_hetero <- function(path) {

  dt <- fread(path)

  if (ncol(dt) != length(hetero_columns)) {
    stop(
      "Unexpected heteroskedasticity file structure: ",
      path
    )
  }

  setnames(dt, hetero_columns)

  dt
}

dt_pos <- read_hetero(positive_heteroskedasticity_file)
dt_neg <- read_hetero(negative_heteroskedasticity_file)

# -----------------------------------------------------------------------------
# 4. Load covariance data
# -----------------------------------------------------------------------------

dt_cov <- fread(covariance_file)

required_cov <- c(
  "probe1",
  "probe2",
  "age_group",
  "covar"
)

if (!all(required_cov %in% names(dt_cov))) {
  stop("Covariance file is missing required columns.")
}

dt_cov[
  ,
  Probe_pair := paste(probe1, probe2, sep = "_")
]

dt_cov[
  ,
  age_group := as.character(age_group)
]

# -----------------------------------------------------------------------------
# 5. Select and prepare illustrative CpG pairs
# -----------------------------------------------------------------------------

prepare_pair_data <- function(hetero_dt, direction) {

  # Select the 10 smallest p-values, matching the original analysis
  selected <- hetero_dt[
    order(p_val)
  ][
    seq_len(min(n_pairs, .N))
  ]

  selected_ids <- selected$Probe_pair

  # Retain selected pairs from the Age1-reference correlation table
  corr_sub <- dt_corr[
    Probe_pair %in% selected_ids
  ]

  if (nrow(corr_sub) == 0L) {
    stop("No selected pairs matched correlation results: ", direction)
  }

  # Correlation trajectories
  corr_long <- melt(
    corr_sub,
    id.vars = "Probe_pair",
    measure.vars = age_levels,
    variable.name = "AgeGroup",
    value.name = "Value"
  )

  corr_long[
    ,
    Analysis := "Correlation"
  ]

  corr_long <- corr_long[
    ,
    .(Probe_pair, AgeGroup, Value, Analysis)
  ]

  # Covariance trajectories for the same CpG pairs
  cov_long <- dt_cov[
    Probe_pair %in% corr_sub$Probe_pair,
    .(
      Probe_pair,
      AgeGroup = age_group,
      Value = covar,
      Analysis = "Covariance"
    )
  ]

  # Combine
  combined <- rbindlist(
    list(corr_long, cov_long),
    use.names = TRUE
  )

  combined[
    ,
    AgeGroup := factor(
      AgeGroup,
      levels = age_levels
    )
  ]

  combined[
    ,
    Analysis := factor(
      Analysis,
      levels = c("Correlation", "Covariance")
    )
  ]

  combined[
    ,
    Direction := direction
  ]

  list(
    selected = selected,
    combined = combined
  )
}

# -----------------------------------------------------------------------------
# 6. Plot correlation and covariance trajectories
# -----------------------------------------------------------------------------

plot_pair_trajectories <- function(plot_data) {

  ggplot(
    plot_data,
    aes(
      x = AgeGroup,
      y = Value,
      group = Probe_pair,
      colour = Probe_pair
    )
  ) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.7) +
    facet_wrap(
      ~ Analysis,
      nrow = 2,
      scales = "free_y"
    ) +
    labs(
      x = "Age group",
      y = NULL,
      colour = "CpG pair"
    ) +
    theme_bw(base_size = 15) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      legend.position = "right",
      strip.text = element_text(face = "bold")
    )
}

# -----------------------------------------------------------------------------
# 7. Generate both figures
# -----------------------------------------------------------------------------

analysis_sets <- list(
  Positive = dt_pos,
  Negative = dt_neg
)

selected_tables <- list()

for (direction in names(analysis_sets)) {

  result <- prepare_pair_data(
    analysis_sets[[direction]],
    direction
  )

  p <- plot_pair_trajectories(
    result$combined
  )

  output_name <- if (direction == "Positive") {
    "Fig_S2c_positive_correlation_covariance.pdf"
  } else {
    "Fig_S2d_negative_correlation_covariance.pdf"
  }

  ggsave(
    filename = file.path(output_dir, output_name),
    plot = p,
    width = 11,
    height = 9
  )

  selected_tables[[direction]] <- copy(result$selected)[
    ,
    Direction := direction
  ]

  message(
    direction,
    ": selected ",
    uniqueN(result$combined$Probe_pair),
    " CpG pairs"
  )
}

# -----------------------------------------------------------------------------
# 8. Export selected pairs
# -----------------------------------------------------------------------------

fwrite(
  rbindlist(
    selected_tables,
    use.names = TRUE,
    fill = TRUE
  ),
  file.path(output_dir, "Fig_S2_selected_pairs.tsv"),
  sep = "\t"
)


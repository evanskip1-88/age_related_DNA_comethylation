#!/usr/bin/env Rscript

# =============================================================================
# Supplementary Figure S2a,b
# Q-Q plots of heteroskedasticity-test p-values
#
# Inputs:
#   Combined CpG-pair correlation results
#   Combined heteroskedasticity-test results
#
# Outputs:
#   Fig_S2a_positive_correlation_QQ.pdf
#   Fig_S2b_negative_correlation_QQ.pdf
#
# Note:
#   The original analysis used a custom qqPlot() function.
#   This version uses an independent implementation in base R.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# 1. Paths
# -----------------------------------------------------------------------------

cor_file <- "/results/03_pairwise_correlations/precorr_gs_qc_chr_all_AgebaseAll_vsOthers.cor"
heterosk_file <- "/results/03_pairwise_correlations/chr_all_AgeBaseAll_age_beta_est_resiq_sq.res"
output_dir <- "/figures"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 2. Load results
# -----------------------------------------------------------------------------

cor_res <- fread(cor_file)
het_res <- fread(heterosk_file)

stopifnot(
  all(c("probe1", "probe2", "AgeBaseCorr") %in% names(cor_res)),
  all(c("probe1", "probe2", "p_val") %in% names(het_res))
)

# -----------------------------------------------------------------------------
# 3. Construct direction-independent CpG-pair identifiers
# -----------------------------------------------------------------------------

make_pair_id <- function(probe1, probe2) {
  paste(
    pmin(probe1, probe2),
    pmax(probe1, probe2),
    sep = "_"
  )
}

cor_res[, probe_pair := make_pair_id(probe1, probe2)]
het_res[, probe_pair := make_pair_id(probe1, probe2)]

# -----------------------------------------------------------------------------
# 4. Separate positive and negative reference correlations
# -----------------------------------------------------------------------------

positive_pairs <- unique(
  cor_res[AgeBaseCorr > 0, probe_pair]
)

negative_pairs <- unique(
  cor_res[AgeBaseCorr < 0, probe_pair]
)

pos_het <- het_res[
  probe_pair %in% positive_pairs
]

neg_het <- het_res[
  probe_pair %in% negative_pairs
]

message("Positive correlation tests: ", nrow(pos_het))
message("Negative correlation tests: ", nrow(neg_het))

# -----------------------------------------------------------------------------
# 5. Q-Q plotting function
# -----------------------------------------------------------------------------

plot_qq <- function(p_values, output_file, plot_title) {

  p_values <- p_values[
    !is.na(p_values) &
      is.finite(p_values) &
      p_values >= 0 &
      p_values <= 1
  ]

  if (length(p_values) == 0L) {
    stop("No valid p-values available for: ", plot_title)
  }

  n <- length(p_values)

  # Sort observed p-values from smallest to largest
  observed <- -log10(sort(pmax(p_values, .Machine$double.xmin)))

  # Expected uniform p-value quantiles
  expected <- -log10(ppoints(n))

  # Both axes must run from smallest to largest
  observed <- rev(observed)
  expected <- rev(expected)

  pdf(
    output_file,
    width = 6,
    height = 6
  )

  plot(
    expected,
    observed,
    pch = 16,
    cex = 0.35,
    col = grDevices::adjustcolor("black", alpha.f = 0.35),
    xlab = expression(Expected~~-log[10](italic(p))),
    ylab = expression(Observed~~-log[10](italic(p))),
    main = plot_title,
    ylim = c(0, 50),
    xlim = c(0, max(expected))
  )

  abline(
    a = 0,
    b = 1,
    col = "red",
    lty = 2,
    lwd = 1.5
  )

  dev.off()

  message("Saved: ", output_file)
}

# -----------------------------------------------------------------------------
# 6. Supplementary Figure S2a: Positive correlations
# -----------------------------------------------------------------------------

plot_qq(
  p_values = pos_het$p_val,
  output_file = file.path(
    output_dir,
    "Fig_S2a_positive_correlation_QQ.pdf"
  ),
  plot_title = "Positive correlations"
)

# -----------------------------------------------------------------------------
# 7. Supplementary Figure S2b: Negative correlations
# -----------------------------------------------------------------------------

plot_qq(
  p_values = neg_het$p_val,
  output_file = file.path(
    output_dir,
    "Fig_S2b_negative_correlation_QQ.pdf"
  ),
  plot_title = "Negative correlations"
)

message("Supplementary Figures S2a and S2b generated.")

#!/usr/bin/env Rscript

# =============================================================================
# Pairwise CpG covariance analysis across eight age groups
#
# For heteroskedasticity-significant CpG pairs identified using Age1
# as the reference, calculate covariance separately in Age1–Age8
# and fit the original covariance ~ age_group ANOVA model.
#
# Inputs:
#   chr_all_AgeBase1_age_beta_est_resiq_sq.res
#   precorr_gs_qc_chr{chr}.txt
#   GS_18869_all.cov
#
# Output:
#   chr{chr}_AgeBase1_cov_anova_test.res
#
# Note:
#   The age-specific covariance output remains commented out,
#   consistent with the original analysis.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(broom)
})

cat("Analysis started:", as.character(Sys.time()), "\n\n")

# -----------------------------------------------------------------------------
# 1. Paths
# -----------------------------------------------------------------------------

dnam_path   <- "..."  # Individual level data
cov_inpath  <- "..."  # Individual level data
probe_inpath <- "/results/03_pairwise_correlations"
output_path <- "/results/03_pairwise_correlations"

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 2. Analysis parameters
# -----------------------------------------------------------------------------

Age <- 1

bonferroni_threshold <- 0.05 / 752722

# Load Age1-reference heteroskedasticity results
probe_pairs <- fread(
  file.path(
    probe_inpath,
    paste0(
      "chr_all_AgeBase",
      Age,
      "_age_beta_est_resiq_sq.res"
    )
  )
)

# Load participant covariates
all_covar <- fread(
  file.path(cov_inpath, "GS_18869_all.cov")
)

# Preserve original removal of the first covariate column
all_covar <- all_covar[, -1, with = FALSE]

# -----------------------------------------------------------------------------
# 3. Analyse each chromosome
# -----------------------------------------------------------------------------

for (chr in 1:22) {

  cat("Analysis for chromosome:", chr, "\n")

  # Select significant heteroskedasticity CpG pairs
  probe_pairs_sig <- probe_pairs[
    p_val < bonferroni_threshold &
    Chr == chr
  ]

  cat(
    "Significant heteroskedasticity probe pairs:",
    nrow(probe_pairs_sig),
    "\n"
  )

  if (nrow(probe_pairs_sig) == 0L) {
    next
  }

  # ---------------------------------------------------------------------------
  # Load residualised DNAm values
  # ---------------------------------------------------------------------------

  dnam <- fread(
    file.path(
      dnam_path,
      paste0("precorr_gs_qc_chr", chr, ".txt")
    )
  )

  # Merge methylation data with participant covariates
  m <- merge(
    all_covar,
    dnam,
    by = "IID_DNAm"
  )

  # ---------------------------------------------------------------------------
  # Assign age groups
  # ---------------------------------------------------------------------------

  m <- m %>%
    mutate(
      age_group = case_when(
        age < 36.0                 ~ 1,
        age >= 36.0 & age < 44.0   ~ 2,
        age >= 44.0 & age < 49.5   ~ 3,
        age >= 49.5 & age < 54.2   ~ 4,
        age >= 54.2 & age < 57.8   ~ 5,
        age >= 57.8 & age < 61.5   ~ 6,
        age >= 61.5 & age < 66.0   ~ 7,
        age >= 66.0                ~ 8
      )
    )

  # ---------------------------------------------------------------------------
  # Calculate covariance and fit ANOVA for each CpG pair
  # ---------------------------------------------------------------------------

  cat("Covariance ANOVA for chromosome", chr, "\n")

  cov_res_list <- vector("list", nrow(probe_pairs_sig))
  anova_covar_list <- vector("list", nrow(probe_pairs_sig))

  for (i in seq_len(nrow(probe_pairs_sig))) {

    probe1 <- probe_pairs_sig$probe1[i]
    probe2 <- probe_pairs_sig$probe2[i]

    # Extract the two CpGs and participant age groups
    m_sub <- m[
      ,
      c(probe1, probe2, "age_group", "IID_DNAm"),
      with = FALSE
    ]

    m_sub[, Chr := chr]

    # Calculate covariance within each age group
    covs <- m_sub[
      ,
      .(
        covar = cov(
          .SD[[1]],
          .SD[[2]]
        )
      ),
      by = age_group,
      .SDcols = c(probe1, probe2)
    ]

    covs[, `:=`(
      probe1 = probe1,
      probe2 = probe2
    )]

    cov_res_list[[i]] <- covs

    # Original ANOVA model: numeric age-group predictor
    anova_model <- aov(
      covar ~ age_group,
      data = covs
    )

    tidy_anova <- broom::tidy(anova_model)

    anova_covar_list[[i]] <- data.table(
      probe1 = probe1,
      probe2 = probe2,
      Chr = chr,
      term = tidy_anova$term[1],
      anova.stats = tidy_anova$statistic[1],
      p.value = tidy_anova$p.value[1]
    )
  }

  # ---------------------------------------------------------------------------
  # Combine chromosome results
  # ---------------------------------------------------------------------------

  cov_all_probe <- rbindlist(
    cov_res_list,
    use.names = TRUE
  )

  cov_all_probe[
    ,
    age_group := paste0("Age", age_group)
  ]

  anova_all_probe <- rbindlist(
    anova_covar_list,
    use.names = TRUE
  )

  # ---------------------------------------------------------------------------
  # Save results
  # ---------------------------------------------------------------------------

  # Age-specific covariance output intentionally disabled.
  #
  # fwrite(
  #   cov_all_probe,
  #   file.path(
  #     output_path,
  #     paste0("chr", chr, "_AgeBase", Age, "_covar.res")
  #   ),
  #   quote = FALSE,
  #   sep = "\t"
  # )

  # Save ANOVA results
  fwrite(
    anova_all_probe,
    file.path(
      output_path,
      paste0(
        "chr",
        chr,
        "_AgeBase",
        Age,
        "_cov_anova_test.res"
      )
    ),
    quote = FALSE,
    sep = "\t"
  )

  cat("Completed chromosome:", chr, "\n\n")
}

cat("Analysis complete!\n")
cat("Analysis ended:", as.character(Sys.time()), "\n")

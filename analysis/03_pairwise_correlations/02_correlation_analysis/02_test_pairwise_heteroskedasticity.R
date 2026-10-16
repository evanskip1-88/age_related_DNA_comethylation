# ============================================================
# Age-related heteroskedasticity of CpG-pair relationships
#
# Primary model:
#   CpG1 ~ CpG2 + age_group
#
# Secondary model:
#   normalised_squared_residuals ~ age_group
#
# age_group is numeric (1-8), preserving the original
# linear trend analysis.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 6L) {
  stop(
    "Usage: Rscript 02_test_pairwise_heteroskedasticity.R ",
    "<chr> <AgeGrp> <dnam_path> <cov_inpath> ",
    "<probe_inpath> <output_path>"
  )
}

chr <- as.integer(args[1])
Age <- as.integer(args[2])

dnam_path <- args[3]
cov_inpath <- args[4]
probe_inpath <- args[5]
output_path <- args[6]

if (is.na(chr) || !(chr %in% 1:22)) {
  stop("Chromosome must be between 1 and 22.")
}

if (is.na(Age) || !(Age %in% 1:8)) {
  stop("Reference age group must be between 1 and 8.")
}

dir.create(
  output_path,
  recursive = TRUE,
  showWarnings = FALSE
)

start_time <- Sys.time()
cat("Analysis started:", as.character(start_time), "\n")

# ------------------------------------------------------------
# 2. Input files
# ------------------------------------------------------------

pair_file <- file.path(
  probe_inpath,
  sprintf("precorr_gs_qc_chr%d_Age%d.cor", chr, Age)
)

dnam_file <- file.path(
  dnam_path,
  sprintf("precorr_gs_qc_chr%d.txt", chr)
)

cov_file <- file.path(
  cov_inpath,
  "GS_18869_all.cov"
)

for (f in c(pair_file, dnam_file, cov_file)) {
  if (!file.exists(f)) {
    stop("Missing input file: ", f)
  }
}

# ------------------------------------------------------------
# 3. Load CpG pairs and residualised DNAm
# ------------------------------------------------------------

probe_pairs <- fread(pair_file)
dnam <- fread(dnam_file)
covar <- fread(cov_file)

required_pair_cols <- c(
  "probe1", "probe2", "Chr", "AgeBase"
)

if (!all(required_pair_cols %in% names(probe_pairs))) {
  stop("Input pair file is missing required columns.")
}

if (!"IID_DNAm" %in% names(dnam) ||
    !"IID_DNAm" %in% names(covar) ||
    !"age" %in% names(covar)) {
  stop("Required sample identifiers or age column missing.")
}

# ------------------------------------------------------------
# 4. Merge methylation and covariates
# ------------------------------------------------------------

m <- merge(
  covar,
  dnam,
  by = "IID_DNAm"
)

# ------------------------------------------------------------
# 5. Assign age groups
# ------------------------------------------------------------

age_breaks <- c(
  17.9, 36.0, 44.0, 49.5,
  54.2, 57.8, 61.5, 66.0, 94.0
)

m[, age_group := as.integer(
  cut(
    age,
    breaks = age_breaks,
    right = FALSE,
    include.lowest = TRUE,
    labels = FALSE
  )
)]

m <- m[!is.na(age_group)]

cat("Samples retained:", nrow(m), "\n")
cat("CpG pairs:", nrow(probe_pairs), "\n")

# ------------------------------------------------------------
# 6. Heteroskedasticity analysis
# ------------------------------------------------------------

results <- vector(
  "list",
  nrow(probe_pairs)
)

for (i in seq_len(nrow(probe_pairs))) {

  probe1 <- probe_pairs$probe1[i]
  probe2 <- probe_pairs$probe2[i]

  if (!(probe1 %in% names(m)) ||
      !(probe2 %in% names(m))) {
    next
  }

  # Work with the two CpGs and numeric age group
  df <- data.frame(
    y = m[[probe1]],
    x = m[[probe2]],
    age_group = m$age_group
  )

  # Complete observations for this pair
  df <- df[complete.cases(df), , drop = FALSE]

  if (nrow(df) < 4L) {
    next
  }

  # Primary model
  fit1 <- lm(
    y ~ x + age_group,
    data = df
  )

  residual_sq <- residuals(fit1)^2
  mean_residual_sq <- mean(residual_sq)

  if (!is.finite(mean_residual_sq) ||
      mean_residual_sq <= 0) {
    next
  }

  # Normalised squared residuals
  df$norm_resid <- residual_sq / mean_residual_sq

  # Secondary model
  fit2 <- lm(
    norm_resid ~ age_group,
    data = df
  )

  coef_table <- summary(fit2)$coefficients

  if (!"age_group" %in% rownames(coef_table)) {
    next
  }

  results[[i]] <- data.table(
    probe1 = probe1,
    probe2 = probe2,
    Chr = probe_pairs$Chr[i],
    AgeBase = probe_pairs$AgeBase[i],
    beta = coef_table["age_group", "Estimate"],
    beta.se = coef_table["age_group", "Std. Error"],
    t_val = coef_table["age_group", "t value"],
    p_val = coef_table["age_group", "Pr(>|t|)"]
  )

  if (i %% 1000L == 0L) {
    cat(
      "Processed", i, "of",
      nrow(probe_pairs), "pairs\n"
    )
  }
}

# ------------------------------------------------------------
# 7. Save results
# ------------------------------------------------------------

results <- Filter(Negate(is.null), results)

if (length(results) == 0L) {
  warning("No valid CpG-pair regression results.")
  result_dt <- data.table(
    probe1 = character(),
    probe2 = character(),
    Chr = integer(),
    AgeBase = character(),
    beta = numeric(),
    beta.se = numeric(),
    t_val = numeric(),
    p_val = numeric()
  )
} else {
  result_dt <- rbindlist(results)
}

output_file <- file.path(
  output_path,
  sprintf(
    "chr%d_AgeBase%d_age_beta_est_resiq_sq.txt",
    chr,
    Age
  )
)

fwrite(
  result_dt,
  output_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("Results saved:", output_file, "\n")
cat("Analysis completed:", as.character(Sys.time()), "\n")
cat(
  "Elapsed minutes:",
  round(
    as.numeric(
      difftime(Sys.time(), start_time, units = "mins")
    ),
    2
  ),
  "\n"
)

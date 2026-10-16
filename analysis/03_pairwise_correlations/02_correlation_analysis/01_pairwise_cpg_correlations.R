# ============================================================
# Pairwise CpG correlation analysis across age groups
#
# Identify CpG pairs with |r| > 0.8 in a reference age group
# and calculate their correlations across all eight age groups.
#
# Input:
#   Residualised methylation: precorr_gs_qc_chr*.txt
#   Covariates: GS_18869_all.cov
#
# Output:
#   One .cor file per chromosome and reference age group
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(HiClimR)
})

# ------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------

dnam_path   <- "..."
cov_inpath  <- "..."
output_path <- "..."

# Values supplied by SLURM or set manually
chr <- 1L
age_groups_to_run <- 1:8

cor_cutoff <- 0.8

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

start_time <- Sys.time()
cat("Analysis started:", as.character(start_time), "\n")

# ------------------------------------------------------------
# 2. Define age groups
# ------------------------------------------------------------

age_groups <- matrix(
  c(
    17.9, 36.0,
    36.0, 44.0,
    44.0, 49.5,
    49.5, 54.2,
    54.2, 57.8,
    57.8, 61.5,
    61.5, 66.0,
    66.0, 94.0
  ),
  ncol = 2,
  byrow = TRUE
)

age_labels <- paste0("Age", seq_len(nrow(age_groups)))

# ------------------------------------------------------------
# 3. Load residualised methylation and covariates
# ------------------------------------------------------------

input_file <- file.path(
  dnam_path,
  paste0("precorr_gs_qc_chr", chr, ".txt")
)

cov_file <- file.path(
  cov_inpath,
  "GS_18869_all.cov"
)

stopifnot(
  file.exists(input_file),
  file.exists(cov_file)
)

dnam <- fread(input_file)
covar <- fread(cov_file)

stopifnot(
  "IID_DNAm" %in% names(dnam),
  "IID_DNAm" %in% names(covar),
  "age" %in% names(covar)
)

# Merge covariates and residualised methylation
merged <- merge(
  covar,
  dnam,
  by = "IID_DNAm"
)

# Identify CpG columns explicitly
probe_ids <- setdiff(
  names(dnam),
  "IID_DNAm"
)

cat("Samples:", nrow(merged), "\n")
cat("CpGs:", length(probe_ids), "\n")

# ------------------------------------------------------------
# 4. Split residualised methylation by age group
# ------------------------------------------------------------

age_data <- vector("list", nrow(age_groups))

for (i in seq_len(nrow(age_groups))) {

  lower <- age_groups[i, 1]
  upper <- age_groups[i, 2]

  samples <- merged[
    age >= lower & age < upper
  ]

  age_data[[i]] <- as.matrix(
    samples[, ..probe_ids]
  )

  storage.mode(age_data[[i]]) <- "numeric"

  cat(
    age_labels[i], ":",
    nrow(age_data[[i]]), "samples\n"
  )
}

# ------------------------------------------------------------
# 5. Identify highly correlated CpG pairs
# ------------------------------------------------------------

for (ref_age in age_groups_to_run) {

  cat(
    "\nReference age group:",
    age_labels[ref_age], "\n"
  )

  ref_data <- age_data[[ref_age]]

  if (nrow(ref_data) < 3L) {
    warning("Insufficient samples: ", age_labels[ref_age])
    next
  }

  # Calculate upper-triangular correlation matrix
  cor_mat <- HiClimR::fastCor(
    ref_data,
    nSplit = 10,
    upperTri = TRUE,
    optBLAS = TRUE,
    verbose = TRUE
  )

  # Retain strongly correlated pairs
  pair_indices <- which(
    abs(cor_mat) > cor_cutoff &
      upper.tri(cor_mat, diag = FALSE),
    arr.ind = TRUE
  )

  if (nrow(pair_indices) == 0L) {
    cat("No strongly correlated pairs found.\n")
    next
  }

  high_cor <- data.table(
    probe1 = colnames(cor_mat)[pair_indices[, 1]],
    probe2 = colnames(cor_mat)[pair_indices[, 2]]
  )

  cat(
    "Selected CpG pairs:",
    nrow(high_cor), "\n"
  )

  # ----------------------------------------------------------
  # 6. Calculate pairwise correlations in every age group
  # ----------------------------------------------------------

  probe1_idx <- pair_indices[, 1]
  probe2_idx <- pair_indices[, 2]

  for (age_i in seq_along(age_data)) {

    x <- age_data[[age_i]][, probe1_idx, drop = FALSE]
    y <- age_data[[age_i]][, probe2_idx, drop = FALSE]

    # Pairwise Pearson correlation for each selected pair
    pair_cor <- vapply(
      seq_along(probe1_idx),
      function(j) {
        cor(
          x[, j],
          y[, j],
          use = "everything",
          method = "pearson"
        )
      },
      numeric(1)
    )

    high_cor[
      ,
      (age_labels[age_i]) := pair_cor
    ]
  }

  # ----------------------------------------------------------
  # 7. Save results
  # ----------------------------------------------------------

  high_cor[, Chr := chr]
  high_cor[, AgeBase := age_labels[ref_age]]

  output_file <- file.path(
    output_path,
    paste0(
      "precorr_gs_qc_chr",
      chr,
      "_Age",
      ref_age,
      ".cor"
    )
  )

  fwrite(
    high_cor,
    output_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  cat("Saved:", output_file, "\n")
}

cat(
  "\nAnalysis completed:",
  as.character(Sys.time()), "\n"
)

cat(
  "Elapsed minutes:",
  round(
    as.numeric(
      difftime(
        Sys.time(),
        start_time,
        units = "mins"
      )
    ),
    2
  ),
  "\n"
)

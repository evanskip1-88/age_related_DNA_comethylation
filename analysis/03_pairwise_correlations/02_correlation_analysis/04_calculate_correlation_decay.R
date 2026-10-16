#!/usr/bin/env Rscript

# =============================================================================
# Correlation decay with genomic distance
#
# Input:
#   Residualised DNAm matrix for one chromosome and age group
#   GS_sub_qc1.opi (CpG annotation)
#
# Output:
#   Chromosome- and age-specific correlation decay summary
#
# Method:
#   1. Identify CpG pairs separated by < 10 kb
#   2. Calculate Pearson correlations
#   3. Divide distances into 200 quantile-based bins
#   4. Summarise positive and negative correlations separately
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Hmisc)
})

# -----------------------------------------------------------------------------
# 1. Command-line arguments
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5L) {
  stop(
    "Usage: Rscript 04_calculate_correlation_decay.R ",
    "<dnam_file> <opi_file> <chr> <age_group> <output_dir>"
  )
}

dnam_file  <- args[1]
opi_file   <- args[2]
chr        <- args[3]
age_group  <- as.integer(args[4])
output_dir <- args[5]

if (is.na(age_group) || !age_group %in% 1:8) {
  stop("age_group must be an integer from 1 to 8.")
}

if (!file.exists(dnam_file)) stop("Missing DNAm file: ", dnam_file)
if (!file.exists(opi_file)) stop("Missing annotation file: ", opi_file)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

window_size <- 10000
n_bins <- 200L

start_time <- Sys.time()

message("Analysis started: ", start_time)
message("Chromosome: ", chr)
message("Age group: ", age_group)

# -----------------------------------------------------------------------------
# 2. Load data
# -----------------------------------------------------------------------------

dnam <- fread(dnam_file)
opi  <- fread(opi_file, header = FALSE)

# Original input format:
# First 13 columns contain sample metadata
# Columns 14 onward contain residualised DNAm values

probe_names <- names(dnam)[14:ncol(dnam)]

# Annotation:
# V2 = CpG probe ID
# V3 = genomic position

opi <- opi[
  V2 %in% probe_names &
    !is.na(V3),
  .(probe = V2, position = as.numeric(V3))
]

opi <- unique(opi, by = "probe")
setorder(opi, position, probe)

if (nrow(opi) < 2L) {
  stop("Fewer than two annotated CpGs found.")
}

# Keep the DNAm columns in exactly the same order as annotation

dnam_mat <- as.matrix(dnam[, ..opi$probe])
storage.mode(dnam_mat) <- "double"

positions <- opi$position
probes <- opi$probe

message("Number of CpGs: ", length(probes))

# -----------------------------------------------------------------------------
# 3. Identify CpG pairs within 10 kb
# -----------------------------------------------------------------------------

# CpGs are sorted by genomic position.
# For each CpG, find subsequent CpGs within the window.
#
# This avoids constructing the full chromosome-wide distance matrix.

pair_i <- vector("list", length(probes))
pair_j <- vector("list", length(probes))

for (i in seq_len(length(probes) - 1L)) {

  j <- i + 1L

  while (
    j <= length(probes) &&
    (positions[j] - positions[i]) < window_size
  ) {
    pair_i[[i]] <- c(pair_i[[i]], i)
    pair_j[[i]] <- c(pair_j[[i]], j)

    j <- j + 1L
  }
}

idx_i <- unlist(pair_i, use.names = FALSE)
idx_j <- unlist(pair_j, use.names = FALSE)

rm(pair_i, pair_j)
gc()

if (length(idx_i) == 0L) {
  stop("No CpG pairs found within 10 kb.")
}

pair_dist <- positions[idx_j] - positions[idx_i]

message("Eligible CpG pairs: ", length(idx_i))

# -----------------------------------------------------------------------------
# 4. Calculate pairwise Pearson correlations
# -----------------------------------------------------------------------------

# Calculate correlations in blocks to limit intermediate memory use.

block_size <- 10000L
n_pairs <- length(idx_i)

pair_corr <- numeric(n_pairs)

for (start in seq.int(1L, n_pairs, by = block_size)) {

  end <- min(start + block_size - 1L, n_pairs)

  for (k in start:end) {

    pair_corr[k] <- suppressWarnings(
      cor(
        dnam_mat[, idx_i[k]],
        dnam_mat[, idx_j[k]],
        use = "everything",
        method = "pearson"
      )
    )
  }
}

# -----------------------------------------------------------------------------
# 5. Construct correlation-distance table
# -----------------------------------------------------------------------------

dt_pairs <- data.table(
  probe1 = probes[idx_i],
  probe2 = probes[idx_j],
  dist = pair_dist,
  corr = pair_corr
)

dt_pairs <- dt_pairs[
  is.finite(dist) &
    is.finite(corr)
]

# -----------------------------------------------------------------------------
# 6. Divide distances into 200 quantile-based bins
# -----------------------------------------------------------------------------

# Preserve the original Hmisc::cut2 approach.
# Bin membership is determined using the pooled distances,
# before separating positive and negative correlations.

dt_pairs[, dist.bin := Hmisc::cut2(
  dist,
  g = n_bins
)]

# -----------------------------------------------------------------------------
# 7. Summarise positive and negative correlations
# -----------------------------------------------------------------------------

dt_pairs[, corr.grp := fifelse(corr < 0, "-ve", "+ve")]

dt_summary <- dt_pairs[
  ,
  .(
    dist.mean.bin = mean(dist),
    corr.mean.bin = mean(corr),
    corr.median.bin = median(corr),
    corr.sd.bin = sd(corr)
  ),
  by = .(dist.bin, corr.grp)
]

dt_summary[, age.grp := paste0("Age", age_group)]
dt_summary[, Chr := chr]

setcolorder(
  dt_summary,
  c(
    "dist.mean.bin",
    "corr.mean.bin",
    "corr.median.bin",
    "corr.sd.bin",
    "corr.grp",
    "age.grp",
    "Chr"
  )
)

setorder(dt_summary, corr.grp, dist.mean.bin)

# -----------------------------------------------------------------------------
# 8. Save results
# -----------------------------------------------------------------------------

output_file <- file.path(
  output_dir,
  paste0(
    "precorr_gs_qc_chr",
    chr,
    "_Age",
    age_group,
    "_decay.cor"
  )
)

fwrite(
  dt_summary,
  output_file,
  sep = "\t",
  quote = FALSE
)

message("Output: ", output_file)
message("Analysis completed: ", Sys.time())
message("Elapsed time: ", Sys.time() - start_time)

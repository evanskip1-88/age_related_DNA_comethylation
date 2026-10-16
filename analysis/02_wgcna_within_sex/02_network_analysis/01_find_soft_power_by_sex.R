# ============================================================
# Identify WGCNA soft-thresholding power within sex
# Reference network: Age group 1
# ============================================================

suppressPackageStartupMessages({
  library(WGCNA)
  library(data.table)
})

# Configure directories
dnam_inpath <- "..."
output_path <- "figures"

# Select sex
sex_id <- "F"  # "M" or "F"

# Load Age1 DNAm residuals
ref_fname <- file.path(
  dnam_inpath,
  paste0("precorr_gs_qc_chr_all_age_grp1_", sex_id, ".txt")
)

cat("Loading reference DNAm data:", ref_fname, "\n")

ref_dnam <- fread(ref_fname)

# Remove covariate columns (columns 2–13)
covs <- names(ref_dnam)[2:13]

ref_dnam_sub <- ref_dnam[, !..covs]

# Remove sample ID column and convert CpGs to matrix
ref_dnam_sub_mat <- as.matrix(
  ref_dnam_sub[, -1, with = FALSE]
)

rownames(ref_dnam_sub_mat) <- ref_dnam$IID_DNAm

cat("Samples:", nrow(ref_dnam_sub_mat), "\n")
cat("CpGs:", ncol(ref_dnam_sub_mat), "\n")

# Candidate soft-threshold powers
powers <- c(1:10, seq(12, 20, by = 2))

# Estimate scale-free topology fit
sft <- pickSoftThreshold(
  ref_dnam_sub_mat,
  dataIsExpr = TRUE,
  weights = NULL,
  RsquaredCut = 0.85,
  powerVector = powers,
  removeFirst = FALSE,
  nBreaks = 10,
  blockSize = NULL,
  corFnc = cor,
  corOptions = list(use = "p"),
  networkType = "signed",
  moreNetworkConcepts = FALSE,
  gcInterval = NULL,
  verbose = 5,
  indent = 0
)

# Plot scale-free topology fit
png(
  file.path(output_path, paste0("scale_free_fit_plot_", sex_id, ".png")),
  width = 800,
  height = 600
)

plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit, signed R^2",
  type = "n"
)

text(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = sft$fitIndices[, 1],
  cex = 0.9,
  col = "red"
)

abline(h = 0.8, col = "blue")

dev.off()

cat("Soft-threshold analysis completed for sex:", sex_id, "\n")

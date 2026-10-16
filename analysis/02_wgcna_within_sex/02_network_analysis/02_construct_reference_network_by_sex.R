# ============================================================
# Construct sex-specific WGCNA reference network
# Reference: Age group 1
# Soft-threshold power: 6 (empirically selected from females)
# ============================================================

suppressPackageStartupMessages({
  library(WGCNA)
  library(data.table)
})

# Configure directories
dnam_inpath <- "..."
output_path <- "/results/02_wgcna_within_sex"

# Reference network settings
sex_id <- "M"  # "M" or "F"
pickPower <- 6
blockSize <- 30000
nClusters <- 1881

deepsplit <- 3
cutHeight <- 0.10
min_mod <- 20

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

# Convert CpG data to matrix
ref_dnam_sub_mat <- as.matrix(
  ref_dnam_sub[, -1, with = FALSE]
)

rownames(ref_dnam_sub_mat) <- ref_dnam$IID_DNAm

cat("Samples:", nrow(ref_dnam_sub_mat), "\n")
cat("CpGs:", ncol(ref_dnam_sub_mat), "\n")
cat("Soft-threshold power:", pickPower, "\n")
cat("Maximum block size:", blockSize, "\n")
cat("Preclustering centres:", nClusters, "\n")

# Construct WGCNA reference network
ref_network <- blockwiseModules(
  datExpr = ref_dnam_sub_mat,
  power = pickPower,
  TOMType = "signed",
  minModuleSize = min_mod,
  networkType = "signed",

  # Block settings
  blocks = NULL,
  maxBlockSize = blockSize,
  blockSizePenaltyPower = 5,
  nPreclusteringCenters = nClusters,

  # Correlation
  corType = "bicor",

  # Module detection
  deepSplit = deepsplit,
  pamRespectsDendro = FALSE,
  useBlocks = TRUE,

  # Module refinement
  reassignThreshold = 0.05,
  mergeCutHeight = cutHeight,

  # Output settings
  numericLabels = TRUE,
  nThreads = 20,
  verbose = 5
)

# Summarise module assignments
print(table(ref_network$colors))

# Extract module labels
df_module_labs <- data.table(
  probe_id = names(ref_network$colors),
  moduleNumID = ref_network$colors,
  moduleColors = labels2colors(ref_network$colors),
  mergeCutHeight = cutHeight,
  deepSplit = deepsplit,
  power = pickPower
)

# Save module assignments
out_file <- file.path(
  output_path,
  paste0(
    "precorr_gs_qc_split", deepsplit,
    "_cut", cutHeight,
    "_", sex_id, ".modLabs"
  )
)

fwrite(
  df_module_labs,
  out_file,
  sep = "\t",
  quote = FALSE
)

cat("Reference network construction completed for sex:", sex_id, "\n")

# ==============================================================================
# Construct the reference DNA co-methylation network
#
# WGCNA is applied to residualised DNA methylation data from the youngest
# age group (Age group 1), which serves as the reference network.
#
# The large DNAm matrix is read in row-wise subsets to reduce memory
# requirements before network construction.
# ==============================================================================

suppressMessages(library(WGCNA))
suppressMessages(library(data.table))


# ------------------------------------------------------------------------------
# WGCNA parameters
# ------------------------------------------------------------------------------

ref_age <- "age1"

soft_power      <- 6
max_block_size  <- 30000
n_precluster    <- 1881
min_module_size <- 20
deep_split      <- 3
merge_cut_height <- 0.10


cat("Reference age group:", ref_age, "\n")
cat("Soft-thresholding power:", soft_power, "\n")
cat("Maximum block size:", max_block_size, "\n")
cat("Pre-clustering centres:", n_precluster, "\n")
cat("Minimum module size:", min_module_size, "\n")
cat("deepSplit:", deep_split, "\n")
cat("Merge cut height:", merge_cut_height, "\n\n")


# ------------------------------------------------------------------------------
# User-defined paths
# ------------------------------------------------------------------------------

dnam_inpath <- "..."
output_path <- "..."


# ------------------------------------------------------------------------------
# Reference DNAm dataset
# ------------------------------------------------------------------------------

ref_fname <- file.path(
    dnam_inpath,
    "precorr_gs_qc_chr_all_age1_roundoff.txt"
)


# ------------------------------------------------------------------------------
# Load reference DNAm data
#
# The matrix is read in three row-wise subsets to reduce memory requirements
# during file import.
# ------------------------------------------------------------------------------

cat("Loading reference DNAm data...\n")


# Read column names
dt_header <- fread(
    ref_fname,
    nrows = 5
)


# Samples 1-300
ref_data_sub1 <- fread(
    ref_fname,
    skip = 0,
    nrows = 300
)


# Samples 301-600
ref_data_sub2 <- fread(
    ref_fname,
    skip = 300,
    nrows = 300
)

setnames(
    ref_data_sub2,
    colnames(dt_header)
)


# Remaining samples
ref_data_sub3 <- fread(
    ref_fname,
    skip = 600
)

setnames(
    ref_data_sub3,
    colnames(dt_header)
)


# ------------------------------------------------------------------------------
# Recombine reference samples
# ------------------------------------------------------------------------------

ref_data <- rbind(
    ref_data_sub1,
    ref_data_sub2,
    ref_data_sub3
)


# First four columns contain sample/covariate information
ref_dnam <- data.matrix(
    ref_data[, -(1:4), with = FALSE]
)

rownames(ref_dnam) <- ref_data$IID_DNAm


cat(
    "Reference DNAm matrix:",
    nrow(ref_dnam), "samples x",
    ncol(ref_dnam), "CpGs\n\n"
)


# ------------------------------------------------------------------------------
# Construct WGCNA reference network
# ------------------------------------------------------------------------------

cat("Running WGCNA blockwiseModules()...\n\n")


ref_network <- blockwiseModules(
    datExpr = ref_dnam,

    # Network construction
    power = soft_power,
    networkType = "signed",
    TOMType = "signed",
    corType = "bicor",

    # Blocking
    blocks = NULL,
    maxBlockSize = max_block_size,
    blockSizePenaltyPower = 5,
    nPreclusteringCenters = n_precluster,

    # Module detection
    minModuleSize = min_module_size,
    deepSplit = deep_split,
    pamRespectsDendro = FALSE,

    # Module adjustment
    reassignThreshold = 0,
    mergeCutHeight = merge_cut_height,

    # Output/computation
    numericLabels = TRUE,
    nThreads = 10,
    verbose = 5
)


# ------------------------------------------------------------------------------
# Summarise module assignments
# ------------------------------------------------------------------------------

cat("\nModule sizes:\n")
print(table(ref_network$colors))


module_assignments <- data.table(
    probe_id = names(ref_network$colors),
    moduleNumID = ref_network$colors,
    moduleColor = labels2colors(ref_network$colors),
    mergeCutHeight = merge_cut_height,
    deepSplit = deep_split,
    power = soft_power
)


# ------------------------------------------------------------------------------
# Save module assignments
# ------------------------------------------------------------------------------

output_file <- file.path(
    output_path,
    paste0(
        "precorr_gs_qc_split",
        deep_split,
        "_cut",
        merge_cut_height,
        "_allsex.modLabs"
    )
)


fwrite(
    module_assignments,
    output_file,
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    col.names = TRUE
)


cat("\nReference network construction complete.\n")

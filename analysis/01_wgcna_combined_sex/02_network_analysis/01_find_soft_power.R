# ==============================================================================
# Select WGCNA soft-thresholding power
#
# The youngest age group (Age group 1 <= 36) was used as the reference dataset for
# constructing the co-methylation network. This script evaluates candidate
# soft-thresholding powers using WGCNA::pickSoftThreshold().
#
# The reference DNAm matrix is read in row-wise subsets to reduce memory
# requirements before being recombined for soft-threshold selection.
# ==============================================================================

suppressMessages(library(data.table))
suppressMessages(library(WGCNA))


# ------------------------------------------------------------------------------
# User-defined paths
# ------------------------------------------------------------------------------

dnam_inpath <- "..."
output_path <- "..."


# ------------------------------------------------------------------------------
# Input file: reference age group
# ------------------------------------------------------------------------------

ref_fname <- file.path(
    dnam_inpath,
    "precorr_gs_qc_chr_all_age1_roundoff.txt"
)


# Candidate soft-thresholding powers
powers <- c(
    seq(1, 10, by = 1),
    seq(12, 20, by = 2)
)


# ------------------------------------------------------------------------------
# Load reference DNAm data
#
# The large reference matrix is read in three row-wise subsets to reduce
# memory requirements during file import.
# ------------------------------------------------------------------------------

cat("Loading reference DNAm data...\n")


# Read header
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
    "Reference dataset loaded:",
    nrow(ref_dnam), "samples and",
    ncol(ref_dnam), "CpGs\n\n"
)


# ------------------------------------------------------------------------------
# Select soft-thresholding power
# ------------------------------------------------------------------------------

cat("Evaluating candidate soft-thresholding powers...\n\n")

sft <- pickSoftThreshold(
    ref_dnam,
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


# ------------------------------------------------------------------------------
# Save soft-thresholding results
# ------------------------------------------------------------------------------

fwrite(
    as.data.table(sft$fitIndices),
    file.path(
        output_path,
        "soft_threshold_fit_indices.txt"
    ),
    sep = "\t"
)


# ------------------------------------------------------------------------------
# Plot scale-free topology fit
# ------------------------------------------------------------------------------

png(
    file.path(
        output_path,
        "scale_free_fit_plot.png"
    ),
    width = 800,
    height = 600
)

plot(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    xlab = "Soft-thresholding power",
    ylab = "Scale-free topology model fit, signed R²",
    type = "n"
)

text(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    labels = sft$fitIndices[, 1],
    cex = 0.9
)

abline(
    h = 0.8,
    lty = 2
)

dev.off()


cat("\nSoft-threshold analysis complete.\n")

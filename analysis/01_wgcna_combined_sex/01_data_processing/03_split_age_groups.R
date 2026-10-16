# ==============================================================================
# Split residualised DNAm data into age groups
#
# Chromosome-specific residualised DNA methylation data are divided into eight
# age groups for downstream co-methylation network analysis.
#
# The chromosome number is supplied by the accompanying SLURM script through
# replacement of the @chr@ placeholder.
# ==============================================================================

suppressMessages(library(data.table))


# ------------------------------------------------------------------------------
# User-defined paths
# ------------------------------------------------------------------------------

dnam_path <- "..."
cov_path  <- "..."
out_path  <- "..."


# ------------------------------------------------------------------------------
# Start analysis
# ------------------------------------------------------------------------------

cat("Analysis started:", as.character(Sys.time()), "\n\n")
ptm <- proc.time()


# ------------------------------------------------------------------------------
# Define age groups
#
# Intervals are left-closed and right-open:
#   lower age <= age < upper age 
#   age-groups chosen to have comparable sample per age-group   
# ------------------------------------------------------------------------------

age.grps <- matrix(
    c(
        17.9, 36.0,  # Age group 1
        36.0, 44.0,  # Age group 2
        44.0, 49.5,  # Age group 3
        49.5, 54.2,  # Age group 4
        54.2, 57.8,  # Age group 5
        57.8, 61.5,  # Age group 6
        61.5, 66.0,  # Age group 7
        66.0, 94.0   # Age group 8
    ),
    ncol = 2,
    byrow = TRUE
)


# ------------------------------------------------------------------------------
# Read chromosome-specific residualised DNAm data
# ------------------------------------------------------------------------------

chr <- @chr@

fname <- file.path(
    dnam_path,
    paste0("precorr_gs_qc_chr", chr, ".txt")
)

cat("Processing chromosome", chr, "\n")

dnam <- fread(fname)


# Round residualised DNAm values to reduce memory requirements
dnam[, 2:ncol(dnam) := lapply(.SD, round, digits = 3),
     .SDcols = 2:ncol(dnam)]


# ------------------------------------------------------------------------------
# Read covariates
# ------------------------------------------------------------------------------

all.covar <- fread(
    file.path(cov_path, "GS_18869_all.cov")
)

# Round numeric covariates used downstream
all.covar[, 7:ncol(all.covar) := lapply(.SD, round, digits = 3),
          .SDcols = 7:ncol(all.covar)]


# ------------------------------------------------------------------------------
# Merge residualised DNAm with sample information
# ------------------------------------------------------------------------------

m <- merge(
    all.covar[, 2:ncol(all.covar)],
    dnam,
    by = "IID_DNAm"
)


# ------------------------------------------------------------------------------
# Split samples into eight age groups
# ------------------------------------------------------------------------------

for (i in seq_len(nrow(age.grps))) {

    lower_age <- age.grps[i, 1]
    upper_age <- age.grps[i, 2]

    dnam.sub.age <- data.table(
        m[
            age >= lower_age &
            age < upper_age
        ]
    )

    cat(
        "Age group", i,
        ":", lower_age, "-", upper_age,
        "| N =", nrow(dnam.sub.age), "\n"
    )

    # Display sex distribution as a basic check
    print(table(dnam.sub.age$sex))

    output_file <- file.path(
        out_path,
        paste0(
            tools::file_path_sans_ext(basename(fname)),
            "_age_grp", i,
            ".txt"
        )
    )

    fwrite(
        dnam.sub.age,
        output_file,
        quote = FALSE,
        row.names = FALSE,
        sep = "\t"
    )
}


# ------------------------------------------------------------------------------
# Done!
# ------------------------------------------------------------------------------

cat("\nAnalysis ended:", as.character(Sys.time()), "\n\n")
proc.time() - ptm

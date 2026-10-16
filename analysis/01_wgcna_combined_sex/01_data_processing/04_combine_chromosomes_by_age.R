# ==============================================================================
# Combine chromosome-specific DNAm matrices within each age group
#
# Residualised DNAm data were previously split into eight age groups separately
# for each chromosome. This script combines chromosomes 1-22 to create one
# whole-genome DNAm matrix for each age group.
#
# Sample covariates are retained from chromosome 1 only. For chromosomes 2-22,
# only the sample identifier and CpG methylation columns are retained.
# ==============================================================================

suppressMessages(library(data.table))


# ------------------------------------------------------------------------------
# User-defined paths
# ------------------------------------------------------------------------------

data_path <- "..."
out_path  <- "..."


# ------------------------------------------------------------------------------
# Age group supplied by the accompanying SLURM script
# ------------------------------------------------------------------------------

age_groups <- @AgeGrp@


# Covariate columns present in each chromosome-specific input file
covariate_cols <- c(
    "IID_DNAm",
    "age",
    "sex",
    "Batch",
    "ever_smoke",
    "pack_years",
    "Bcell",
    "CD4T",
    "CD8T",
    "Eos",
    "Mono",
    "Neu",
    "NK"
)


# ------------------------------------------------------------------------------
# Combine chromosomes for each requested age group
# ------------------------------------------------------------------------------

for (age in age_groups) {

    message("Processing age group: ", age)


    # --------------------------------------------------------------------------
    # Identify chromosome-specific files
    # --------------------------------------------------------------------------

    pattern <- sprintf(
        "precorr_gs_qc_chr\\d+_age_grp%d\\.txt",
        age
    )

    files <- list.files(
        path = data_path,
        pattern = pattern
    )

    if (length(files) == 0) {
        stop("No chromosome files found for age group ", age)
    }


    # --------------------------------------------------------------------------
    # Order files by chromosome number
    # --------------------------------------------------------------------------

    chr_number <- function(f) {
        as.integer(
            sub(
                ".*chr([0-9]+)_.*",
                "\\1",
                f
            )
        )
    }

    files <- files[order(chr_number(files))]


    # Confirm that all 22 autosomes are present
    if (length(files) != 22) {
        warning(
            "Expected 22 chromosome files for age group ",
            age,
            " but found ",
            length(files)
        )
    }


    # --------------------------------------------------------------------------
    # Read chromosome 1
    #
    # Covariates are retained from chromosome 1 and are not duplicated when
    # subsequent chromosomes are added.
    # --------------------------------------------------------------------------

    chr1_file <- files[chr_number(files) == 1]

    if (length(chr1_file) != 1) {
        stop(
            "Expected exactly one chromosome 1 file for age group ",
            age
        )
    }

    dt_merged <- fread(
        file.path(data_path, chr1_file)
    )


    # --------------------------------------------------------------------------
    # Add DNAm data from chromosomes 2-22
    # --------------------------------------------------------------------------

    other_files <- files[chr_number(files) != 1]

    for (file in other_files) {

        message("  Adding ", file)

        dt <- fread(
            file.path(data_path, file)
        )

        # Remove repeated covariates while retaining IID_DNAm and CpG columns
        covariates_to_remove <- setdiff(
            covariate_cols,
            "IID_DNAm"
        )

        dt[, (covariates_to_remove) := NULL]

        # Merge by sample identifier
        dt_merged <- merge(
            dt_merged,
            dt,
            by = "IID_DNAm",
            all = TRUE
        )
    }


    # --------------------------------------------------------------------------
    # Save whole-genome age-group matrix
    # --------------------------------------------------------------------------

    out_file <- sprintf(
        "precorr_gs_qc_chr_all_age_grp%d.txt",
        age
    )

    fwrite(
        dt_merged,
        file.path(out_path, out_file),
        sep = "\t"
    )

    message("Created: ", out_file, "\n")
}

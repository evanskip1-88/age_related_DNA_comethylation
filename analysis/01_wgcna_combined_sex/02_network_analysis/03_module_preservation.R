# ==============================================================================
# Assess preservation of reference co-methylation modules across age groups
#
# Modules defined in the youngest age group (Age group 1) are evaluated in
# Age groups 2-8 using WGCNA::modulePreservation().
#
# The grey module (module 0; unassigned CpGs) is excluded.
#
# Preservation statistics are calculated without permutations
# (nPermutations = 0); therefore, the output contains observed preservation
# statistics rather than permutation-based Z-summary statistics.
# ==============================================================================

suppressMessages(library(WGCNA))
suppressMessages(library(data.table))


# ------------------------------------------------------------------------------
# User-defined paths
# ------------------------------------------------------------------------------

dnam_inpath <- "..." # Individual level restricted DNAm data
output_path <- "/results/01_wgcna_combined_sex"


# ------------------------------------------------------------------------------
# Input files
# ------------------------------------------------------------------------------

ref_fname <- file.path(
    dnam_inpath,
    "precorr_gs_qc_chr_all_age1.txt"
)

module_file <- file.path(
    output_path,
    "precorr_gs_qc_chr_all_refage1_wg.modLabs"
)


# ------------------------------------------------------------------------------
# Load reference module assignments
# ------------------------------------------------------------------------------

dt_modules <- fread(module_file)


# Count CpGs per module
dt_module_counts <- as.data.table(
    table(dt_modules$moduleNumID)
)

colnames(dt_module_counts) <- c(
    "moduleNumID",
    "probeCounts"
)


# Exclude grey/unassigned module
dt_select_mods <- dt_modules[
    moduleNumID != 0
]


cat(
    "Testing",
    length(unique(dt_select_mods$moduleNumID)),
    "reference modules across Age groups 2-8\n\n"
)


# ------------------------------------------------------------------------------
# Load reference DNAm matrix
# ------------------------------------------------------------------------------

ref_dnam <- fread(
    ref_fname,
    select = dt_select_mods$probe_id
)

ref_dnam <- as.matrix(ref_dnam)


# ------------------------------------------------------------------------------
# Assess module preservation across age groups
# ------------------------------------------------------------------------------

preservation_results <- list()


for (age_grp in 2:8) {

    cat("Testing preservation in Age group", age_grp, "...\n")

    test_fname <- file.path(
        dnam_inpath,
        paste0(
            "precorr_gs_qc_chr_all_age",
            age_grp,
            ".txt"
        )
    )


    # Load CpGs belonging to reference modules
    test_dnam <- fread(
        test_fname,
        select = dt_select_mods$probe_id
    )

    test_dnam <- as.matrix(test_dnam)


    # --------------------------------------------------------------------------
    # WGCNA module preservation
    # --------------------------------------------------------------------------

    preservation <- modulePreservation(
        multiData = list(
            ref = list(data = ref_dnam),
            test = list(data = test_dnam)
        ),
        multiColor = list(
            ref = dt_select_mods$moduleColors
        ),
        dataIsExpr = TRUE,
        nPermutations = 0,
        verbose = 5
    )

    # --------------------------------------------------------------------------
    # Extract observed preservation statistics
    # --------------------------------------------------------------------------

    observed <- data.table(
        moduleColor = as.factor(
            rownames(
                preservation$observed$ref$inter$ref_vs_test
            )
        ),
        refAge = "Age1",
        testAge = paste0("Age", age_grp),

        cbind(
            preservation$observed$ref$quality[[2]][, 1:2],

            preservation$observed$ref$intra$ref_vs_test[
                , c(2, 1)
            ],

            preservation$observed$ref$inter$ref_vs_test[
                , c(1, 2, 4)
            ]
        )
    )


    # Remove unused preservation statistics
    observed[, c(
        "propVarExplained.qual",
        "meanSignAwareKME.pres"
    ) := NULL]


    # Rename preservation variance-explained statistic
    setnames(
        observed,
        old = names(observed)[5],
        new = "propVarExplained"
    )


    preservation_results[[age_grp]] <- observed

    cat("Done.\n\n")
}

# ------------------------------------------------------------------------------
# Combine age-group results
# ------------------------------------------------------------------------------

preservation_results <- rbindlist(
    preservation_results,
    use.names = TRUE,
    fill = TRUE
)


# Round preservation statistics
numeric_cols <- 5:ncol(preservation_results)

preservation_results[
    ,
    (numeric_cols) := lapply(.SD, round, digits = 4),
    .SDcols = numeric_cols
]

# ------------------------------------------------------------------------------
# Save results
# ------------------------------------------------------------------------------

fwrite(
    preservation_results,
    file.path(
        output_path,
        "precorr_gs_qc_chr_all_refage1_wg.obsStat"
    ),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t",
    col.names = TRUE
)


cat("\nModule preservation analysis complete.\n")

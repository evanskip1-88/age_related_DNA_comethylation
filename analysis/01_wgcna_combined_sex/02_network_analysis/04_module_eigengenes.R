#==============================================================================
#
# Calculate WGCNA module eigengenes across age groups
#
# Modules were defined in the youngest reference age group (Age1).
# The same reference module assignments are used to calculate module
# eigengenes and average module profiles in each test age group.
#
#==============================================================================

suppressPackageStartupMessages({
  library(WGCNA)
  library(data.table)
})

rm(list = ls())

#------------------------------------------------------------------------------
# Paths
#------------------------------------------------------------------------------

dnam_inpath <- "..." # Individual level restricted DNAm data
output_path <- "results/01_wgcna_combined_sex"

# Reference module assignments
module_file <- file.path(
  output_path,
  "precorr_gs_qc_chr_all_refage1_wg.modLabs"
)

mod_labs <- fread(module_file)

#------------------------------------------------------------------------------
# Calculate module eigengenes
#------------------------------------------------------------------------------

for (ageGrp in 1:8) {

  cat("\nProcessing age group", ageGrp, "\n")

  dnam_file <- file.path(
    dnam_inpath,
    paste0(
      "precorr_gs_qc_chr_all_age",
      ageGrp,
      "_roundoff.txt"
    )
  )

  cat("Loading DNAm data...\n")

  # The large DNAm matrix was read in row blocks to reduce memory usage.
  dt_header <- fread(dnam_file, nrows = 5)

  dnam_part1 <- fread(
    dnam_file,
    skip = 0,
    nrows = 300
  )

  dnam_part2 <- fread(
    dnam_file,
    skip = 300,
    nrows = 300
  )
  setnames(dnam_part2, colnames(dt_header))

  dnam_part3 <- fread(
    dnam_file,
    skip = 600
  )
  setnames(dnam_part3, colnames(dt_header))

  # Recombine samples.
  dnam_data <- rbind(
    dnam_part1,
    dnam_part2,
    dnam_part3
  )

  # The first four columns contain sample/covariate information.
  dnam_matrix <- data.matrix(
    dnam_data[, -c(1:4), with = FALSE]
  )

  rownames(dnam_matrix) <- dnam_data$IID_DNAm

  cat(
    "Loaded",
    nrow(dnam_matrix),
    "samples and",
    ncol(dnam_matrix),
    "CpGs.\n"
  )

  #--------------------------------------------------------------------------
  # Module eigengenes
  #--------------------------------------------------------------------------

  me <- moduleEigengenes(
    dnam_matrix,
    colors = mod_labs$moduleColors,
    excludeGrey = TRUE,
    verbose = 5
  )

  #--------------------------------------------------------------------------
  # Save results
  #--------------------------------------------------------------------------

  fwrite(
    me$eigengenes,
    file.path(
      output_path,
      paste0(
        "precorr_gs_qc_chr_all_age",
        ageGrp,
        "_roundoff.eigengenes"
      )
    ),
    quote = FALSE,
    col.names = TRUE,
    row.names = TRUE,
    sep = "\t"
  )

  fwrite(
    me$averageExpr,
    file.path(
      output_path,
      paste0(
        "precorr_gs_qc_chr_all_age",
        ageGrp,
        "_roundoff.averageExpr"
      )
    ),
    quote = FALSE,
    col.names = TRUE,
    row.names = TRUE,
    sep = "\t"
  )

  cat("Completed age group", ageGrp, "\n")
}

cat("\nAnalysis complete.\n")

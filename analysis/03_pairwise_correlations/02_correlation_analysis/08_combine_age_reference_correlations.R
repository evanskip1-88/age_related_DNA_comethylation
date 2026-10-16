
#!/usr/bin/env Rscript

# =============================================================================
# Combine chromosome-level CpG pairwise correlations across reference ages
#
# Input:
#   precorr_gs_qc_chr{1..22}_Age{1..8}.cor
#
# Outputs:
#   1. precorr_gs_qc_chr_all_AgebaseAll_vsOthers.cor
#   2. precorr_gs_qc_chr_all_Ages_vsOthers_long.cor
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# 1. Configuration
# -----------------------------------------------------------------------------

input_dir  <- "/results/03_pairwise_correlations"
output_dir <- "/results/03_pairwise_correlations"


dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

age_groups <- paste0("Age", 1:8)
chromosomes <- 1:22

wide_list <- vector(
  "list",
  length(age_groups) * length(chromosomes)
)

k <- 0L

# -----------------------------------------------------------------------------
# 2. Read and standardise chromosome-level correlation files
# -----------------------------------------------------------------------------

for (age in age_groups) {

  comparison_ages <- setdiff(age_groups, age)

  for (chr in chromosomes) {

    input_file <- file.path(
      input_dir,
      sprintf("precorr_gs_qc_chr%d_%s.cor", chr, age)
    )

    if (!file.exists(input_file)) {
      stop("Missing input file: ", input_file)
    }

    dt <- fread(input_file)

    required_columns <- c(
      "probe1", "probe2",
      age_groups,
      "Chr", "AgeBase"
    )

    if (!setequal(names(dt), required_columns) ||
        length(names(dt)) != length(required_columns)) {
      stop("Unexpected columns in: ", input_file)
    }

    if (anyNA(dt$AgeBase) || any(dt$AgeBase != age)) {
      stop("AgeBase mismatch in: ", input_file)
    }

    if (anyNA(dt$Chr) ||
        any(as.character(dt$Chr) != as.character(chr))) {
      stop("Chromosome mismatch in: ", input_file)
    }

    # Rename the reference correlation
    setnames(dt, age, "AgeBaseCorr")

    # Ensure comparison columns are in ascending age order
    setcolorder(
      dt,
      c(
        "probe1",
        "probe2",
        "AgeBaseCorr",
        comparison_ages,
        "Chr",
        "AgeBase"
      )
    )

    # Standardise comparison column names for row binding
    setnames(
      dt,
      comparison_ages,
      paste0("Comparison", 1:7)
    )

    k <- k + 1L
    wide_list[[k]] <- dt

    message(
      age, " | chr", chr,
      " | CpG pairs: ", nrow(dt)
    )
  }
}

# -----------------------------------------------------------------------------
# 3. Wide-format correlation table
# -----------------------------------------------------------------------------

dt_wide <- rbindlist(
  wide_list,
  use.names = TRUE
)

setnames(
  dt_wide,
  paste0("Comparison", 1:7),
  paste0("Age", 2:8)
)

setcolorder(
  dt_wide,
  c(
    "probe1",
    "probe2",
    "AgeBaseCorr",
    paste0("Age", 2:8),
    "Chr",
    "AgeBase"
  )
)

wide_file <- file.path(
  output_dir,
  "precorr_gs_qc_chr_all_AgebaseAll_vsOthers.cor"
)

fwrite(
  dt_wide,
  wide_file,
  sep = "\t",
  quote = FALSE
)

message("\nSaved wide-format table: ", wide_file)

# -----------------------------------------------------------------------------
# 4. Generate historical long-format correlation table
# -----------------------------------------------------------------------------

# Melt only the seven comparison columns.

dt_long <- melt(
  dt_wide,
  id.vars = c(
    "probe1",
    "probe2",
    "AgeBaseCorr",
    "AgeBase"
  ),
  measure.vars = paste0("Age", 2:8),
  variable.name = "Comparison",
  value.name = "corr"
)

dt_long[
  ,
  Age.grp := {
    comparison_ages <- setdiff(age_groups, AgeBase[1])

    comparison_ages[
      match(
        as.character(Comparison),
        paste0("Age", 2:8)
      )
    ]
  },
  by = AgeBase
]

# Remove the temporary positional comparison label
dt_long[, Comparison := NULL]

# Match the original six-column structure
setcolorder(
  dt_long,
  c(
    "probe1",
    "probe2",
    "AgeBaseCorr",
    "AgeBase",
    "Age.grp",
    "corr"
  )
)

long_file <- file.path(
  output_dir,
  "precorr_gs_qc_chr_all_Ages_vsOthers_long.cor"
)

fwrite(
  dt_long,
  long_file,
  sep = "\t",
  quote = FALSE
)

message("Saved long-format table: ", long_file)



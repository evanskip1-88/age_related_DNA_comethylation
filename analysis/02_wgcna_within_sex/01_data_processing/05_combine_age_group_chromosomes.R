# ============================================================
# Combine chromosome-level DNAm data within age groups
# Retain covariates from chromosome 1 only
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

# Configure directories
data_path <- "..."
out_path  <- "..."

# Age group placeholder substituted by SLURM
age_groups <- @AgeGrp@

# Sex to analyse
sex <- "M"

# Covariate columns retained from chromosome 1
covariate_cols <- c(
  "IID_DNAm", "age", "sex", "Batch", "ever_smoke",
  "pack_years", "Bcell", "CD4T", "CD8T",
  "Eos", "Mono", "Neu", "NK"
)

for (age in age_groups) {

  message("Processing age group: ", age, " | Sex: ", sex)

  # Identify chromosome files for this age group and sex
  pattern <- sprintf(
    "^precorr_gs_qc_chr[0-9]+_%s_age_grp%d\\.txt$",
    sex, age
  )

  files <- list.files(data_path, pattern = pattern)

  if (length(files) != 22L) {
    stop("Expected 22 chromosome files; found ", length(files))
  }

  # Sort by chromosome number
  chr_numbers <- as.integer(
    sub("^precorr_gs_qc_chr([0-9]+)_.*$", "\\1", files)
  )

  files <- files[order(chr_numbers)]
  chr_numbers <- sort(chr_numbers)

  if (!identical(chr_numbers, 1:22)) {
    stop("Missing or duplicated chromosomes for age group ", age)
  }

  # Read chromosome 1, retaining its covariates
  dt_merged <- fread(file.path(data_path, files[1]))

  # Merge remaining chromosomes using sample IDs
  for (file in files[-1]) {

    message("Merging ", file)

    dt <- fread(file.path(data_path, file))

    # Keep sample ID and CpG columns only
    dt <- dt[, !covariate_cols[
      covariate_cols != "IID_DNAm"
    ], with = FALSE]

    dt_merged <- merge(
      dt_merged,
      dt,
      by = "IID_DNAm",
      all = TRUE
    )
  }

  # Write combined age-group dataset
  out_file <- sprintf(
    "precorr_gs_qc_chr_all_age_grp%d_%s.txt",
    age, sex
  )

  fwrite(
    dt_merged,
    file.path(out_path, out_file),
    sep = "\t",
    quote = FALSE
  )

  message("Completed: ", out_file)
}

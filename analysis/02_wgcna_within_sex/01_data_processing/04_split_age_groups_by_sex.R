# ============================================================
# Split pre-corrected DNAm data into age groups within sex
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

# Configure input and output directories
dnam_path <- "..."
cov_path  <- "..."
out_path  <- "..."

cat("Analysis started:", as.character(Sys.time()), "\n")
ptm <- proc.time()

# Age group boundaries: lower inclusive, upper exclusive
age.grps <- matrix(
  c(
    17.9, 36.0,
    36.0, 44.0,
    44.0, 49.5,
    49.5, 54.2,
    54.2, 57.8,
    57.8, 61.5,
    61.5, 66.0,
    66.0, 94.0
  ),
  ncol = 2,
  byrow = TRUE
)

# Chromosome placeholder substituted by SLURM
fname <- file.path(
  dnam_path,
  paste0("precorr_gs_qc_chr", @chr@, "_M.txt")
)

# Read pre-corrected DNAm residuals
dnam <- fread(fname)

# Round DNAm residuals to three decimal places
dnam[, (2:ncol(dnam)) := lapply(.SD, round, digits = 3),
     .SDcols = 2:ncol(dnam)]

# Read covariates
all.covar <- fread(file.path(cov_path, "GS_18869_all.cov"))

# Round numeric covariates to three decimal places
numeric_cols <- names(all.covar)[
  7:ncol(all.covar)
]

numeric_cols <- numeric_cols[
  vapply(all.covar[, ..numeric_cols], is.numeric, logical(1))
]

all.covar[, (numeric_cols) := lapply(.SD, round, digits = 3),
          .SDcols = numeric_cols]

# Merge DNAm residuals with covariates
m <- merge(
  all.covar[, 2:ncol(all.covar), with = FALSE],
  dnam,
  by = "IID_DNAm"
)

# Split into eight age groups
for (i in 1:8) {

  lower <- age.grps[i, 1]
  upper <- age.grps[i, 2]

  dnam.sub.age <- m[age >= lower & age < upper]

  cat("Age group:", i, "\n")
  print(table(dnam.sub.age$sex))

  out_file <- file.path(
    out_path,
    paste0(
      tools::file_path_sans_ext(basename(fname)),
      "_age_grp", i, ".txt"
    )
  )

  fwrite(
    dnam.sub.age,
    out_file,
    quote = FALSE,
    sep = "\t",
    na = "NA"
  )
}

cat("Analysis ended:", as.character(Sys.time()), "\n")
print(proc.time() - ptm)

# Pre-correction of DNAm data within sex
# Linear regression residualisation by chromosome and CpG block

suppressPackageStartupMessages({
  library(data.table)
  library(magrittr)
  library(WGCNA)
  library(tidyverse)
})

# Input directories: configure for your environment
dnam_inpath <- "..." # Individual level data
cov_inpath  <- "..."

cat("Analysis started:", as.character(Sys.time()), "\n")
ptm <- proc.time()

# Load covariates once
all.covar <- fread(file.path(cov_inpath, "GS_18869_all.cov"))

all.covar <- all.covar %>%
  separate(Batch, into = c("BatchGrp", "BatchInt"), remove = FALSE)

all.covar$BatchGrp <- case_when(
  all.covar$BatchGrp == "W1" ~ "1",
  all.covar$BatchGrp == "W2" ~ "2",
  all.covar$BatchGrp == "W3" ~ "3",
  all.covar$BatchGrp == "W4" ~ "4"
)

for (chr in 1:22) {

  fname <- file.path(dnam_inpath, paste0("chr", chr, "_gs_qc_tr.txt"))

  if (!file.exists(fname)) next

  cat("Chromosome:", chr, "\n")

  dnam <- fread(fname)

  # Block range is substituted by the SLURM script
  dnam.blck <- dnam[@block@, ][!is.na(ID), ]

  if (nrow(dnam.blck) == 0) next

  dnam.blck.t <- transposeBigData(
    dnam.blck,
    blocksize = 1000
  ) %>%
    tibble::rownames_to_column("ID")

  colnames(dnam.blck.t) <- dnam.blck.t[1, ]
  dnam.blck.t <- dnam.blck.t[-1, ]

  m <- merge(
    all.covar[, 2:ncol(all.covar)],
    dnam.blck.t,
    by.x = "IID_DNAm",
    by.y = "ID"
  )

  m$ever_smoke <- ifelse(
    is.na(m$ever_smoke),
    median(m$ever_smoke, na.rm = TRUE),
    m$ever_smoke
  )

  idx <- ncol(all.covar[, 2:ncol(all.covar)]) + 1

  # Residualise separately within males and females
  for (sex_id in c("M", "F")) {

    cat("Analysis for:", sex_id, "\n")

    m.sub <- as.data.table(m[m$sex == sex_id, ])

    if (nrow(m.sub) == 0) next

    fwrite(
      m.sub[, 1],
      paste0("gs_qc_tr_lm_chr", chr, "_", sex_id, ".IDs"),
      sep = "\t"
    )

    lm_fun <- function(df) {

      lm.res <- vector("list", ncol(df) - idx + 1)

      for (i in idx:ncol(df)) {

        fit <- lm(
          df[[i]] ~ age + ever_smoke + BatchGrp +
            Bcell + CD4T + CD8T + Eos + Mono + Neu + NK,
          data = df,
          na.action = na.exclude
        )

        lm.res[[i - idx + 1]] <- residuals(fit)
      }

      lm.res
    }

    res.final <- as.data.table(
      do.call(cbind, lm_fun(m.sub))
    )

    setnames(
      res.final,
      names(m.sub)[idx:ncol(m.sub)]
    )

    fwrite(
      res.final,
      paste0("pre_corr_gs_@block@_chr", chr, "_", sex_id, ".res"),
      sep = "\t",
      na = "NA"
    )
  }
}

cat("Analysis ended:", as.character(Sys.time()), "\n")
print(proc.time() - ptm)

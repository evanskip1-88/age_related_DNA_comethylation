# ==============================================================================
# Residualise DNA methylation data prior to downstream analysis
#
# For each CpG, DNA methylation levels are adjusted for age, sex, smoking,
# batch group and estimated blood cell proportions.
#
# The script is designed for chromosome-specific DNAm files and block-wise
# execution on an HPC system. The placeholder @block@ is replaced by the
# corresponding CpG row range by the accompanying SLURM script.
# ==============================================================================

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(WGCNA))


# ------------------------------------------------------------------------------
# User-defined paths
# ------------------------------------------------------------------------------

dnam_inpath <- "..."
cov_inpath  <- "..."


# ------------------------------------------------------------------------------
# Start analysis
# ------------------------------------------------------------------------------

cat("Analysis started:", as.character(Sys.time()), "\n\n")
ptm <- proc.time()


# ------------------------------------------------------------------------------
# Read covariates
# ------------------------------------------------------------------------------

all.covar <- fread(file.path(cov_inpath, "GS_18869_all.cov"))

# Create batch-group variable from the original batch identifier
all.covar <- all.covar %>%
    separate(
        Batch,
        into = c("BatchGrp", "BatchInt"),
        remove = FALSE
    )

all.covar$BatchGrp <- case_when(
    all.covar$BatchGrp == "W1" ~ "1",
    all.covar$BatchGrp == "W2" ~ "2",
    all.covar$BatchGrp == "W3" ~ "3",
    all.covar$BatchGrp == "W4" ~ "4"
)


# ------------------------------------------------------------------------------
# Residualise DNAm chromosome by chromosome
# ------------------------------------------------------------------------------

for (chr in 1:22) {

    fname <- file.path(
        dnam_inpath,
        paste0("chr", chr, "_gs_qc_tr.txt")
    )

    # Skip chromosome if input file is unavailable
    if (!file.exists(fname)) {
        next
    }

    cat("Processing chromosome", chr, "\n")

    dnam <- fread(fname)


    # --------------------------------------------------------------------------
    # Select CpGs belonging to the current block
    # --------------------------------------------------------------------------

    dnam.blck <- dnam[@block@, ][!is.na(ID), ]

    if (nrow(dnam.blck) == 0) {
        next
    }


    # --------------------------------------------------------------------------
    # Transpose DNAm data
    # --------------------------------------------------------------------------

    dnam.blck.t <- transposeBigData(
        dnam.blck,
        blocksize = 1000
    ) %>%
        tibble::rownames_to_column("ID")

    colnames(dnam.blck.t) <- dnam.blck.t[1, ]
    dnam.blck.t <- dnam.blck.t[-1, ]


    # --------------------------------------------------------------------------
    # Merge DNAm data with sample covariates
    # --------------------------------------------------------------------------

    m <- merge(
        all.covar[, 2:ncol(all.covar)],
        dnam.blck.t,
        by.x = "IID_DNAm",
        by.y = "ID"
    )

    # Encode sex as male = 1, female = 0
    m[, sex := ifelse(sex == "M", 1, 0)]

    # Median-impute missing values for ever-smoking status
    m$ever_smoke <- ifelse(
        is.na(m$ever_smoke),
        median(m$ever_smoke, na.rm = TRUE),
        m$ever_smoke
    )

    # First column containing CpG methylation values
    idx <- ncol(all.covar[, 2:ncol(all.covar)]) + 1


    # Save sample IDs corresponding to the residual matrix
    fwrite(
        m[, 1],
        paste0("gs_qc_tr_lm_chr", chr, ".IDs"),
        quote = FALSE,
        row.names = FALSE,
        na = NA,
        col.names = TRUE,
        sep = "\t"
    )


    # --------------------------------------------------------------------------
    # Residualise each CpG
    #
    # Model:
    # DNAm ~ age + sex + ever_smoke + BatchGrp +
    #        Bcell + CD4T + CD8T + Eos + Mono + Neu + NK
    # --------------------------------------------------------------------------

    lm_fun <- function(df) {

        lm.res <- list()

        progress <- ceiling(
            quantile(
                idx:ncol(df),
                seq(0.05, 1, 0.05)
            )
        )

        for (i in idx:ncol(df)) {

            if (i %in% progress) {
                print(
                    paste(
                        "...",
                        names(progress[which(progress == i)])
                    )
                )
            }

            fit <- lm(
                df[[i]] ~
                    age +
                    sex +
                    ever_smoke +
                    BatchGrp +
                    Bcell +
                    CD4T +
                    CD8T +
                    Eos +
                    Mono +
                    Neu +
                    NK,
                na.action = na.exclude,
                data = df
            )

            lm.res[[i]] <- residuals(fit)
        }

        return(lm.res)
    }


    # --------------------------------------------------------------------------
    # Combine and save residualised CpGs
    # --------------------------------------------------------------------------

    res.final <- setDT(
        data.frame(
            Reduce("cbind", lm_fun(m))
        )
    )

    colnames(res.final) <- colnames(m)[idx:ncol(m)]

    cat("Writing residualised chromosome", chr, "results...\n")

    fwrite(
        res.final,
        paste0("pre_corr_gs_@block@_chr", chr, ".res"),
        quote = FALSE,
        row.names = FALSE,
        na = NA,
        col.names = TRUE,
        sep = "\t"
    )
}


# ------------------------------------------------------------------------------
# Done
# ------------------------------------------------------------------------------

cat("\nAnalysis complete.\n")
cat("Analysis ended:", as.character(Sys.time()), "\n\n")

proc.time() - ptm

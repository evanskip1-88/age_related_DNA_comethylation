#==============================================================================
# Subset raw DNAm beta values for co-methylation CpG sets
#
# Extracts CpGs belonging to the age-stable and age-variable pairwise and
# WGCNA sets from chromosome-specific Generation Scotland DNAm M-value files.
#
# M-values are converted to beta values for subsequent OSCA/OREML analysis
# of smoking-related variance explained.
#==============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

rm(list = ls())

#------------------------------------------------------------------------------
# Paths
#------------------------------------------------------------------------------

raw_cpgs_dat <- "..." # Restricted individual level data
output_path <- "results/01_wgcna_combined_sex"
wgcna_cpgs_path <- "results/01_wgcna_combined_sex"
opi_bg_dat_path <- "data/gsDNAm/"

#==============================================================================
# Load CpG annotation
#==============================================================================

opi <- fread(
  paste0(
    opi_bg_dat_path,
    "/GS_sub_qc1.opi"
  )
)

colnames(opi) <- c(
  "chr",
  "probe_id",
  "position",
  "gene",
  "strand"
)

opi <- opi[
  position != "NA"
]

opi[
  ,
  chr := paste0(
    "chr",
    chr
  )
]

#==============================================================================
# Load CpG sets
#==============================================================================

cpg_sets <- readRDS(
  paste0(
    wgcna_cpgs_path,
    "/heteroskid_wgcna_cpg_list.rds"
  )
)

#==============================================================================
# Function to extract raw DNAm values
#==============================================================================

subset_raw_beta <- function(analysis) {

  cat(
    "Processing:",
    analysis,
    "\n"
  )

  cpg_hits <- unique(
    cpg_sets[[analysis]]
  )

  cpg_hits_opi <- opi[
    probe_id %in% cpg_hits
  ]

  chromosomes <- unique(
    cpg_hits_opi$chr
  )

  raw_dnam_list <- list()

  for (chromosome in chromosomes) {

    cat(
      "Subsetting:",
      chromosome,
      "\n"
    )

    chr_cpgs <- unique(
      cpg_hits_opi[
        chr == chromosome,
        probe_id
      ]
    )

    if (length(chr_cpgs) == 0) {
      next
    }

    raw_dnam_list[[chromosome]] <- fread(
      paste0(
        raw_cpgs_dat,
        "/",
        chromosome,
        "_gs_qc.txt"
      ),
      select = c(
        "FID",
        "IID",
        chr_cpgs
      )
    )
  }

  #---------------------------------------------------------------------------
  # Combine chromosome-specific matrices
  #---------------------------------------------------------------------------

  sample_ids <- raw_dnam_list[[1]][
    ,
    .(
      FID,
      IID
    )
  ]

  dnam_values <- lapply(
    raw_dnam_list,
    function(dt) {
      dt[
        ,
        !c(
          "FID",
          "IID"
        ),
        with = FALSE
      ]
    }
  )

  dnam_mvalues <- cbind(
    sample_ids,
    do.call(
      cbind,
      dnam_values
    )
  )

  cpg_ids <- setdiff(
    colnames(dnam_mvalues),
    c(
      "FID",
      "IID"
    )
  )

  #---------------------------------------------------------------------------
  # Convert M-values to beta values
  #---------------------------------------------------------------------------

  dnam_beta <- copy(
    dnam_mvalues
  )

  dnam_beta[
    ,
    (cpg_ids) := lapply(
      .SD,
      function(x) {
        2^x / (2^x + 1)
      }
    ),
    .SDcols = cpg_ids
  ]

  #---------------------------------------------------------------------------
  # Save
  #---------------------------------------------------------------------------

  fwrite(
    dnam_beta,
    paste0(
      output_path,
      "/res_raw_beta_",
      analysis,
      ".txt"
    ),
    sep = "\t",
    quote = FALSE
  )

  cat(
    "Saved:",
    analysis,
    " | samples:",
    nrow(dnam_beta),
    " | CpGs:",
    length(cpg_ids),
    "\n"
  )
}

#==============================================================================
# Run all CpG sets
#==============================================================================

for (analysis in names(cpg_sets)) {

  subset_raw_beta(
    analysis
  )
}

cat("Raw beta-value extraction complete.\n")

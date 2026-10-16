#==============================================================================
# GoDMC meQTL enrichment analysis
#
# Tests enrichment of age-related co-methylation CpG sets among cis- and
# trans-meQTL-associated CpGs from GoDMC.
#
# Because GoDMC was generated using the Illumina HumanMethylation450 array,
# study CpG sets and the analysis background are restricted to CpGs represented
# on both the GS EPIC dataset and the 450K array.
#
# Enrichment is assessed using a one-sided hypergeometric test.
#==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

rm(list = ls())

#------------------------------------------------------------------------------
# Paths
#------------------------------------------------------------------------------

input_path <- "results/01_wgcna_combined_sex"
meqtls_path <- "data/goDMC"
opi_cpgs <- "data/gsDNAm"
output_path <- "analysis/01_wgcna_combined_sex"

#==============================================================================
# Define 450K-compatible analysis background
#==============================================================================

# Illumina HumanMethylation450 manifest
manifest_450k <- fread(
  paste0(
    meqtls_path,
    "/humanmethylation450_15017482_v1-2.csv"
  ),
  skip = 7,
  fill = TRUE
)[, 1:3]

# CpGs retained after QC in the GS EPIC dataset
epic_750k <- fread(
  paste0(
    opi_cpgs,
    "/GS_sub_qc1_h19_cpg.list"
  ),
  header = FALSE
)

# Restrict the background to GS CpGs represented on the 450K array
background_cpgs <- epic_750k[
  V1 %in% manifest_450k$IlmnID
]

background_cpg_vec <- unique(
  background_cpgs$V1
)

whole_cpgs <- length(
  background_cpg_vec
)

#==============================================================================
# Load age-related co-methylation CpG sets
#==============================================================================

cpgs_hit_sets <- readRDS(
  paste0(
    input_path,
    "/heteroskid_wgcna_cpg_list.rds"
  )
)

#------------------------------------------------------------------------------
# Add unclassified WGCNA modules ("Others")
#------------------------------------------------------------------------------

wgcna_mods <- setDT(
  read_xlsx(
    paste0(
      input_path,
      "/01a_selected_mods_summary_stats_wgcna.xlsx"
    )
  )
)

corkME_Other_cpgs <- unique(
  wgcna_mods[
    presStats == "corkME" &
      Category == "Others",
    probe_id
  ]
)

propVar_Other_cpgs <- unique(
  wgcna_mods[
    presStats == "propVarExplained" &
      Category == "Others",
    probe_id
  ]
)

propVar_corkME_Other_cpgs <- intersect(
  corkME_Other_cpgs,
  propVar_Other_cpgs
)

cpgs_hit_sets$WGCNA_corkME_Others <-
  corkME_Other_cpgs

cpgs_hit_sets$WGCNA_propVar_Others <-
  propVar_Other_cpgs

cpgs_hit_sets$WGCNA_propVar_corkME_Others <-
  propVar_corkME_Other_cpgs

# Restrict every study CpG set to the 450K-compatible background
cpgs_hit_sets_bg <- lapply(
  cpgs_hit_sets,
  function(x) {
    intersect(
      unique(x),
      background_cpg_vec
    )
  }
)

#==============================================================================
# Load GoDMC meQTL results
#==============================================================================

godmc_meqtls <- setDT(
  read_xlsx(
    paste0(
      meqtls_path,
      "/Josine_et_alEMS131971-supplement-Supplementary_Data_1.xlsx"
    ),
    skip = 1
  )
)

colnames(godmc_meqtls)[2] <- "DNAm_site"

# Restrict GoDMC results to the eligible CpG background
godmc_meqtls <- godmc_meqtls[
  DNAm_site %in% background_cpg_vec
]

#------------------------------------------------------------------------------
# Select top 10,000 cis- and trans-meQTL associations
#------------------------------------------------------------------------------

dt_meqtls_filtered_cis <- godmc_meqtls[
  Cis == "TRUE"
][
  order(Replication_P)
][
  1:10000
]

godmc_cis_meqtls_cpgs <- unique(
  dt_meqtls_filtered_cis$DNAm_site
)

dt_meqtls_filtered_trans <- godmc_meqtls[
  Cis == "FALSE"
][
  order(Replication_P)
][
  1:10000
]

godmc_trans_meqtls_cpgs <- unique(
  dt_meqtls_filtered_trans$DNAm_site
)

#==============================================================================
# Enrichment function
#==============================================================================

perform_enrichment <- function(
    study_hit_list,
    pub_cpg_list,
    bg_set,
    eqtl_type,
    study_data) {

  results <- list()

  for (study_set in names(study_hit_list)) {

    study_hits <- study_hit_list[[study_set]]

    overlaps <- intersect(
      study_hits,
      pub_cpg_list
    )

    # Contingency-table counts
    N <- bg_set
    K <- length(pub_cpg_list)
    n <- length(study_hits)
    k <- length(overlaps)

    # Odds ratio and 95% confidence interval
    if (
      k == 0 ||
      K == k ||
      n == k ||
      (N - K - (n - k)) == 0
    ) {

      OR <- NA
      CI_lower <- NA
      CI_upper <- NA

    } else {

      OR <- (k / (K - k)) /
        ((n - k) / (N - K - (n - k)))

      SE_log_OR <- sqrt(
        1 / k +
        1 / (n - k) +
        1 / (K - k) +
        1 / (N - K - (n - k))
      )

      log_OR <- log(OR)

      CI_lower <- exp(
        log_OR - 1.96 * SE_log_OR
      )

      CI_upper <- exp(
        log_OR + 1.96 * SE_log_OR
      )
    }

    results[[study_set]] <- data.frame(
      analysis = study_set,
      study_data = study_data,
      eqtl_type = eqtl_type,
      publishedHits = K,
      study_hits = n,
      overlap = k,
      pvalue = phyper(
        k - 1,
        K,
        N - K,
        n,
        lower.tail = FALSE
      ),
      OR = OR,
      CI_lower = CI_lower,
      CI_upper = CI_upper
    )
  }

  setDT(
    do.call(
      rbind,
      results
    )
  )
}

#==============================================================================
# Run GoDMC enrichment
#==============================================================================

# Cis-meQTL enrichment
res_cis_godmc <- perform_enrichment(
  study_hit_list = cpgs_hit_sets_bg,
  pub_cpg_list = godmc_cis_meqtls_cpgs,
  bg_set = whole_cpgs,
  eqtl_type = "Cis",
  study_data = "GoDMC"
)

# Trans-meQTL enrichment
res_trans_godmc <- perform_enrichment(
  study_hit_list = cpgs_hit_sets_bg,
  pub_cpg_list = godmc_trans_meqtls_cpgs,
  bg_set = whole_cpgs,
  eqtl_type = "Trans",
  study_data = "GoDMC"
)

# Combine cis- and trans-meQTL results
res_cis_trans_godmc <- rbind(
  res_cis_godmc,
  res_trans_godmc
)

#==============================================================================
# Save results
#==============================================================================

fwrite(
  res_cis_trans_godmc,
  paste0(
    output_path,
    "/01a_res_by_study_meqtl_enrichments.tsv"
  ),
  sep = "\t",
  quote = FALSE
)

cat("GoDMC meQTL enrichment analysis complete.\n")

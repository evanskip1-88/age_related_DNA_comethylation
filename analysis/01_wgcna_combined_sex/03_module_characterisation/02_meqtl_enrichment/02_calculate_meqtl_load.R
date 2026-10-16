#==============================================================================
# Calculate cis- and trans-meQTL load across WGCNA modules
#
# For each age-stable or age-variable WGCNA module, calculates the proportion
# of module CpGs represented among the top 10,000 cis- and trans-meQTL
# associations from GoDMC.
#
# NB:
# This analysis reflects the module-level meQTL load calculation used in the
# submitted manuscript. Unlike the primary GoDMC enrichment analysis, module
# CpGs were not restricted to the 450K-compatible CpG background here.
# Future analyses should harmonise the CpG universe before calculating load.
#==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

rm(list = ls())

#------------------------------------------------------------------------------
# Paths
#------------------------------------------------------------------------------

output_path <- "results/01_wgcna_combined_sex/"
wgcna_cpgs <- "results/01_wgcna_combined_sex"
opi_cpgs <- "data/gsDNAm"
meqtls_path <- "data/goDMC"

#==============================================================================
# Define 450K-compatible CpG universe
#==============================================================================

manifest_450k <- fread(
  paste0(
    meqtls_path,
    "/humanmethylation450_15017482_v1-2.csv"
  ),
  skip = 7,
  fill = TRUE
)[, 1:3]

epic_cpgs <- fread(
  paste0(
    opi_cpgs,
    "/GS_sub_qc1_h19_cpg.list"
  ),
  header = FALSE
)

background_cpgs <- epic_cpgs[
  V1 %in% manifest_450k$IlmnID
]

background_cpg_vec <- unique(
  background_cpgs$V1
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

# Restrict to the eligible 450K-compatible CpG universe
godmc_meqtls <- godmc_meqtls[
  DNAm_site %in% background_cpg_vec
]

#------------------------------------------------------------------------------
# Select top 10,000 cis- and trans-meQTL associations
#------------------------------------------------------------------------------

godmc_cis <- godmc_meqtls[
  Cis == "TRUE"
][
  order(Replication_P)
][
  1:10000
]

godmc_trans <- godmc_meqtls[
  Cis == "FALSE"
][
  order(Replication_P)
][
  1:10000
]

cis_set <- unique(
  godmc_cis$DNAm_site
)

trans_set <- unique(
  godmc_trans$DNAm_site
)

#==============================================================================
# Load selected WGCNA modules
#==============================================================================

wgcna_mods <- setDT(
  read_xlsx(
    paste0(
      wgcna_cpgs,
      "/01a_selected_mods_summary_stats_wgcna.xlsx"
    )
  )
)

# Retain age-stable and age-variable modules
wgcna_mods <- wgcna_mods[
  Category %in% c(
    "Stable",
    "Changing"
  )
]

# Restrict module CpGs to the same eligible CpG universe
wgcna_mods <- wgcna_mods[
  probe_id %in% background_cpg_vec
]

# Create unique module labels
wgcna_mods[
  ,
  ModuleLabel := paste0(
    Category,
    "_",
    moduleNumID
  )
]

# Split CpGs by module
module_cpgs <- split(
  wgcna_mods$probe_id,
  wgcna_mods$ModuleLabel
)

module_cpgs <- lapply(
  module_cpgs,
  unique
)

#==============================================================================
# Calculate meQTL load per module
#==============================================================================

meqtl_load <- rbindlist(
  lapply(
    names(module_cpgs),
    function(module_name) {

      cpgs <- module_cpgs[[module_name]]

      module_size <- length(cpgs)

      module_info <- strsplit(
        module_name,
        "_",
        fixed = TRUE
      )[[1]]

      category <- module_info[1]
      module_id <- module_info[2]

      n_cis <- sum(
        cpgs %in% cis_set
      )

      n_trans <- sum(
        cpgs %in% trans_set
      )

      data.table(
        ModuleID = module_id,
        Category = category,
        ModuleSize_450K = module_size,
        n_cis = n_cis,
        n_trans = n_trans,
        L_cis = n_cis / module_size,
        L_trans = n_trans / module_size,
        TransDominance =
          n_trans / (n_trans + n_cis + 1e-6)
      )
    }
  )
)

#==============================================================================
# Add module-level preservation statistics
#==============================================================================

wgcna_mods$moduleNumID <- as.character(
  wgcna_mods$moduleNumID
)

module_summary <- unique(
  wgcna_mods[
    ,
    .(
      moduleColor,
      moduleNumID,
      presStats,
      moduleSize,
      change_presStats,
      gam_p_value,
      Category
    )
  ]
)

meqtl_load <- merge(
  meqtl_load,
  module_summary,
  by.x = c(
    "ModuleID",
    "Category"
  ),
  by.y = c(
    "moduleNumID",
    "Category"
  )
)

#==============================================================================
# Save results
#==============================================================================

fwrite(
  meqtl_load,
  paste0(
    output_path,
    "/01b_res_all_per_module_cis_trans_geneticLoad.tsv"
  ),
  sep = "\t",
  quote = FALSE
)

cat("Module-level meQTL load analysis complete.\n")

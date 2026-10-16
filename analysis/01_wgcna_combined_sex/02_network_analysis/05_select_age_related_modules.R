
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Select age-stable and age-variable WGCNA co-methylation modules
#
# Purpose:
#   1. Classify Age1-reference WGCNA modules as aSCMs, aVCMs or Others.
#   2. Annotate CpGs and export module-preservation statistics.
#   3. Combine selected WGCNA and heteroskedasticity CpG/gene sets.
#   4. Save the combined lists for downstream enrichment comparisons.
#
# Notes:
#   - Module-size exclusions: >=1,000 CpGs and >100 CpGs.
#   - Input/output paths are placeholders for public release.
#
# Reference:
#   Langfelder et al. (2011), PLoS Comput Biol 7:e1001057.
#   Module-size thresholds are specific to this analysis.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

suppressPackageStartupMessages({
  library(mgcv)
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(ggplot2)
  library(writexl)
})

rm(list = ls())

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

input_dir      <- "results/01_wgcna_combined_sex"
heteroskid_dir <- "results/03_pairwise_correlations"
output_dir     <- "figures/01_wgcna_combined_sex"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load WGCNA preservation statistics and module assignments
# ---------------------------------------------------------------------------

pres_stats <- fread(
  file.path(input_dir, "precorr_gs_qc_chr_all_refage1_wg.obsStat")
)

module_labs <- fread(
  file.path(input_dir, "precorr_gs_qc_chr_all_refage1_wg.modLabs")
)

pres_stats[, testAgeNum := as.numeric(sub("Age", "", testAge))]

setnames(
  pres_stats,
  old = names(pres_stats)[6:8],
  new = c("corkIM", "corkME", "corCor")
)

setnames(
  module_labs,
  old = names(module_labs)[3],
  new = "moduleColor"
)

module_counts <- module_labs[
  ,
  .(moduleSize = .N),
  by = .(moduleNumID, moduleColor)
][order(-moduleSize)]

if (!"moduleSize" %in% names(pres_stats)) {
  pres_stats <- merge(
    pres_stats,
    module_counts[, .(moduleColor, moduleSize)],
    by = "moduleColor",
    all.x = TRUE
  )
}

# ---------------------------------------------------------------------------
# Function: classify modules using preservation trajectories
# ---------------------------------------------------------------------------

select_top_modules <- function(
    pres_obs_data,
    pres_stats,
    connectionCutoff,
    stableSDcutoff,
    gamPvalCutoff,
    changeQuantCutoff,
    changePercentCutoff
) {

  gam_results <- list()

  for (mod in unique(pres_obs_data$moduleColor)) {

    module_data <- pres_obs_data %>%
      dplyr::filter(moduleColor == mod)

    gam_model <- mgcv::gam(
      as.formula(paste0(pres_stats, " ~ s(testAgeNum, k = 4)")),
      data = module_data
    )

    gam_summary <- summary(gam_model)

    gam_results[[mod]] <- data.frame(
      moduleColor = mod,
      edf = gam_summary$s.table[1, "edf"],
      F_value = gam_summary$s.table[1, "F"],
      gam_p_value = gam_summary$s.table[1, "p-value"],
      R2_adj = gam_summary$r.sq
    )
  }

  gam_stats <- dplyr::bind_rows(gam_results) %>%
    dplyr::mutate(
      sigGroup = dplyr::case_when(
        gam_p_value < gamPvalCutoff ~ "aVCMs",
        TRUE ~ "aSCMs"
      )
    )

  pres_obs_data <- merge(
    pres_obs_data,
    gam_stats[, c("moduleColor", "sigGroup", "gam_p_value")],
    by = "moduleColor",
    all.x = TRUE
  )

  setDT(pres_obs_data)

  # Preservation variability across test age groups.
  sd_values <- pres_obs_data[
    ,
    .(sd_presStats = sd(get(pres_stats), na.rm = TRUE)),
    by = moduleColor
  ]

  # Absolute preservation change from Age2 to Age8.
  extreme_changes <- pres_obs_data[
    testAge %in% c("Age2", "Age8"),
    .(change_presStats = abs(diff(get(pres_stats)))),
    by = moduleColor
  ]

  pres_obs_data <- merge(
    pres_obs_data,
    sd_values,
    by = "moduleColor",
    all.x = TRUE
  )

  pres_obs_data <- merge(
    pres_obs_data,
    extreme_changes,
    by = "moduleColor",
    all.x = TRUE
  )

  sd_cutoff <- quantile(
    pres_obs_data$sd_presStats,
    stableSDcutoff,
    na.rm = TRUE
  )

  change_cutoff <- max(
    quantile(
      pres_obs_data$change_presStats,
      changeQuantCutoff,
      na.rm = TRUE
    ),
    changePercentCutoff
  )

  connected_mods <- pres_obs_data[
    testAge == "Age2" & get(pres_stats) > connectionCutoff,
    unique(moduleColor)
  ]

  pres_obs_data <- pres_obs_data %>%
    dplyr::mutate(
      sigGroup = dplyr::case_when(
        sigGroup == "aSCMs" &
          sd_presStats >= sd_cutoff ~ "Others",

        sigGroup == "aSCMs" &
          !(moduleColor %in% connected_mods) ~ "Others",

        sigGroup == "aVCMs" &
          change_presStats < change_cutoff ~ "Others",

        sigGroup == "aVCMs" &
          !(moduleColor %in% connected_mods) ~ "Others",

        TRUE ~ sigGroup
      )
    )

  pres_obs_data_long <- pres_obs_data %>%
    tidyr::pivot_longer(
      cols = c(propVarExplained, corkIM, corkME, corCor),
      names_to = "presStats",
      values_to = "pres_obs_value"
    ) %>%
    as.data.table()

  pres_obs_data_long <- pres_obs_data_long[
    presStats == pres_stats
  ]

  return(pres_obs_data_long)
}

# ---------------------------------------------------------------------------
# Classify modules using proportion of variance explained
# ---------------------------------------------------------------------------

res_propVarExplained <- select_top_modules(
  pres_obs_data = copy(pres_stats),
  pres_stats = "propVarExplained",
  connectionCutoff = 0.40,
  stableSDcutoff = 0.40,
  gamPvalCutoff = 0.01,
  changeQuantCutoff = 0.50,
  changePercentCutoff = 0.10
)

# ---------------------------------------------------------------------------
# Classify modules using correlation of module membership
# ---------------------------------------------------------------------------

res_corkME <- select_top_modules(
  pres_obs_data = copy(pres_stats),
  pres_stats = "corkME",
  connectionCutoff = 0.40,
  stableSDcutoff = 0.40,
  gamPvalCutoff = 0.01,
  changeQuantCutoff = 0.50,
  changePercentCutoff = 0.05
)

# ---------------------------------------------------------------------------
# Combine preservation statistics
# ---------------------------------------------------------------------------

sel_mods_c <- rbind(
  res_propVarExplained,
  res_corkME
)

# ---------------------------------------------------------------------------
# Exclude modules containing >=1,000 CpGs
# ---------------------------------------------------------------------------

# WGCNA modulePreservation() subsampled modules exceeding default
# maxModuleSizev(1000 Cpgs). Analysis excluded modules containing
# >=1,000 CpGs to avoid relying on preservation estimates potentially
# affected by random subsampling.

sel_mods_c <- sel_mods_c[moduleSize < 1000]

# ---------------------------------------------------------------------------
# Identify conflicting preservation classifications
# ---------------------------------------------------------------------------

# Modules classified as aSCMs using propVarExplained but aVCMs
# using corkME are reclassified as Others.

misclassified_mods <- intersect(
  res_propVarExplained[
    sigGroup == "aSCMs",
    unique(moduleColor)
  ],
  res_corkME[
    sigGroup == "aVCMs",
    unique(moduleColor)
  ]
)

sel_mods_c[
  moduleColor %in% misclassified_mods,
  sigGroup := "Others"
]

# Analysis excluded unusually large candidate modules
# (>100 CpGs) from the final aSCM/aVCM classification based on visual
# of module size distribution
#
# These modules remain in the complete table as Others.

dt_module_size <- unique(
  sel_mods_c[
    sigGroup != "Others",
    .(moduleColor, moduleSize, sigGroup)
  ]
)

module_size_outliers <- unique(
  dt_module_size[
    moduleSize > 100,
    moduleColor
  ]
)

# ---------------------------------------------------------------------------
# Assign final module-level categories
# ---------------------------------------------------------------------------

sel_mods_c_final <- copy(sel_mods_c)

sel_mods_c_final[
  ,
  Category := dplyr::case_when(
    moduleColor[1] %in% module_size_outliers ~ "Others",
    any(sigGroup == "aVCMs") ~ "aVCMs",
    any(sigGroup == "aSCMs") ~ "aSCMs",
    TRUE ~ "Others"
  ),
  by = moduleColor
]

sel_mods_unique <- unique(
  sel_mods_c_final[
    ,
    .(moduleColor, Category)
  ]
)

# ---------------------------------------------------------------------------
# Plot preservation trajectories
# ---------------------------------------------------------------------------

p_preservation <- ggplot(
  sel_mods_c_final,
  aes(
    x = testAgeNum,
    y = pres_obs_value,
    color = Category,
    group = moduleColor
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~ presStats, scales = "free_y") +
  scale_color_manual(
    values = c(
      aSCMs = "blue",
      aVCMs = "red",
      Others = "grey70"
    )
  ) +
  labs(
    x = "Test age group",
    y = "Preservation value",
    color = "Module category"
  ) +
  theme_bw()

ggsave(
  file.path(output_dir, "wgcna_preservation_trajectories.pdf"),
  p_preservation,
  width = 10,
  height = 6
)

# ---------------------------------------------------------------------------
# Merge CpG assignments with preservation statistics
# ---------------------------------------------------------------------------

sel_mods_c_final_m <- merge(
  module_labs,
  sel_mods_c_final,
  by = "moduleColor",
  allow.cartesian = TRUE
)

# ---------------------------------------------------------------------------
# Load genomic and functional annotations
# ---------------------------------------------------------------------------

opi <- fread(
  file.path(input_dir, "GS_sub_qc1.opi")
)

setnames(
  opi,
  c("chr", "probe_id", "position", "gene", "strand")
)

opi_rgreat <- fread(
  file.path(
    input_dir,
    "GS_sub_qc1_gene_cpgs_assoc_rGREAT_h19_all_CpGs.opi"
  )
)

setnames(
  opi_rgreat,
  old = names(opi_rgreat)[4],
  new = "gene_rgreat"
)

opi_rgreat <- opi_rgreat[position != "NA"]

dt_cpgs <- fread(
  file.path(input_dir, "GS_sub_qc1_h19_cpgs.annotr")
)

dt_cpgs[
  ,
  annot.type := sub("hg19_cpg_", "", annot.type)
]

dt_genic <- fread(
  file.path(input_dir, "GS_sub_qc1_h19_genic.annotr")
)

dt_genic[
  ,
  annot.type := sub("hg19_genes_", "", annot.type)
]

# ---------------------------------------------------------------------------
# Annotate WGCNA CpGs
# ---------------------------------------------------------------------------

sel_mods_c_final_m2 <- merge(
  sel_mods_c_final_m,
  opi,
  by = "probe_id"
)

sel_mods_c_final_m3 <- merge(
  sel_mods_c_final_m2,
  opi_rgreat[
    ,
    .(probe_id, gene_rgreat, dist_to_TSS)
  ],
  by = "probe_id"
)

sel_mods_c_final_m4 <- merge(
  sel_mods_c_final_m3,
  dt_cpgs[, .(probe_id, annot.type)],
  by = "probe_id",
  all.x = TRUE
)

setnames(
  sel_mods_c_final_m4,
  "annot.type",
  "CpG_annotation"
)

sel_mods_c_final_m5 <- merge(
  sel_mods_c_final_m4,
  dt_genic[, .(probe_id, annot.type)],
  by = "probe_id",
  all.x = TRUE,
  allow.cartesian = TRUE
)

setnames(
  sel_mods_c_final_m5,
  "annot.type",
  "Genic_annotation"
)

# ---------------------------------------------------------------------------
# Collapse overlapping genic annotations
# ---------------------------------------------------------------------------

sel_mods_c_final_m5 <- sel_mods_c_final_m5[
  ,
  .(
    Genic_annotation = paste(
      unique(na.omit(Genic_annotation)),
      collapse = "|"
    )
  ),
  by = setdiff(
    names(sel_mods_c_final_m5),
    "Genic_annotation"
  )
]

# ---------------------------------------------------------------------------
# Validate annotated CpG-level table
# ---------------------------------------------------------------------------

duplicate_records <- sel_mods_c_final_m5[
  ,
  .N,
  by = .(
    probe_id,
    moduleColor,
    testAge,
    presStats
  )
][N > 1]

stopifnot(nrow(duplicate_records) == 0L)

print(
  sel_mods_c_final_m5[
    ,
    .(
      n_records = .N,
      n_cpgs = uniqueN(probe_id),
      n_modules = uniqueN(moduleColor)
    )
  ]
)

# ---------------------------------------------------------------------------
# Export supplementary WGCNA results
# ---------------------------------------------------------------------------

writexl::write_xlsx(
  sel_mods_c_final_m5,
  file.path(
    output_dir,
    "01a_selected_mods_summary_stats_wgcna.xlsx"
  )
)

fwrite(
  sel_mods_unique,
  file.path(output_dir, "wgcna_module_classifications.tsv"),
  sep = "\t"
)

# ---------------------------------------------------------------------------
# Load heteroskedasticity-selected CpGs
# ---------------------------------------------------------------------------

# The top-1,000 stable/changing CpG sets
# from the heteroskedasticity analysis.

heteroskid_Stable <- fread(
  file.path(
    heteroskid_dir,
    "chr_all_AgeBaseAll_stablePosNeg_heteroskid_top1000_cpgs_list.tsv"
  )
)

heteroskid_Changing <- fread(
  file.path(
    heteroskid_dir,
    "chr_all_AgeBaseAll_changingPosNeg_heteroskid_top1000_cpgs_list.tsv"
  )
)

stopifnot(
  "probe_id" %in% names(heteroskid_Stable),
  "probe_id" %in% names(heteroskid_Changing)
)

# ---------------------------------------------------------------------------
# Extract unique semicolon-separated gene symbols
# ---------------------------------------------------------------------------

extract_genes <- function(x) {

  x <- x[!is.na(x) & nzchar(x)]

  if (length(x) == 0L) {
    return(character(0))
  }

  genes <- trimws(
    unlist(strsplit(as.character(x), ";"), use.names = FALSE)
  )

  unique(genes[!is.na(genes) & nzchar(genes)])
}

# ---------------------------------------------------------------------------
# Construct WGCNA gene and CpG sets
# ---------------------------------------------------------------------------

wgcna_stable <- sel_mods_c_final_m5[
  Category == "aSCMs"
]

wgcna_changing <- sel_mods_c_final_m5[
  Category == "aVCMs"
]

wgcna_gene_list <- list(
  WGCNA_Stable = extract_genes(wgcna_stable$gene_rgreat),
  WGCNA_Changing = extract_genes(wgcna_changing$gene_rgreat)
)

wgcna_cpgs_list <- list(
  WGCNA_Stable = unique(wgcna_stable$probe_id),
  WGCNA_Changing = unique(wgcna_changing$probe_id)
)

# ---------------------------------------------------------------------------
# Construct heteroskedasticity gene and CpG sets
# ---------------------------------------------------------------------------

heteroskid_gene_list <- list(
  Heteroskid_Stable = extract_genes(
    opi[
      probe_id %in% heteroskid_Stable$probe_id,
      gene
    ]
  ),
  Heteroskid_Changing = extract_genes(
    opi[
      probe_id %in% heteroskid_Changing$probe_id,
      gene
    ]
  )
)

heteroskid_cpgs_list <- list(
  Heteroskid_Stable = unique(heteroskid_Stable$probe_id),
  Heteroskid_Changing = unique(heteroskid_Changing$probe_id)
)

# ---------------------------------------------------------------------------
# Combine WGCNA and heteroskedasticity sets
# ---------------------------------------------------------------------------

heteroskid_wgcna_gene_list <- c(
  wgcna_gene_list,
  heteroskid_gene_list
)

heteroskid_wgcna_cpgs_list <- c(
  wgcna_cpgs_list,
  heteroskid_cpgs_list
)

expected_groups <- c(
  "WGCNA_Stable",
  "WGCNA_Changing",
  "Heteroskid_Stable",
  "Heteroskid_Changing"
)

stopifnot(
  identical(names(heteroskid_wgcna_gene_list), expected_groups),
  identical(names(heteroskid_wgcna_cpgs_list), expected_groups)
)

# ---------------------------------------------------------------------------
# Save combined gene and CpG lists as RDS files
# ---------------------------------------------------------------------------

# These two files are required for downstream enrichment comparisons,
# gene-set overlaps and WGCNA-vs-heteroskedasticity analyses.

gene_rds <- file.path(
  output_dir,
  "heteroskid_wgcna_gene_list.rds"
)

cpg_rds <- file.path(
  output_dir,
  "heteroskid_wgcna_cpg_list.rds"
)

saveRDS(
  heteroskid_wgcna_gene_list,
  file = gene_rds
)

saveRDS(
  heteroskid_wgcna_cpgs_list,
  file = cpg_rds
)

# ---------------------------------------------------------------------------
# Verify saved RDS files
# ---------------------------------------------------------------------------

stopifnot(
  file.exists(gene_rds),
  file.exists(cpg_rds)
)

gene_list_check <- readRDS(gene_rds)
cpg_list_check <- readRDS(cpg_rds)

stopifnot(
  identical(gene_list_check, heteroskid_wgcna_gene_list),
  identical(cpg_list_check, heteroskid_wgcna_cpgs_list)
)

# ---------------------------------------------------------------------------
# Summarise gene and CpG sets
# ---------------------------------------------------------------------------

set_summary <- data.table(
  Analysis = expected_groups,
  CpGs = lengths(heteroskid_wgcna_cpgs_list),
  Genes = lengths(heteroskid_wgcna_gene_list)
)

print(set_summary)

fwrite(
  set_summary,
  file.path(output_dir, "wgcna_heteroskid_set_summary.tsv"),
  sep = "\t"
)

# <<<<<<<< End of script: e.cheruiyot@uq.edu.au >>>>>>>



#!/usr/bin/env Rscript

# =============================================================================
# Functional enrichment of age-related co-methylation modules using rGREAT
#
# Analyses:
#   - GO Molecular Function (MF)
#   - GO Biological Process (BP)
#   - Enriched genes
#   - Mouse and human phenotype annotations
#
# Background: QC-filtered CpGs with genomic positions (hg19)
# Significance: BH-adjusted P < 0.01
# Output: Top 30 enriched terms per category
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(rGREAT)
  library(readxl)
})

# -----------------------------------------------------------------------------
# 1. Paths and settings
# -----------------------------------------------------------------------------

opi_dir  <- "/data/gsDNAm"
cpgs_dir <- "/results/01_wgcna_combined_sex"
output_dir <- "/results/01_wgcna_combined_sex/functional_enrichment"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

bg_file <- file.path(
  opi_dir,
  "GS_sub_qc1_gene_cpgs_assoc_rGREAT_h19_all_CpGs.opi"
)

module_file <- file.path(
  cpgs_dir,
  "01a_selected_mods_summary_stats.xlsx"
)

# Select "Stable" or "Changing"
modGrp <- "Stable"

statsName <- "cor.kME_n_propVar"

# -----------------------------------------------------------------------------
# 2. Load background CpGs
# -----------------------------------------------------------------------------

bg.beta <- fread(bg_file)

# Preserve original GREAT background construction
bg.beta <- na.omit(
  bg.beta[, .(chr, position, position)]
)

# -----------------------------------------------------------------------------
# 3. Load selected modules
# -----------------------------------------------------------------------------

dt <- as.data.table(read_xlsx(module_file))

dt <- dt[Category == modGrp]

if (nrow(dt) == 0L) {
  stop("No modules found for Category: ", modGrp)
}

# Combined and individual module analyses
module_names <- c("All", unique(dt$moduleColor))

# -----------------------------------------------------------------------------
# 4. Helper: select significant top 30 terms
# -----------------------------------------------------------------------------

select_top_terms <- function(tbl, group_name) {

  res <- as.data.table(tbl)

  res <- res[
    order(Hyper_Raw_PValue)
  ][
    Hyper_Adjp_BH < 0.01
  ][
    1:min(.N, 30)
  ]

  res[, Group := group_name]

  res
}

# -----------------------------------------------------------------------------
# 5. Run rGREAT enrichment
# -----------------------------------------------------------------------------

for (mod in module_names) {

  if (mod == "All") {
    dt_mod <- dt
    modCol <- "All"
  } else {
    dt_mod <- dt[moduleColor == mod]
    modCol <- mod
  }

  message(
    "Running rGREAT: ",
    modGrp, " / ", modCol,
    " (", nrow(dt_mod), " CpG records)"
  )

  bed_mod <- unique(
    dt_mod[
      order(chr),
      .(chr, position, position)
    ]
  )

  # Original GREAT submission settings
  job <- submitGreatJob(
    bed_mod,
    bg = bg.beta,
    request_interval = 60,
    help = FALSE
  )

  # Retrieve enrichment tables
  tb.go <- getEnrichmentTables(
    job,
    category = "GO"
  )

  tb.pheno <- getEnrichmentTables(
    job,
    category = "Phenotype"
  )

  tb.genes <- getEnrichmentTables(
    job,
    category = "Genes"
  )

  # ---------------------------------------------------------------------------
  # 6. Enriched genes
  # ---------------------------------------------------------------------------

  Genes <- select_top_terms(
    tb.genes[[1]],
    "EnsGenes"
  )

  # ---------------------------------------------------------------------------
  # 7. GO enrichment
  # ---------------------------------------------------------------------------

  MF <- select_top_terms(
    tb.go[[1]],
    "MF"
  )

  BP <- select_top_terms(
    tb.go[[2]],
    "BP"
  )

  GO_res <- rbindlist(
    list(MF, BP),
    use.names = TRUE,
    fill = TRUE
  )

  GO_res <- GO_res[
    ,
    .(
      ID,
      name,
      Hyper_Fold_Enrichment,
      Hyper_Raw_PValue,
      Hyper_Adjp_BH,
      Group
    )
  ]

  # ---------------------------------------------------------------------------
  # 8. Phenotype enrichment
  # ---------------------------------------------------------------------------

  mp <- select_top_terms(
    tb.pheno[[1]],
    "Mouse pheno"
  )

  mp_sKO <- select_top_terms(
    tb.pheno[[2]],
    "Mouse pheno sKO"
  )

  hp <- select_top_terms(
    tb.pheno[[3]],
    "Human pheno"
  )

  pheno_res <- rbindlist(
    list(mp, mp_sKO, hp),
    use.names = TRUE,
    fill = TRUE
  )

  pheno_res <- pheno_res[
    ,
    .(
      ID,
      name,
      Hyper_Fold_Enrichment,
      Hyper_Raw_PValue,
      Hyper_Adjp_BH,
      Group
    )
  ]

  # ---------------------------------------------------------------------------
  # 9. Add module metadata
  # ---------------------------------------------------------------------------

  for (res_name in c("Genes", "GO_res", "pheno_res")) {

    res <- get(res_name)

    # Original metadata fields
    res[, sigGroup := paste(unique(dt_mod$Group), collapse = ";")]
    res[, moduleColor := modCol]
    res[, statsPres := statsName]

    assign(res_name, res)
  }

  # ---------------------------------------------------------------------------
  # 10. Save results
  # ---------------------------------------------------------------------------

  output_prefix <- file.path(
    output_dir,
    paste0(
      "01c_", modGrp, "_", modCol,
      "_", statsName, "_rGREAT"
    )
  )

  fwrite(
    Genes,
    paste0(output_prefix, "_ensgenes.out"),
    sep = "\t",
    quote = FALSE
  )

  fwrite(
    pheno_res,
    paste0(output_prefix, "_pheno.out"),
    sep = "\t",
    quote = FALSE
  )

  fwrite(
    GO_res,
    paste0(output_prefix, "_GO.out"),
    sep = "\t",
    quote = FALSE
  )

  message("Completed: ", modCol)
}

message("Functional enrichment analysis completed.")


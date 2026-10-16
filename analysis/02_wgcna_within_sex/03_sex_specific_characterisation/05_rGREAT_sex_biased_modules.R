# ============================================================
# GREAT enrichment of sex-biased co-methylation modules
#
# Comparison: Age2 vs Age4
# Reference modules: Male Age1 network
#
# Enrichment categories:
#   - Genes
#   - Gene Ontology (MF and BP)
#   - Mouse and human phenotypes
#
# Significance: BH-adjusted hypergeometric p < 0.01
# Retain up to 30 enriched terms per category
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(rGREAT)
})

options(scipen = 999)

# ------------------------------------------------------------
# 1. Configure directories and parameters
# ------------------------------------------------------------

input_path <- "/results/02_wgcna_within_sex"
output_path <- "/results/02_wgcna_within_sex"
bgProbes_path <- "/data/gsDNAm"

modGrp <- "change_Female"
mod_ref <- "M"

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Load background CpGs
# ------------------------------------------------------------

opi <- fread(
  file.path(bgProbes_path, "GS_sub_qc1_h19.opi")
)

# GREAT input: chromosome, start, end, probe ID
bg.beta <- opi[, .(V1, V3, V3, V2)]

# ------------------------------------------------------------
# 3. Load sex-biased modules from Age2–Age4 analysis
# ------------------------------------------------------------

# Output from 03_identify_sex_biased_modules.R
dt_summary <- fread(
  file.path(
    input_path,
    paste0(
      "01a_sex_specific_mods_modref_",
      mod_ref, ".tsv"
    )
  )
)

# Select modules showing sex-biased midlife changes
selected_modules <- unique(
  dt_summary[
    sexChangeGroup2v4 == modGrp,
    moduleColor
  ]
)

# Load reference module membership to recover CpG IDs
dt_mod_labs <- fread(
  file.path(
    input_path,
    paste0(
      "precorr_gs_qc_chr_all_refage1_MinMod20_",
      mod_ref, ".modLabs"
    )
  )
)

dt <- dt_mod_labs[
  moduleColors %in% selected_modules
]

if (nrow(dt) == 0L) {
  stop("No CpGs found for the selected sex-biased modules.")
}

# ------------------------------------------------------------
# 4. Helper: retain significant GREAT enrichment results
# ------------------------------------------------------------

extract_enrichment <- function(tbl, group_name) {

  res <- as.data.table(tbl)

  res <- res[
    !is.na(Hyper_Adjp_BH) & Hyper_Adjp_BH < 0.01
  ]

  setorder(res, Hyper_Raw_PValue)

  res <- head(res, 30)

  res[, Group := group_name]

  return(res)
}

# ------------------------------------------------------------
# 5. Run GREAT for each selected module
# ------------------------------------------------------------

for (mod in unique(dt$moduleColors)) {

  cat("\nAnnotating module:", mod, "\n")

  dt_mod <- dt[moduleColors == mod]

  # Retrieve CpG genomic positions from background
  bed_mod <- bg.beta[
    V2 %in% dt_mod$probe_id
  ][
    order(V1)
  ]

  if (nrow(bed_mod) == 0L) {
    warning("No background CpGs found for module: ", mod)
    next
  }

  # Submit GREAT job
  job <- submitGreatJob(
    bed_mod,
    bg = bg.beta,
    request_interval = 30,
    help = FALSE
  )

  # Retrieve enrichment categories
  tb.go <- getEnrichmentTables(job, category = "GO")
  tb.pheno <- getEnrichmentTables(job, category = "Phenotype")
  tb.genes <- getEnrichmentTables(job, category = "Genes")

  # ----------------------------------------------------------
  # 6. Enriched genes
  # ----------------------------------------------------------

  Genes <- extract_enrichment(
    tb.genes[[1]],
    "EnsGenes"
  )

  Genes[, sigGroup := modGrp]
  Genes[, moduleColor := mod]

  # ----------------------------------------------------------
  # 7. Gene Ontology enrichment
  # ----------------------------------------------------------

  MF <- extract_enrichment(tb.go[[1]], "MF")
  BP <- extract_enrichment(tb.go[[2]], "BP")

  GO_results <- rbindlist(
    list(MF, BP),
    use.names = TRUE,
    fill = TRUE
  )

  GO_results <- GO_results[
    , .(
      ID,
      name,
      Hyper_Fold_Enrichment,
      Hyper_Raw_PValue,
      Hyper_Adjp_BH,
      Group
    )
  ]

  GO_results[, sigGroup := modGrp]
  GO_results[, moduleColor := mod]

  # ----------------------------------------------------------
  # 8. Phenotype enrichment
  # ----------------------------------------------------------

  mp <- extract_enrichment(
    tb.pheno[[1]],
    "Mouse pheno"
  )

  mp_sKO <- extract_enrichment(
    tb.pheno[[2]],
    "Mouse pheno sKO"
  )

  hp <- extract_enrichment(
    tb.pheno[[3]],
    "Human pheno"
  )

  pheno_results <- rbindlist(
    list(mp, mp_sKO, hp),
    use.names = TRUE,
    fill = TRUE
  )

  pheno_results <- pheno_results[
    , .(
      ID,
      name,
      Hyper_Fold_Enrichment,
      Hyper_Raw_PValue,
      Hyper_Adjp_BH,
      Group
    )
  ]

  pheno_results[, sigGroup := modGrp]
  pheno_results[, moduleColor := mod]

  # ----------------------------------------------------------
  # 9. Save results
  # ----------------------------------------------------------

  output_prefix <- file.path(
    output_path,
    paste0(
      "01a_", modGrp, "_", mod,
      "modref_", mod_ref, "_rGREAT"
    )
  )

  fwrite(
    Genes,
    paste0(output_prefix, "_ensgenes.out"),
    sep = "\t",
    quote = FALSE
  )

  fwrite(
    GO_results,
    paste0(output_prefix, "_GO.out"),
    sep = "\t",
    quote = FALSE
  )

  fwrite(
    pheno_results,
    paste0(output_prefix, "_pheno.out"),
    sep = "\t",
    quote = FALSE
  )

  cat("Completed module:", mod, "\n")
}

cat("\nGREAT enrichment analysis completed.\n")

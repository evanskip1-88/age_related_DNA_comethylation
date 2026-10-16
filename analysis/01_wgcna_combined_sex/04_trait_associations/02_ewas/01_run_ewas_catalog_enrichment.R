
# =============================================================================
# Purpose: Test gene-level enrichment of EWAS Catalog associations among
#          age-stable and age-variable DNA co-methylation groups.
#
# Inputs:  WGCNA gene lists, GS background rGREAT annotations, and
#          EWAS Catalog top-1,000 CpG rGREAT annotations.
# Outputs: EWAS enrichment statistics and selected top associations.
# =============================================================================

library(data.table)
library(writexl)

# -------------------------------------------------------------------------
# 1. Define input and output paths
# -------------------------------------------------------------------------

wgcna_dir <- "/results/01_wgcna_combined_sex"
annotation_dir <- "/data/gsDNAm"
ewas_dir <- "/data/ewas_catalog"
output_dir <- "/results/01_wgcna_combined_sex"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# 2. Load module genes and background annotations
# -------------------------------------------------------------------------

heteroskid_wgcna_gene_list <- readRDS(
  file.path(wgcna_dir, "heteroskid_wgcna_gene_list.rds")
)

opi_rgreat <- fread(
  file.path(
    annotation_dir,
    "GS_sub_qc1_gene_cpgs_assoc_rGREAT_h19_all_CpGs.opi"
  )
)

# Background universe: unique rGREAT-annotated genes
whole_bg_genes <- unique(
  na.omit(unlist(strsplit(opi_rgreat$rgreat_genes, ";")))
)

# -------------------------------------------------------------------------
# 3. Load EWAS Catalog top-1,000 CpG annotations
# -------------------------------------------------------------------------

dt_topN_rgreat <- fread(
  file.path(
    ewas_dir,
    "ewas_catalog_top1k_gene_cpgs_assoc_basalExt_h38.out"
  )
)

# Expand semicolon-separated gene annotations
dt_topN_sub_rgreat <- unique(
  dt_topN_rgreat[
    ,
    .(
      rgreat_genes_long = unlist(strsplit(annotated_genes, ";"))
    ),
    by = setdiff(names(dt_topN_rgreat), "annotated_genes")
  ]
)

# Remove redundant genomic annotation columns
dt_topN_sub_rgreat[
  ,
  c("seqnames.y", "start.y") := NULL
]

# -------------------------------------------------------------------------
# 4. Retain EWAS studies with more than 100 annotated genes
# -------------------------------------------------------------------------

gene_per_study <- unique(
  dt_topN_sub_rgreat[, .(study, rgreat_genes_long)]
)[
  ,
  .N,
  by = study
]

eligible_studies <- gene_per_study[N > 100, study]

dt_topN_sub_rgreat <- dt_topN_sub_rgreat[
  study %in% eligible_studies
]

# Create a gene list for each EWAS study
ewas_study_hit_list <- split(
  dt_topN_sub_rgreat$rgreat_genes_long,
  dt_topN_sub_rgreat$study
)

# -------------------------------------------------------------------------
# 5. Hypergeometric enrichment function
# -------------------------------------------------------------------------

perform_enrichment_fun <- function(ewas_study_hit_list, module_genes) {

  N <- length(whole_bg_genes)
  K <- length(module_genes)

  res_list <- vector("list", length(ewas_study_hit_list))

  i <- 1

  for (studySet in names(ewas_study_hit_list)) {

    study_hits <- unique(ewas_study_hit_list[[studySet]])

    n <- length(study_hits)

    overlaps <- intersect(study_hits, module_genes)
    k <- length(overlaps)

    # Contingency table
    a <- k
    b <- K - k
    c <- n - k
    d <- N - K - (n - k)

    # Odds ratio and approximate 95% confidence interval
    if (a == 0 || b == 0 || c == 0 || d == 0) {

      OR <- NA_real_
      CI_lower <- NA_real_
      CI_upper <- NA_real_

    } else {

      OR <- (a * d) / (b * c)

      SE <- sqrt(1/a + 1/b + 1/c + 1/d)

      logOR <- log(OR)

      CI_lower <- exp(logOR - 1.96 * SE)
      CI_upper <- exp(logOR + 1.96 * SE)
    }

    # One-sided hypergeometric enrichment test
    pval <- phyper(
      k - 1,
      K,
      N - K,
      n,
      lower.tail = FALSE
    )

    res_list[[i]] <- data.table(
      ewas_study = studySet,
      study_hits = K,
      ewas_hits = n,
      overlap = k,
      pvalue = pval,
      OR = OR,
      CI_lower = CI_lower,
      CI_upper = CI_upper
    )

    i <- i + 1
  }

  res <- rbindlist(res_list)

  # Benjamini-Hochberg correction within each analysis group
  res[, FDR := p.adjust(pvalue, method = "fdr")]

  return(res)
}

# -------------------------------------------------------------------------
# 6. Run enrichment across all study-defined gene groups
# -------------------------------------------------------------------------

all_module_results <- list()

for (mod in names(heteroskid_wgcna_gene_list)) {

  module_genes <- heteroskid_wgcna_gene_list[[mod]]

  message("Running EWAS enrichment: ", mod)

  res_mod <- perform_enrichment_fun(
    ewas_study_hit_list = ewas_study_hit_list,
    module_genes = module_genes
  )

  res_mod[, analysis := mod]

  all_module_results[[mod]] <- res_mod
}

# Combine enrichment results
res_all <- rbindlist(all_module_results, fill = TRUE)

# Confidence interval width
res_all[, CI_range := CI_upper - CI_lower]

# -------------------------------------------------------------------------
# 7. Select top EWAS associations
# -------------------------------------------------------------------------

# Preserve original criteria:
#   - At least 100 EWAS-associated genes
#   - Nominal p-value < 0.05
#   - Confidence interval width < 5
#   - Top 20 associations per analysis

sub_top <- res_all[
  ewas_hits >= 100 &
    pvalue < 0.05 &
    CI_range < 5
][
  ,
  head(.SD[order(pvalue)], 20),
  by = analysis
]

# -------------------------------------------------------------------------
# 8. Export results
# -------------------------------------------------------------------------

# Complete enrichment results
fwrite(
  res_all,
  file.path(output_dir, "12a_ewas_catalog_all_overlaps.tsv"),
  sep = "\t",
  quote = FALSE
)

# Selected associations used for downstream visualisation
write_xlsx(
  sub_top,
  file.path(output_dir, "12a_ewas_catalog_top_overlaps.xlsx")
)

message("EWAS Catalog enrichment completed.")


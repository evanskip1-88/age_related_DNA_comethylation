# =============================================================================
# Purpose: Test enrichment of methylation-associated genes among top MAGMA
#          GWAS genes using hypergeometric tests.
# Inputs:  WGCNA gene lists, CpG-to-gene mappings, MAGMA gene annotations
#          and gene-level association results.
# Outputs: GWAS enrichment statistics and overlapping gene lists.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# Run from the repository root. Adjust paths to match the repository's data layout.
gwas_annotated_dir <- "results/trait_associations/gwas/gene_annotation"
gwas_gene_results_dir <- "results/trait_associations/gwas/magma_gene_analysis"
cpg_gene_mapping_file <- "data/processed/gwas/magma/opi_build37_cpg5kb_up_1kdown_window_merged_genes.txt"
wgcna_gene_list_file <- "data/processed/wgcna/heteroskid_wgcna_gene_list.rds"
output_dir <- "results/trait_associations/gwas/gene_enrichment"

excludeMHC <- "yes"
gwasHitpercentCutoff <- 0.05  # Original analysis also considered 0.025

if (!excludeMHC %in% c("yes", "no")) {
  stop("excludeMHC must be either 'yes' or 'no'")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Discover MAGMA annotations; preserve the original excluded trait.
file_list <- list.files(gwas_annotated_dir, pattern = "\\.genes\\.annot$")
file_list <- file_list[file_list != "GIANT-heightYengo.genes.annot"]
if (length(file_list) == 0L) stop("No MAGMA gene annotation files found")

heteroskid_wgcna_gene_list <- readRDS(wgcna_gene_list_file)
dt_annotate_opi <- fread(cpg_gene_mapping_file)

perform_gwas_enrich <- function(analysis) {
  res_pvals_lst <- list()
  res_top_gene_hits <- list()

  for (i in seq_along(file_list)) {
    fname <- file_list[i]
    filePrefix <- sub("\\.genes\\.annot$", "", fname)
    cat("Analysing", filePrefix, "(", i, "/", length(file_list), ")\n")

    # Genes associated with the selected methylation group.
    gene_hit_list <- heteroskid_wgcna_gene_list[[analysis]]
    dt_study_hit_sub <- dt_annotate_opi[gene_build37 %in% gene_hit_list]
    uniq_gene_hits <- unique(dt_study_hit_sub$geneID)

    # MAGMA gene annotations and their overlap with the CpG annotation.
    # MAGMA .genes.annot files contain gene ID, position, and SNP IDs.
    annotation_file <- file.path(gwas_annotated_dir, fname)
    dt_annoted_gwas <- fread(annotation_file, header = FALSE, select = 1:3)
    setnames(dt_annoted_gwas, c("geneID", "position", "snp_id"))
    dt_opi_gwas_gene_overlap <- dt_annoted_gwas[geneID %in% dt_annotate_opi$geneID]

    # Select the top fraction of MAGMA genes before excluding MHC genes,
    # exactly as in the original analysis.
    gwas_genes_out <- fread(file.path(gwas_gene_results_dir,
                                     paste0(filePrefix, ".genes.out")))
    sorted_data <- gwas_genes_out[order(P)]
    top_rows <- ceiling(gwasHitpercentCutoff * nrow(sorted_data))
    top_N_percent_gwas_genes <- sorted_data[seq_len(top_rows)]

    if (excludeMHC == "yes") {
      top_N_percent_gwas_genes <- top_N_percent_gwas_genes[
        !(CHR == 6 & START <= 34000000 & STOP >= 25000000)
      ]
    }

    top_N_percent_gwas_uniq_genes <- unique(top_N_percent_gwas_genes$GENE)
    overlapping_genes <- intersect(top_N_percent_gwas_uniq_genes,
                                   uniq_gene_hits)

    # Export the identities of overlapping genes.
    sub_top_overlap_genes <- top_N_percent_gwas_genes[GENE %in% uniq_gene_hits]
    sub_top_overlap_genes_m <- unique(merge(
      dt_annotate_opi[, .(gene_opi, geneID)],
      sub_top_overlap_genes,
      by.x = "geneID", by.y = "GENE"
    ))
    sub_top_overlap_genes_m[, `:=`(TraitID = filePrefix, analysis = analysis)]
    res_top_gene_hits[[filePrefix]] <- sub_top_overlap_genes_m

    # Preserve original hypergeometric universe/count definitions.
    # Audit separately: N counts annotation rows, not unique gene IDs.
    N <- length(dt_opi_gwas_gene_overlap$geneID)
    K <- length(top_N_percent_gwas_uniq_genes)
    n <- length(uniq_gene_hits)
    k <- length(overlapping_genes)

    if (k == 0 || K == k || n == k || (N - K - (n - k)) == 0) {
      OR <- NA_real_
      CI_lower <- NA_real_
      CI_upper <- NA_real_
    } else {
      OR <- (k / (K - k)) / ((n - k) / (N - K - (n - k)))
      SE_log_OR <- sqrt(1/k + 1/(n-k) + 1/(K-k) +
                        1/(N-K-(n-k)))
      log_OR <- log(OR)
      CI_lower <- exp(log_OR - 1.96 * SE_log_OR)
      CI_upper <- exp(log_OR + 1.96 * SE_log_OR)
    }

    res_pvals_lst[[fname]] <- data.frame(
      analysis = analysis,
      trait = filePrefix,
      bg_genes = N,
      topgwas_genes = K,
      study_hit_genes = n,
      overlap = k,
      gwasHitpercentCutoff = gwasHitpercentCutoff,
      pvalue = phyper(k - 1, K, N - K, n, lower.tail = FALSE),
      OR = OR,
      CI_lower = CI_lower,
      CI_upper = CI_upper
    )
  }

  list(
    pvals = rbindlist(res_pvals_lst, fill = TRUE),
    genes = rbindlist(res_top_gene_hits, fill = TRUE)
  )
}

res_final_list <- list()
res_top_gene_hits_final <- list()

for (list_name in names(heteroskid_wgcna_gene_list)) {
  cat("Analysis for:", list_name, "\n")
  res_enrich <- perform_gwas_enrich(list_name)
  res_final_list[[list_name]] <- res_enrich$pvals
  res_top_gene_hits_final[[list_name]] <- res_enrich$genes
}

res_gwas_enrich_final <- rbindlist(res_final_list, fill = TRUE)
res_gwas_enrich_final[, MHC := excludeMHC]
res_gwas_gene_overlap <- rbindlist(res_top_gene_hits_final, fill = TRUE)

output_suffix <- paste0("_", gwasHitpercentCutoff * 100, "percent.tsv")
fwrite(
  res_gwas_enrich_final[order(pvalue)],
  file.path(output_dir, paste0("01c_mhc_", excludeMHC,
                               "_gwas_gene_overlap", output_suffix)),
  sep = "\t", quote = FALSE
)
fwrite(
  res_gwas_gene_overlap,
  file.path(output_dir, paste0("01c_mhc_", excludeMHC,
                               "_gwas_overlapping_genes", output_suffix)),
  sep = "\t", quote = FALSE
)

cat("Done!\n")


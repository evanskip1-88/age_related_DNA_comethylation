#==============================================================================
# Transcription factor binding-site (TFBS) enrichment
#
# Tests enrichment of age-stable and age-variable CpG sets in ENCODE
# transcription factor binding regions (hg19).
#
# Statistical test: one-sided hypergeometric test
#==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

#------------------------------------------------------------------------------
# Paths
#------------------------------------------------------------------------------

output_path <- "/results/01_wgcna_combined_sex"
wgcna_cpgs <- "/results/01_wgcna_combined_sex"
opi_bg_dat_path <- "/data/gsDNAm"
tfbs_path <- "/data/encode_tfbs"

#------------------------------------------------------------------------------
# Load background CpGs (hg19)
#------------------------------------------------------------------------------

background_cpgs <- fread(
  file.path(opi_bg_dat_path, "GS_sub_qc1.opi")
)

setnames(
  background_cpgs,
  c("chr_num", "probe_id", "position", "gene", "strand")
)

background_cpgs[
  ,
  chr := paste0("chr", chr_num)
]

background_cpgs <- background_cpgs[
  !is.na(position)
]

#------------------------------------------------------------------------------
# Load CpG sets
#------------------------------------------------------------------------------

cpg_sets <- readRDS(
  file.path(wgcna_cpgs, "heteroskid_wgcna_cpg_list.rds")
)

#------------------------------------------------------------------------------
# Load ENCODE TFBS annotations
#------------------------------------------------------------------------------

tfbs <- fread(
  file.path(tfbs_path, "encRegTfbsClusteredWithCells.hg19.bed")
)

setnames(
  tfbs,
  c("chr", "start", "end", "tf", "score", "cell_type")
)

#------------------------------------------------------------------------------
# Identify CpGs overlapping TFBS regions
#------------------------------------------------------------------------------

tfbs_gr <- GRanges(
  seqnames = tfbs$chr,
  ranges = IRanges(
    start = tfbs$start,
    end = tfbs$end
  ),
  TF = tfbs$tf,
  CellType = tfbs$cell_type
)

cpg_gr <- GRanges(
  seqnames = background_cpgs$chr,
  ranges = IRanges(
    start = background_cpgs$position,
    end = background_cpgs$position
  ),
  probe_id = background_cpgs$probe_id,
  gene = background_cpgs$gene
)

cat("Identifying CpG overlaps with ENCODE TFBS...\n")

overlaps <- findOverlaps(cpg_gr, tfbs_gr)

overlapping_cpgs <- data.table(
  probe_id = cpg_gr$probe_id[queryHits(overlaps)],
  gene = cpg_gr$gene[queryHits(overlaps)],
  TF = tfbs_gr$TF[subjectHits(overlaps)],
  CellType = tfbs_gr$CellType[subjectHits(overlaps)]
)

#------------------------------------------------------------------------------
# Filter transcription factors by number of associated genes
#------------------------------------------------------------------------------

filtered_cpgs <- overlapping_cpgs[
  !is.na(gene)
]

filtered_cpgs <- filtered_cpgs[
  ,
  .(gene = unlist(strsplit(gene, ";"))),
  by = .(probe_id, TF)
]

filtered_cpgs <- unique(filtered_cpgs)

tf_gene_counts <- filtered_cpgs[
  ,
  .(num_genes = uniqueN(gene)),
  by = TF
]

valid_tfs <- tf_gene_counts[
  num_genes >= 5 & num_genes <= 2000,
  TF
]

filtered_cpgs_fnl <- filtered_cpgs[
  TF %in% valid_tfs
]

cat("TFs retained:", uniqueN(filtered_cpgs_fnl$TF), "\n")

#------------------------------------------------------------------------------
# TFBS enrichment function
#------------------------------------------------------------------------------

tf_enrichment <- function(tfbs_data, hit_data, sig_group, analysis) {

  unique_tfs <- unique(tfbs_data$TF)
  results <- vector("list", length(unique_tfs))

  N <- nrow(background_cpgs)
  study_hit_list <- unique(hit_data)
  n <- length(study_hit_list)

  for (i in seq_along(unique_tfs)) {

    tf <- unique_tfs[i]

    pub_cpg_list <- unique(
      tfbs_data[TF == tf, probe_id]
    )

    K <- length(pub_cpg_list)

    k <- length(
      intersect(study_hit_list, pub_cpg_list)
    )

    # One-sided hypergeometric enrichment
    pval <- phyper(
      k - 1,
      K,
      N - K,
      n,
      lower.tail = FALSE
    )

    # Contingency-table counts
    a <- k
    b <- n - k
    c <- K - k
    d <- N - K - n + k

    if (all(c(a, b, c, d) > 0)) {

      OR <- (a * d) / (b * c)

      SE_log_OR <- sqrt(
        1 / a + 1 / b + 1 / c + 1 / d
      )

      CI_lower <- exp(log(OR) - 1.96 * SE_log_OR)
      CI_upper <- exp(log(OR) + 1.96 * SE_log_OR)

    } else {

      OR <- NA_real_
      CI_lower <- NA_real_
      CI_upper <- NA_real_

    }

    results[[i]] <- data.table(
      TF = tf,
      SigGroup = sig_group,
      Analysis = analysis,
      pub_cpg_hits = K,
      study_cpg_hits = n,
      overlap = k,
      pval = pval,
      OR = OR,
      CI_lower = CI_lower,
      CI_upper = CI_upper
    )
  }

  rbindlist(results)
}

#------------------------------------------------------------------------------
# Define analysis groups
#------------------------------------------------------------------------------

analysis_list <- list(
  list(
    hit_data = cpg_sets[["Heteroskid_Stable"]],
    sig_group = "Stable",
    analysis = "Heteroskedasticity",
    output_name = "het_Stable"
  ),
  list(
    hit_data = cpg_sets[["Heteroskid_Changing"]],
    sig_group = "Changing",
    analysis = "Heteroskedasticity",
    output_name = "het_Changing"
  ),
  list(
    hit_data = cpg_sets[["WGCNA_Stable"]],
    sig_group = "Stable",
    analysis = "WGCNA",
    output_name = "wgcna_Stable"
  ),
  list(
    hit_data = cpg_sets[["WGCNA_Changing"]],
    sig_group = "Changing",
    analysis = "WGCNA",
    output_name = "wgcna_Changing"
  )
)

#------------------------------------------------------------------------------
# Run enrichment
#------------------------------------------------------------------------------

all_results <- vector("list", length(analysis_list))

for (i in seq_along(analysis_list)) {

  params <- analysis_list[[i]]

  cat("Running:", params$output_name, "\n")

  result <- tf_enrichment(
    tfbs_data = filtered_cpgs_fnl,
    hit_data = params$hit_data,
    sig_group = params$sig_group,
    analysis = params$analysis
  )

  result[
    ,
    Dataset := params$output_name
  ]

  all_results[[i]] <- result
}

#------------------------------------------------------------------------------
# Save results
#------------------------------------------------------------------------------

combined_results <- rbindlist(
  all_results,
  use.names = TRUE,
  fill = TRUE
)

setorder(combined_results, -CI_lower)

dir.create(
  output_path,
  recursive = TRUE,
  showWarnings = FALSE
)

output_file <- file.path(
  output_path,
  "01a_tfbs_enrich_no_perm.tsv"
)

fwrite(
  combined_results,
  output_file,
  sep = "\t"
)

cat("TFBS enrichment completed.\n")
cat("Results saved to:", output_file, "\n")

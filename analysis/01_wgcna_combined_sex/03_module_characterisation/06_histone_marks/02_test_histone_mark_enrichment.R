#==============================================================================
# Histone-mark enrichment analysis
#
# Tests enrichment of age-stable and age-variable CpG sets in published
# histone-mark regions across selected tissues and cell types.
#
# Statistical test: one-sided hypergeometric enrichment
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
histone_data <- "/data/RoadMapHistoneMark"
opi_bg_dat_path <- "/data/gsDNAm"

#------------------------------------------------------------------------------
# Load background CpGs (hg19)
#------------------------------------------------------------------------------

background_cpgs <- fread(
  file.path(
    opi_bg_dat_path,
    "GS_sub_qc1_gene_cpgs_assoc_rGREAT_h19_all_CpGs.opi"
  )
)

setnames(
  background_cpgs,
  c(2, 4, 7),
  c("chr", "gene", "orig_gene")
)

#------------------------------------------------------------------------------
# Load selected tissues and CpG sets
#------------------------------------------------------------------------------

selected_tissues <- fread(
  file.path(histone_data, "histone_marks_final_meta_sel_tiss.txt")
)

cpg_sets <- readRDS(
  file.path(wgcna_cpgs, "heteroskid_wgcna_cpg_list.rds")
)

#------------------------------------------------------------------------------
# Construct background genomic ranges
#------------------------------------------------------------------------------

background_gr <- GRanges(
  seqnames = background_cpgs$chr,
  ranges = IRanges(
    start = background_cpgs$position,
    end = background_cpgs$position
  ),
  probe_id = background_cpgs$probe_id
)

background_ids <- unique(background_cpgs$probe_id)
N <- length(background_ids)

#------------------------------------------------------------------------------
# Histone-mark enrichment function
#------------------------------------------------------------------------------

perform_histone_enrichment <- function(
    study_cpgs,
    analysis,
    sig_group,
    top_n = 10000
) {

  results <- vector("list", nrow(selected_tissues))

  study_cpgs <- intersect(
    unique(study_cpgs),
    background_ids
  )

  n <- length(study_cpgs)

  for (i in seq_len(nrow(selected_tissues))) {

    histone_meta <- selected_tissues[i]

    fname <- as.character(selected_tissues[[1]][i])

    cat(
      "Analysing histone-mark file:",
      fname,
      "\n"
    )

    # Load published histone-mark regions
    hist_marks <- fread(
      file.path(histone_data, fname)
    )

    # Retain first 10,000 regions
    hist_marks_top <- head(hist_marks, top_n)

    hist_marks_gr <- GRanges(
      seqnames = hist_marks_top$V1,
      ranges = IRanges(
        start = hist_marks_top$V2,
        end = hist_marks_top$V3
      )
    )

    # Identify overlapping background CpGs
    overlaps <- findOverlaps(
      background_gr,
      hist_marks_gr
    )

    overlapping_ids <- unique(
      background_gr$probe_id[queryHits(overlaps)]
    )

    K <- length(overlapping_ids)

    # Study CpGs overlapping histone-mark regions
    k <- length(
      intersect(study_cpgs, overlapping_ids)
    )

    # Hypergeometric enrichment
    p_value <- phyper(
      k - 1,
      K,
      N - K,
      n,
      lower.tail = FALSE
    )

    # Contingency table
    a <- k
    b <- n - k
    c <- K - k
    d <- N - K - n + k

    if (all(c(a, b, c, d) > 0)) {

      OR <- (a * d) / (b * c)

      SE_log_OR <- sqrt(
        1 / a +
        1 / b +
        1 / c +
        1 / d
      )

      CI_lower <- exp(log(OR) - 1.96 * SE_log_OR)
      CI_upper <- exp(log(OR) + 1.96 * SE_log_OR)

    } else {

      OR <- NA_real_
      CI_lower <- NA_real_
      CI_upper <- NA_real_

    }

    results[[i]] <- cbind(
      histone_meta,
      data.table(
        Analysis = analysis,
        sigGroup = sig_group,
        cpg_hm_whole = K,
        cpgs_study_hits = n,
        overlap = k,
        p_value = p_value,
        OR = OR,
        CI_lower = CI_lower,
        CI_upper = CI_upper
      )
    )
  }

  rbindlist(results, use.names = TRUE, fill = TRUE)
}

#------------------------------------------------------------------------------
# Run enrichment for four CpG sets
#------------------------------------------------------------------------------

enrichment_sets <- list(
  WGCNA_Stable = list(
    data = cpg_sets[["WGCNA_Stable"]],
    Analysis = "WGCNA",
    sigGroup = "Stable"
  ),

  WGCNA_Changing = list(
    data = cpg_sets[["WGCNA_Changing"]],
    Analysis = "WGCNA",
    sigGroup = "Changing"
  ),

  Heteroskid_Stable = list(
    data = cpg_sets[["Heteroskid_Stable"]],
    Analysis = "Heteroskedasticity",
    sigGroup = "Stable"
  ),

  Heteroskid_Changing = list(
    data = cpg_sets[["Heteroskid_Changing"]],
    Analysis = "Heteroskedasticity",
    sigGroup = "Changing"
  )
)

all_results <- list()

for (set_name in names(enrichment_sets)) {

  cat("Running enrichment for:", set_name, "\n")

  set_info <- enrichment_sets[[set_name]]

  all_results[[set_name]] <- perform_histone_enrichment(
    study_cpgs = set_info$data,
    analysis = set_info$Analysis,
    sig_group = set_info$sigGroup,
    top_n = 10000
  )
}

#------------------------------------------------------------------------------
# Save results
#------------------------------------------------------------------------------

final_results <- rbindlist(
  all_results,
  use.names = TRUE,
  fill = TRUE
)

dir.create(
  output_path,
  recursive = TRUE,
  showWarnings = FALSE
)

fwrite(
  final_results,
  file.path(
    output_path,
    "01a_histone_enrichment_res_combined_mods_n_het.tsv"
  ),
  sep = "\t"
)

cat("Histone-mark enrichment completed.\n")

#==============================================================================
# Chromatin-state enrichment analysis
#
# Tests enrichment of age-stable and age-variable CpG sets in published
# human ChromHMM chromatin states (hg19).
#
# CpG sets:
#   WGCNA_Stable       - aSCMs
#   WGCNA_Changing     - aVCMs
#   Heteroskid_Stable  - aSCPs
#   Heteroskid_Changing- aVCPs
#
# Statistical test: one-sided hypergeometric enrichment
#==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

rm(list = ls())

#------------------------------------------------------------------------------
# Paths
#------------------------------------------------------------------------------

output_path <- "/results/01_wgcna_combined_sex"
pub_dat_path <- "/data/vuChromStates"
# Published chrom states can be downloaded from Vu & Ernest 2022:
# https://public.hoffman2.idre.ucla.edu/ernst/2k9rs/full_stack/full_stack_annotation_public_release/hg19/
cpg_sets_path <- "/results/01_wgcna_combined_sex"
opi_bg_dat_path <- "/data/gsDNAm"

#------------------------------------------------------------------------------
# Load background CpGs (hg19)
#------------------------------------------------------------------------------

opi <- fread(
  file.path(opi_bg_dat_path, "GS_sub_qc1.opi")
)

setnames(
  opi,
  c("chr_num", "probe_id", "position", "gene", "strand")
)

opi <- opi[!is.na(position)]

opi[, chr := paste0("chr", chr_num)]

# Restrict to unique background CpGs
opi <- unique(
  opi[, .(chr, position, probe_id)],
  by = "probe_id"
)

background_ids <- unique(opi$probe_id)
genome_size <- length(background_ids)

#------------------------------------------------------------------------------
# Load stable and changing CpG sets
#------------------------------------------------------------------------------

cpg_sets <- readRDS(
  file.path(cpg_sets_path, "heteroskid_wgcna_cpg_list.rds")
)

#------------------------------------------------------------------------------
# Load published chromatin-state segmentation (hg19)
#------------------------------------------------------------------------------

chrom_states <- fread(
  file.path(pub_dat_path, "hg19_genome_100_segments.bed")
)

setnames(
  chrom_states,
  c("chr", "binStart", "binEnd", "chromState")
)

# BED coordinates are 0-based, half-open.
# Convert to 1-based coordinates for GRanges.

chrom_gr <- GRanges(
  seqnames = chrom_states$chr,
  ranges = IRanges(
    start = chrom_states$binStart + 1L,
    end = chrom_states$binEnd
  ),
  chromState = chrom_states$chromState
)

#------------------------------------------------------------------------------
# Identify background CpGs overlapping chromatin states
#------------------------------------------------------------------------------

opi_gr <- GRanges(
  seqnames = opi$chr,
  ranges = IRanges(
    start = opi$position,
    end = opi$position
  ),
  probe_id = opi$probe_id
)

overlaps <- findOverlaps(opi_gr, chrom_gr)

overlap_table <- unique(
  data.table(
    probe_id = opi_gr$probe_id[queryHits(overlaps)],
    chrom_state = chrom_gr$chromState[subjectHits(overlaps)]
  )
)

chrom_state_list <- unique(overlap_table$chrom_state)

#------------------------------------------------------------------------------
# Hypergeometric enrichment
#------------------------------------------------------------------------------

perform_chromatin_enrichment <- function(
    study_cpgs,
    chrom_state,
    background_ids,
    overlap_table
) {

  # CpGs assigned to this chromatin state
  state_cpgs <- unique(
    overlap_table[
      chrom_state == ..chrom_state,
      probe_id
    ]
  )

  # Study CpGs represented in the background
  study_cpgs <- intersect(
    unique(study_cpgs),
    background_ids
  )

  N <- length(background_ids)
  K <- length(state_cpgs)
  n <- length(study_cpgs)
  k <- length(intersect(study_cpgs, state_cpgs))

  # One-sided hypergeometric enrichment
  p_value <- phyper(
    k - 1,
    K,
    N - K,
    n,
    lower.tail = FALSE
  )

  # 2 x 2 contingency table
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

  data.table(
    chrom_state = chrom_state,
    N_cpgs_chromS = K,
    N_study_cpgs = n,
    overlap = k,
    pvalue = p_value,
    OR = OR,
    CI_lower = CI_lower,
    CI_upper = CI_upper
  )
}

#------------------------------------------------------------------------------
# Run enrichment for all CpG sets and chromatin states
#------------------------------------------------------------------------------

results <- list()

for (analysis_name in names(cpg_sets)) {

  cat("Analysing:", analysis_name, "\n")

  study_cpgs <- cpg_sets[[analysis_name]]

  for (state_name in chrom_state_list) {

    result <- perform_chromatin_enrichment(
      study_cpgs = study_cpgs,
      chrom_state = state_name,
      background_ids = background_ids,
      overlap_table = overlap_table
    )

    result[, analysis := analysis_name]

    results[[length(results) + 1L]] <- result
  }
}

res_final <- rbindlist(results)

setcolorder(
  res_final,
  c(
    "analysis",
    "chrom_state",
    "N_cpgs_chromS",
    "N_study_cpgs",
    "overlap",
    "pvalue",
    "OR",
    "CI_lower",
    "CI_upper"
  )
)

#------------------------------------------------------------------------------
# Save results
#------------------------------------------------------------------------------

dir.create(
  output_path,
  recursive = TRUE,
  showWarnings = FALSE
)

fwrite(
  res_final,
  file.path(
    output_path,
    "01a_chrom_state_enrich_combined_mods_n_het.tsv"
  ),
  sep = "\t"
)

cat("Chromatin-state enrichment completed.\n")


#!/usr/bin/env Rscript

# =============================================================================
# TAD boundary and chromatin loop enrichment
#
# Reference: Yue Lab 3D Genome
# https://3dgenome.fsm.northwestern.edu/publications.html
#
# Inputs:
#   - CpG annotation (hg19)
#   - Selected co-methylation CpG sets
#   - TAD coordinates (*TADs.txt)
#   - Chromatin loop coordinates (*.merged.loops)
#
# Outputs:
#   01a_chrom_loop_enrich.tsv
#   01b_tad_boundery_enrich.tsv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

# -----------------------------------------------------------------------------
# 1. Paths
# -----------------------------------------------------------------------------

cpgs_dir <- "/results/01_wgcna_combined_sex"
opi_dir <- "/data/gsDNAm"
loop_dir <- "/data/chromatin_loops"
tad_dir <- "/data/tads"
output_dir <- "/results/01_wgcna_combined_sex"

opi_file <- file.path(
  opi_dir,
  "GS_sub_qc1_gene_cpgs_assoc_rGREAT_h19_all_CpGs.opi"
)

cpg_list_file <- file.path(
  cpgs_dir,
  "heteroskid_wgcna_cpg_list.rds"
)

# -----------------------------------------------------------------------------
# 2. Load CpG background and selected sets
# -----------------------------------------------------------------------------

opi <- fread(opi_file)

setnames(
  opi,
  old = names(opi)[c(2, 4, 7)],
  new = c("chr", "gene", "orig_gene")
)

# Original background definition
cpgSize <- uniqueN(opi[position != "NA", probe_id])

heteroskid_wgcna_cpgs_list <- readRDS(cpg_list_file)

# Remove missing genomic positions before creating GRanges
opi <- opi[!is.na(position) & position != "NA"]

opi[, position := as.integer(position)]

opi_gr <- GRanges(
  seqnames = opi$chr,
  ranges = IRanges(
    start = opi$position,
    end = opi$position
  ),
  probe_id = opi$probe_id
)

# -----------------------------------------------------------------------------
# 3. Shared enrichment calculation
# -----------------------------------------------------------------------------

calculate_enrichment <- function(
    overlap_cpgs,
    study_cpgs,
    background_size
) {

  # Preserve original counting definitions
  N <- background_size
  K <- length(overlap_cpgs)
  n <- length(study_cpgs)
  k <- sum(overlap_cpgs %in% study_cpgs)

  OR <- NA_real_
  CI_lower <- NA_real_
  CI_upper <- NA_real_

  if (
    k > 0 &&
    K > k &&
    n > k &&
    (N - K - (n - k)) > 0
  ) {

    OR <- (k / (K - k)) /
      ((n - k) / (N - K - (n - k)))

    SE_log_OR <- sqrt(
      1 / k +
      1 / (n - k) +
      1 / (K - k) +
      1 / (N - K - (n - k))
    )

    CI_lower <- exp(log(OR) - 1.96 * SE_log_OR)
    CI_upper <- exp(log(OR) + 1.96 * SE_log_OR)
  }

  list(
    K = K,
    n = n,
    k = k,
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

# -----------------------------------------------------------------------------
# 4. TAD boundary enrichment
# -----------------------------------------------------------------------------

tad_files <- list.files(
  tad_dir,
  pattern = "TADs\\.txt$",
  full.names = TRUE
)

tad_results <- list()

for (file in tad_files) {

  tissue <- tools::file_path_sans_ext(basename(file))
  message("Processing TADs: ", tissue)

  dt_tad <- fread(file)
  setnames(dt_tad, c("chr", "start", "end"))

  # Preserve original 100-kb boundary windows
  dt_tad_boundaries <- rbind(
    dt_tad[
      ,
      .(
        chr,
        start = pmax(start - 100000, 0),
        end = start
      )
    ],
    dt_tad[
      ,
      .(
        chr,
        start = end,
        end = end + 100000
      )
    ]
  )

  tad_gr <- GRanges(
    seqnames = dt_tad_boundaries$chr,
    ranges = IRanges(
      start = dt_tad_boundaries$start,
      end = dt_tad_boundaries$end
    )
  )

  overlaps <- findOverlaps(opi_gr, tad_gr)

  # CpGs overlapping the TAD boundary windows
  overlapping_cpgs <- opi_gr$probe_id[queryHits(overlaps)]

  for (analysis in names(heteroskid_wgcna_cpgs_list)) {

    stats <- calculate_enrichment(
      overlap_cpgs = overlapping_cpgs,
      study_cpgs = heteroskid_wgcna_cpgs_list[[analysis]],
      background_size = cpgSize
    )

    tad_results[[length(tad_results) + 1L]] <- data.table(
      analysis = analysis,
      tissue = tissue,
      tad_type = "TAD_boundary",
      all_tads = stats$K,
      cpgs_hits = stats$n,
      overlap = stats$k,
      pvalue = stats$pvalue,
      OR = stats$OR,
      CI_lower = stats$CI_lower,
      CI_upper = stats$CI_upper
    )
  }
}

res_tad_enrich_fnl <- rbindlist(
  tad_results,
  use.names = TRUE
)

fwrite(
  res_tad_enrich_fnl,
  file.path(output_dir, "01b_tad_boundery_enrich.tsv"),
  sep = "\t",
  quote = FALSE
)

# -----------------------------------------------------------------------------
# 5. Chromatin loop enrichment
# -----------------------------------------------------------------------------

loop_files <- list.files(
  loop_dir,
  pattern = "merged\\.loops$",
  full.names = TRUE
)

loop_results <- list()

for (file in loop_files) {

  tissue <- tools::file_path_sans_ext(basename(file))
  message("Processing chromatin loops: ", tissue)

  dt_loop <- fread(file)

  # Both chromatin loop anchors
  anchor1 <- dt_loop[, .(chr = V1, start = V2, end = V3)]
  anchor2 <- dt_loop[, .(chr = V4, start = V5, end = V6)]

  dt_anchors <- unique(rbind(anchor1, anchor2))

  loop_gr <- GRanges(
    seqnames = dt_anchors$chr,
    ranges = IRanges(
      start = dt_anchors$start,
      end = dt_anchors$end
    )
  )

  overlaps <- findOverlaps(opi_gr, loop_gr)

  overlapping_cpgs <- opi_gr$probe_id[queryHits(overlaps)]

  for (analysis in names(heteroskid_wgcna_cpgs_list)) {

    stats <- calculate_enrichment(
      overlap_cpgs = overlapping_cpgs,
      study_cpgs = heteroskid_wgcna_cpgs_list[[analysis]],
      background_size = cpgSize
    )

    loop_results[[length(loop_results) + 1L]] <- data.table(
      analysis = analysis,
      tissue = tissue,
      all_chromLoop = stats$K,
      cpgs_hits = stats$n,
      overlap = stats$k,
      pvalue = stats$pvalue,
      OR = stats$OR,
      CI_lower = stats$CI_lower,
      CI_upper = stats$CI_upper
    )
  }
}

res_chrom_loop_enrich_fnl <- rbindlist(
  loop_results,
  use.names = TRUE
)

fwrite(
  res_chrom_loop_enrich_fnl,
  file.path(output_dir, "01a_chrom_loop_enrich.tsv"),
  sep = "\t",
  quote = FALSE
)

message("TAD boundary and chromatin loop enrichment completed.")


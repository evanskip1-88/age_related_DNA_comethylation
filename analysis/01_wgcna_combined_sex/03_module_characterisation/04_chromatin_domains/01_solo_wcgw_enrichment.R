#==============================================================================
# Enrichment of co-methylated CpG sets in PMDs and solo-WCGW regions
#
# Tests overlap of age-stable and age-variable CpG sets with published
# human PMD and solo-WCGW genomic regions (hg19).
#
# Enrichment is assessed using a one-sided hypergeometric test, with
# odds ratios and 95% confidence intervals.
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
wgcna_cpgs <- "/results/01_wgcna_combined_sex"
solo_wcgw_path <- "data/ZhouPMD"
# Published PMD and solo-WCGW reference datasets (hg19)
# Source: https://zwdzwd.github.io/pmd
# Download the required BED files and place them in this directory.
opi_bg_dat_path <- "/data/gsDNAm"

#------------------------------------------------------------------------------
# Load Generation Scotland CpG background (hg19)
#------------------------------------------------------------------------------

background_cpgs <- fread(
  paste0(opi_bg_dat_path, "/GS_sub_qc1_h19.bed")
)

setnames(
  background_cpgs,
  c("chr", "start", "end", "probe_id")
)

#------------------------------------------------------------------------------
# Load co-methylation CpG sets
#------------------------------------------------------------------------------

cpgs_hit_list <- readRDS(
  paste0(wgcna_cpgs, "/heteroskid_wgcna_cpg_list.rds")
)

#------------------------------------------------------------------------------
# Load published PMD and solo-WCGW regions
#------------------------------------------------------------------------------

cat("Loading PMD and solo-WCGW data\n")

# Published human PMD coordinates
pmd_regions <- fread(
  paste0(solo_wcgw_path, "/PMD_coordinates_hg19.bed.gz")
)

setnames(
  pmd_regions,
  c("chr", "window_start", "window_end",
    "dnam_beta", "type", "pmd_group")
)

# Solo-WCGW sites within common PMDs
solo_common <- fread(
  paste0(solo_wcgw_path, "/solo_WCGW_inCommonPMDs_hg19.bed.gz")
)

# Genome-wide solo-WCGW sites
solo_genomewide <- fread(
  paste0(solo_wcgw_path, "/solo_WCGW_hg19.bed.gz")
)

# Group solo-WCGW sites into 100-kb windows
window_size <- 100000

solo_common[
  ,
  `:=`(
    window_start = floor(V2 / window_size) * window_size,
    window_end = floor(V2 / window_size) * window_size +
      window_size - 1,
    pmd_group = "soloinCommonPMDs"
  )
]

setnames(solo_common, "V1", "chr")

solo_genomewide[
  ,
  `:=`(
    window_start = floor(V2 / window_size) * window_size,
    window_end = floor(V2 / window_size) * window_size +
      window_size - 1,
    pmd_group = "soloGenomewide"
  )
]

setnames(solo_genomewide, "V1", "chr")

# Combine genomic regions
all_regions <- rbindlist(
  list(
    pmd_regions[
      , .(chr, window_start, window_end, pmd_group)
    ],
    solo_common[
      , .(chr, window_start, window_end, pmd_group)
    ],
    solo_genomewide[
      , .(chr, window_start, window_end, pmd_group)
    ]
  )
)

#------------------------------------------------------------------------------
# Enrichment function
#------------------------------------------------------------------------------

perform_pmd_enrich <- function(
    cpg_hit_dat,
    regions,
    background_cpgs,
    pmd_group
) {

  region_dat <- unique(
    regions[
      pmd_group == ..pmd_group,
      .(chr, window_start, window_end)
    ]
  )

  region_gr <- GRanges(
    seqnames = region_dat$chr,
    ranges = IRanges(
      start = region_dat$window_start,
      end = region_dat$window_end
    )
  )

  cpg_gr <- GRanges(
    seqnames = background_cpgs$chr,
    ranges = IRanges(
      start = background_cpgs$start,
      end = background_cpgs$end
    )
  )

  overlaps <- findOverlaps(cpg_gr, region_gr)

  overlapping_cpgs <- unique(
    background_cpgs$probe_id[queryHits(overlaps)]
  )

  background_ids <- unique(background_cpgs$probe_id)

  N <- length(background_ids)
  K <- length(overlapping_cpgs)

  results <- lapply(names(cpg_hit_dat), function(cpg_set) {

    study_cpgs <- intersect(
      unique(cpg_hit_dat[[cpg_set]]),
      background_ids
    )

    n <- length(study_cpgs)
    k <- sum(study_cpgs %in% overlapping_cpgs)

    # One-sided hypergeometric enrichment
    p_value <- phyper(
      k - 1,
      K,
      N - K,
      n,
      lower.tail = FALSE
    )

    # Contingency table counts
    a <- k
    b <- n - k
    c <- K - k
    d <- N - K - b

    OR <- (a * d) / (b * c)

    SE_log_OR <- sqrt(
      1 / a + 1 / b + 1 / c + 1 / d
    )

    CI_lower <- exp(log(OR) - 1.96 * SE_log_OR)
    CI_upper <- exp(log(OR) + 1.96 * SE_log_OR)

    data.table(
      analysis = cpg_set,
      pmd_group = pmd_group,
      dataset = "Human",
      cpg_solo_whole = K,
      cpgs_study_hits = n,
      overlap = k,
      p_value = p_value,
      OR = OR,
      CI_lower = CI_lower,
      CI_upper = CI_upper
    )
  })

  rbindlist(results)
}

#------------------------------------------------------------------------------
# Run enrichment across genomic-region groups
#------------------------------------------------------------------------------

results <- lapply(
  unique(all_regions$pmd_group),
  function(group) {

    cat("Processing:", group, "\n")

    perform_pmd_enrich(
      cpg_hit_dat = cpgs_hit_list,
      regions = all_regions,
      background_cpgs = background_cpgs,
      pmd_group = group
    )
  }
)

res_final <- rbindlist(results)

#------------------------------------------------------------------------------
# Save results
#------------------------------------------------------------------------------

fwrite(
  res_final[order(-CI_lower)],
  paste0(output_path, "/01a_PMD_enrichments_human.tsv"),
  sep = "\t",
  quote = FALSE
)

cat("PMD and solo-WCGW enrichment analysis completed.\n")

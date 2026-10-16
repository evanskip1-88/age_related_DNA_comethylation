
#!/usr/bin/env Rscript

# =============================================================================
# Lift over DNA methylation CpG coordinates from hg19 to mm10
#
# Purpose:
# Map human CpG coordinates to the mouse genome for overlap analysis
# with age-associated DARs from Patrick et al. (2024).
#
# Input:
#   GS_sub_qc1_h19.bed
#   hg19ToMm10.over.chain
#
# Output:
#   GS_cpgs_h19_mapped_to_mm10.bed
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
})

# -----------------------------------------------------------------------------
# 1. Input and output paths
# -----------------------------------------------------------------------------

opi_dir <- "/data/gsDNAm"
output_dir <- "/data/gsDNAm"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

chain_file <- file.path(opi_dir, "hg19ToMm10.over.chain")
cpg_file <- file.path(opi_dir, "GS_sub_qc1_h19.bed")

output_file <- file.path(
  output_dir,
  "GS_cpgs_h19_mapped_to_mm10.bed"
)

# -----------------------------------------------------------------------------
# 2. Load chain and CpG coordinates
# -----------------------------------------------------------------------------

chain <- import.chain(chain_file)

cpgs_data <- fread(cpg_file, header = FALSE)

setnames(
  cpgs_data,
  c("chr", "start", "end", "probe_id")
)

cpgs_data[, `:=`(
  start = as.integer(start),
  end = as.integer(end)
)]

# -----------------------------------------------------------------------------
# 3. Convert CpGs to GRanges
# -----------------------------------------------------------------------------

cpgs_gr <- GRanges(
  seqnames = cpgs_data$chr,
  ranges = IRanges(
    start = cpgs_data$start,
    end = cpgs_data$end
  ),
  probe_id = cpgs_data$probe_id
)

# -----------------------------------------------------------------------------
# 4. Lift over hg19 CpGs to mm10
# -----------------------------------------------------------------------------

cpgs_mm10 <- liftOver(cpgs_gr, chain)

# Preserve the original unlisting approach
cpgs_mm10 <- unlist(cpgs_mm10)

dt_cpgs_mm10 <- as.data.table(
  as.data.frame(cpgs_mm10)
)

# -----------------------------------------------------------------------------
# 5. Aggregate lifted-over coordinates by CpG
# -----------------------------------------------------------------------------

# Original analysis used the minimum start and maximum end
# for each CpG probe.

lift_over_res <- dt_cpgs_mm10[
  ,
  .(
    chr = unique(seqnames),
    start = min(start),
    end = max(end)
  ),
  by = probe_id
]

setcolorder(
  lift_over_res,
  c("probe_id", "chr", "start", "end")
)

# -----------------------------------------------------------------------------
# 6. Save results
# -----------------------------------------------------------------------------

fwrite(
  lift_over_res,
  output_file,
  sep = "\t",
  quote = FALSE
)


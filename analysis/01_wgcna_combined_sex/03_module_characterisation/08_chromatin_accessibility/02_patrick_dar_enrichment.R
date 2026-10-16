
#!/usr/bin/env Rscript

# =============================================================================
# Chromatin accessibility enrichment: Patrick et al. (2024)
#
# Test enrichment of age-associated co-methylation CpGs within
# differentially accessible regions (DARs) reported by Patrick et al.
#
# Reference:
# https://www.cell.com/cell-metabolism/fulltext/S1550-4131(24)00231-6
#
# Input:
#   - hg19-to-mm10 lifted-over CpG coordinates
#   - Selected co-methylation CpG sets
#   - Patrick et al. supplementary DAR data (mmc2.xlsx)
#
# Output:
#   02a_DAR_chrom_access_patrick_et_al.tsv
#
# Note:
#   The original enrichment calculations are preserved for reproducibility.
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
  library(GenomicRanges)
})

# -----------------------------------------------------------------------------
# 1. Paths
# -----------------------------------------------------------------------------

input_dir <- "/results/01_wgcna_combined_sex"
chr_access_dir <- "/data/chrom_accessibility"
opi_dir <- "/data/gsDNAm"
output_dir <- "/results/01_wgcna_combined_sex"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

opi_file <- file.path(
  input_dir,
  "GS_sub_qc1_gene_cpgs_assoc_rGREAT_h19_all_CpGs.opi"
)

liftover_file <- file.path(
  opi_dir,
  "GS_cpgs_h19_mapped_to_mm10.bed"
)

cpg_list_file <- file.path(
  input_dir,
  "heteroskid_wgcna_cpg_list.rds"
)

patrick_file <- file.path(
  chr_access_dir,
  "mmc2.xlsx"
)

output_file <- file.path(
  output_dir,
  "02a_DAR_chrom_access_patrick_et_al.tsv"
)

# -----------------------------------------------------------------------------
# 2. Load background CpGs and study CpG sets
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

# Lifted-over CpG coordinates (hg19 to mm10)
bg_cpgs_mm10 <- fread(liftover_file)

# Preserve the original column selection and ordering
bg_cpgs_mm10_bed <- bg_cpgs_mm10[, c(2, 3, 4, 1)]

# -----------------------------------------------------------------------------
# 3. DAR enrichment function
# -----------------------------------------------------------------------------

perform_atac_enrichment <- function(
    study_hit_list,
    dar_group,
    topN_dar = 10000
) {

  if (!dar_group %in% c("Opening", "Closing")) {
    stop("dar_group must be 'Opening' or 'Closing'.")
  }

  res_list <- list()
  sheet_names <- excel_sheets(patrick_file)

  for (i in 2:23) {

    sheet_name <- sheet_names[i]

    message("Processing tissue: ", sheet_name, " | ", dar_group)

    dt_dar <- as.data.table(
      read_xlsx(patrick_file, sheet = i)
    )

    # Extract genomic coordinates from Peak_coordinate
    dt_dar[, c("chr", "start", "end", "strand") :=
      tstrsplit(Peak_coordinate, "[:-]")]

    dt_dar[, `:=`(
      start = as.integer(start),
      end = as.integer(end)
    )]

    # Preserve original age-opening and age-closing definitions
    if (dar_group == "Opening") {
      dt_dar_sub <- dt_dar[log2FoldChange < 0]
    } else {
      dt_dar_sub <- dt_dar[log2FoldChange >= 0]
    }

    # Select top DARs ranked by adjusted P-value
    dt_dar_top <- dt_dar_sub[order(padj)][1:topN_dar]

    dt_dar_top_bed <- dt_dar_top[
      !is.na(chr),
      .(chr, start, end, Peak_name, log2FoldChange, padj)
    ]

    dt_dar_top_bed_gr <- GRanges(
      seqnames = dt_dar_top_bed$chr,
      ranges = IRanges(
        start = dt_dar_top_bed$start,
        end = dt_dar_top_bed$end
      ),
      peak_name = dt_dar_top_bed$Peak_name
    )

    sheet_res_list <- list()

    for (geneSet in names(study_hit_list)) {

      cpg_hit_names <- study_hit_list[[geneSet]]

      bg_cpgs_mm10_bed_sub_hits <- bg_cpgs_mm10_bed[
        probe_id %in% cpg_hit_names
      ]

      bg_cpgs_mm10_bed_sub_gr <- GRanges(
        seqnames = bg_cpgs_mm10_bed_sub_hits$chr,
        ranges = IRanges(
          start = bg_cpgs_mm10_bed_sub_hits$start,
          end = bg_cpgs_mm10_bed_sub_hits$end
        ),
        probe_id = bg_cpgs_mm10_bed_sub_hits$probe_id
      )

      # Identify CpGs overlapping Patrick et al. DARs
      overlaps <- findOverlaps(
        bg_cpgs_mm10_bed_sub_gr,
        dt_dar_top_bed_gr
      )

      overlapping_cpgs <- unique(data.table(
        probe_id = bg_cpgs_mm10_bed_sub_gr$probe_id[
          queryHits(overlaps)
        ]
      ))

      # Preserve original enrichment definitions
      N <- cpgSize
      K <- nrow(dt_dar_top_bed)
      n <- length(cpg_hit_names)
      k <- nrow(overlapping_cpgs)

      # Odds ratio and Wald 95% confidence interval
      if (
        k == 0 ||
        K == k ||
        n == k ||
        (N - K - (n - k)) == 0
      ) {

        OR <- NA_real_
        CI_lower <- NA_real_
        CI_upper <- NA_real_

      } else {

        OR <- (k / (K - k)) /
          ((n - k) / (N - K - (n - k)))

        SE_log_OR <- sqrt(
          1 / k +
          1 / (n - k) +
          1 / (K - k) +
          1 / (N - K - (n - k))
        )

        log_OR <- log(OR)

        CI_lower <- exp(log_OR - 1.96 * SE_log_OR)
        CI_upper <- exp(log_OR + 1.96 * SE_log_OR)
      }

      sheet_res_list[[geneSet]] <- data.table(
        analysis = geneSet,
        tissue = sheet_name,
        DAR_group = dar_group,
        anal_set = "DARs",
        sigDAR_peaks = K,
        cpgs_hits = n,
        overlap = k,
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

    res_list[[i]] <- rbindlist(
      sheet_res_list,
      use.names = TRUE
    )
  }

  rbindlist(res_list, use.names = TRUE)
}

# -----------------------------------------------------------------------------
# 4. Run enrichment for closing and opening DARs
# -----------------------------------------------------------------------------

res_dar_enrichment_closing <- perform_atac_enrichment(
  study_hit_list = heteroskid_wgcna_cpgs_list,
  dar_group = "Closing",
  topN_dar = 10000
)

res_dar_enrichment_opening <- perform_atac_enrichment(
  study_hit_list = heteroskid_wgcna_cpgs_list,
  dar_group = "Opening",
  topN_dar = 10000
)

# -----------------------------------------------------------------------------
# 5. Combine and save results
# -----------------------------------------------------------------------------

res_dar_enrichment_c <- rbindlist(
  list(
    res_dar_enrichment_opening,
    res_dar_enrichment_closing
  ),
  use.names = TRUE
)

fwrite(
  res_dar_enrichment_c,
  output_file,
  sep = "\t",
  quote = FALSE
)

message("DAR enrichment completed: ", output_file)


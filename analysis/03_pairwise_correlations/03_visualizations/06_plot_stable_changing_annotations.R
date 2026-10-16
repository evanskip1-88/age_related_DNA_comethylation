#!/usr/bin/env Rscript

# =============================================================================
# Supplementary Figure S4c
# Genomic annotations of stable and changing CpG pairs
#
# Inputs:
#   annotatr CpG annotation tables
#   annotatr genic annotation tables
#   GREAT nearest-gene annotation tables
#
# Output:
#   Fig_S4c_stable_changing_annotations.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyr)
  library(wesanderson)
  library(ComplexHeatmap)
  library(grid)
})

# -----------------------------------------------------------------------------
# 1. Configure paths
# -----------------------------------------------------------------------------

annot_dir <- "/results/03_pairwise_correlations"
great_dir <- "/results/03_pairwise_correlations"
output_dir <- "/figures"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 2. Input files
# -----------------------------------------------------------------------------

cpg_files <- c(
  Stable = "chr_all_AgeBaseAll_stablePairsPosNeg_heteroskid_top1000_cpgs.annotr",
  Changing = "chr_all_AgeBaseAll_changingPairsPosNeg_heteroskid_top1000_cpgs.annotr"
)

genic_files <- c(
  Stable = "chr_all_AgeBaseAll_stablePairsPosNeg_heteroskid_top1000_genic.annotr",
  Changing = "chr_all_AgeBaseAll_changingPairsPosNeg_heteroskid_top1000_genic.annotr"
)

great_files <- c(
  Stable = "chr_all_AgeBaseAll_stablePairsPosNeg_heteroskid_top1000_cpgsGREATout.txt",
  Changing = "chr_all_AgeBaseAll_changingPairsPos_heteroskid_top1000_cpgsGREATout.txt"
)

group_labels <- c("Stable", "Changing")

# -----------------------------------------------------------------------------
# 3. Helper: summarise annotatr categories
# -----------------------------------------------------------------------------

summarise_annotation <- function(files, directory, prefix, categories,
                                 deduplicate = FALSE) {

  results <- lapply(names(files), function(group) {

    dat <- fread(file.path(directory, files[[group]]))

    required <- c("seqnames", "start", "end", "annot.type")

    if (!all(required %in% names(dat))) {
      stop("Missing annotation columns for group: ", group)
    }

    dat <- dat[, ..required]

    if (deduplicate) {
      dat <- unique(dat)
    }

    dat[, category := sub(paste0("^", prefix), "", annot.type)]

    counts <- dat[
      category %in% categories,
      .N,
      by = category
    ]

    # Denominator includes all annotation records, matching
    # the original frequency-table calculation.
    total <- nrow(dat)

    values <- setNames(
      rep(0, length(categories)),
      categories
    )

    values[counts$category] <- counts$N / total

    values
  })

  mat <- do.call(rbind, results)

  rownames(mat) <- names(files)
  colnames(mat) <- categories

  mat[group_labels, , drop = FALSE]
}

# -----------------------------------------------------------------------------
# 4. CGI annotation
# -----------------------------------------------------------------------------

cpg_categories <- c(
  "cpg_islands",
  "cpg_shores",
  "cpg_shelves",
  "cpg_inter"
)

cpg_annot <- summarise_annotation(
  files = cpg_files,
  directory = annot_dir,
  prefix = "hg19_",
  categories = cpg_categories,
  deduplicate = FALSE
)

# -----------------------------------------------------------------------------
# 5. Genic annotation
# -----------------------------------------------------------------------------

genic_categories <- c(
  "1to5kb",
  "promoters",
  "5UTRs",
  "exons",
  "introns",
  "intergenic",
  "3UTRs"
)

genic_annot <- summarise_annotation(
  files = genic_files,
  directory = annot_dir,
  prefix = "hg19_genes_",
  categories = genic_categories,
  deduplicate = TRUE
)

# -----------------------------------------------------------------------------
# 6. GREAT distance-to-TSS annotation
# -----------------------------------------------------------------------------

tss_categories <- c(
  "<1kb",
  "1kb~5kb",
  "5kb~10kb",
  ">10kb"
)

summarise_tss <- function(filename) {

  dat <- fread(
    file.path(great_dir, filename),
    skip = 1,
    header = FALSE
  )

  if (!"V2" %in% names(dat)) {
    stop("Expected V2 column in GREAT output.")
  }

  # Original GREAT field contains gene and signed distance
  # e.g. GENE (-2500)

  distance_text <- sub(
    ".*\\(([-+]?[0-9]+)\\).*",
    "\\1",
    dat$V2
  )

  distance_bp <- suppressWarnings(as.numeric(distance_text))

  abs_distance <- abs(distance_bp)

  categories <- fifelse(
    abs_distance < 1000,
    "<1kb",
    fifelse(
      abs_distance > 1000 & abs_distance <= 5000,
      "1kb~5kb",
      fifelse(
        abs_distance > 5000 & abs_distance <= 10000,
        "5kb~10kb",
        fifelse(abs_distance > 10000, ">10kb", NA_character_)
      )
    )
  )

  counts <- table(
    factor(categories, levels = tss_categories)
  )

  if (sum(counts) == 0L) {
    stop("No valid GREAT distance categories found: ", filename)
  }

  as.numeric(counts) / sum(counts)
}

dist_TSS <- rbind(
  Stable = summarise_tss(great_files[["Stable"]]),
  Changing = summarise_tss(great_files[["Changing"]])
)

colnames(dist_TSS) <- tss_categories

# -----------------------------------------------------------------------------
# 7. Colours
# -----------------------------------------------------------------------------

cpgs_col <- c(
  "#008000",
  "#C2E105",
  "#8A91D0",
  "#CD5C5C"
)

dist_tss_col <- c(
  "#FF0000",
  "#FF7352",
  "#FFB299",
  "#FFD9CB"
)

genic_col <- wes_palette(
  "Zissou1",
  10,
  type = "continuous"
)[seq_len(ncol(genic_annot))]

# -----------------------------------------------------------------------------
# 8. Construct annotations
# -----------------------------------------------------------------------------

anno_width <- unit(4, "cm")

ht_label <- rowAnnotation(
  text = anno_text(
    group_labels,
    location = unit(1, "npc"),
    just = "right",
    gp = gpar(fontsize = 14)
  )
)

ht_cpgs <- rowAnnotation(
  cpgs_annot = anno_barplot(
    cpg_annot,
    bar_width = 1,
    gp = gpar(fill = cpgs_col),
    width = anno_width
  ),
  show_annotation_name = FALSE
)

ht_tss <- rowAnnotation(
  dist_tss = anno_barplot(
    dist_TSS,
    bar_width = 1,
    gp = gpar(fill = dist_tss_col),
    width = anno_width
  ),
  show_annotation_name = FALSE
)

ht_genic <- rowAnnotation(
  genic_annot = anno_barplot(
    genic_annot,
    bar_width = 1,
    gp = gpar(fill = genic_col),
    width = anno_width
  ),
  show_annotation_name = FALSE
)

ht_list <- ht_label + ht_cpgs + ht_tss + ht_genic

# -----------------------------------------------------------------------------
# 9. Legends
# -----------------------------------------------------------------------------

cgi_lgd <- Legend(
  labels = cpg_categories,
  title = "CGI annotation",
  legend_gp = gpar(fill = cpgs_col)
)

tss_lgd <- Legend(
  labels = tss_categories,
  title = "Distance to TSS",
  legend_gp = gpar(fill = dist_tss_col)
)

genic_lgd <- Legend(
  labels = genic_categories,
  title = "Gene annotation",
  legend_gp = gpar(fill = genic_col)
)

packed_legend <- packLegend(
  cgi_lgd,
  tss_lgd,
  genic_lgd,
  direction = "horizontal"
)

# -----------------------------------------------------------------------------
# 10. Export Supplementary Figure S4c
# -----------------------------------------------------------------------------

output_file <- file.path(
  output_dir,
  "Fig_S4c_stable_changing_annotations.pdf"
)

pdf(
  output_file,
  width = 14,
  height = 5
)

draw(
  ht_list,
  padding = unit(c(40, 5, 20, 2), "mm")
)

draw(
  packed_legend,
  x = unit(0, "npc"),
  y = unit(0, "npc"),
  just = c("left", "bottom")
)

anno_titles <- c(
  cpgs_annot = "CGI\nannotation",
  dist_tss = "Distance\nto TSS",
  genic_annot = "Gene\nannotation"
)

for (annotation_name in names(anno_titles)) {

  decorate_annotation(annotation_name, {
    grid.text(
      anno_titles[[annotation_name]],
      y = unit(1, "npc") + unit(3, "mm"),
      just = "bottom"
    )
  })
}

dev.off()

message("Supplementary Figure S4c saved: ", output_file)

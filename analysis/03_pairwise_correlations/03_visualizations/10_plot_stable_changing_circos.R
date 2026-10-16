#!/usr/bin/env Rscript

# =============================================================================
# Supplementary Figure S4c
# Circos representation of genomic links between CpG pairs
# classified as stable or changing by heteroskedasticity analysis.
#
# Inputs:
#   1. Top stable/changing CpG pairs
#   2. CpG genomic coordinates (hg19)
#   3. Chromosome coordinate ranges
#
# Output:
#   Fig_S4c_stable_changing_circos.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(circlize)
  library(viridis)
})

# -----------------------------------------------------------------------------
# 1. Input/output paths
# -----------------------------------------------------------------------------

pairs_file <- "..."
probe_annotation_file <- "..."
gene_location_file <- "..."
output_dir <- "..."

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 2. Load genomic annotation
# -----------------------------------------------------------------------------

opi <- fread(probe_annotation_file)

required_opi <- c("probe_id", "chr", "position")

if (!all(required_opi %in% names(opi))) {
  stop("CpG annotation must contain: probe_id, chr, position")
}

opi <- unique(
  opi[, .(probe_id, chr, position)],
  by = "probe_id"
)

opi[, chr := sub("^chr", "", as.character(chr))]
opi[, position := as.numeric(position)]

# -----------------------------------------------------------------------------
# 3. Load chromosome coordinate information
# -----------------------------------------------------------------------------

geneLoc <- fread(gene_location_file)

# Original file structure:
# V2 = chromosome
# V4 = genomic coordinate range

setnames(geneLoc, "V2", "Chr")

geneLoc[, Chr := as.character(Chr)]

geneLoc <- geneLoc[
  Chr %in% as.character(1:22)
]

geneLoc[, Chr := factor(Chr, levels = as.character(1:22))]

setorder(geneLoc, Chr)

# -----------------------------------------------------------------------------
# 4. Load stable/changing CpG pairs
# -----------------------------------------------------------------------------

dt_pairs <- fread(pairs_file)

required_pairs <- c(
  "probe1",
  "probe2",
  "sigGroup"
)

if (!all(required_pairs %in% names(dt_pairs))) {
  stop("Pair file must contain: probe1, probe2, sigGroup")
}

# -----------------------------------------------------------------------------
# 5. Prepare CpG coordinates
# -----------------------------------------------------------------------------

prepare_links <- function(dt, annotation) {

  dt <- copy(dt)

  ann1 <- copy(annotation)
  setnames(
    ann1,
    c("probe_id", "chr", "position"),
    c("probe1", "chr1", "pos1")
  )

  ann2 <- copy(annotation)
  setnames(
    ann2,
    c("probe_id", "chr", "position"),
    c("probe2", "chr2", "pos2")
  )

  dt <- merge(
    dt,
    ann1,
    by = "probe1",
    all.x = TRUE,
    sort = FALSE
  )

  dt <- merge(
    dt,
    ann2,
    by = "probe2",
    all.x = TRUE,
    sort = FALSE
  )

  dt <- dt[
    !is.na(chr1) &
    !is.na(chr2) &
    !is.na(pos1) &
    !is.na(pos2)
  ]

  # Retain autosomal pairs
  dt <- dt[
    chr1 %in% as.character(1:22) &
    chr2 %in% as.character(1:22)
  ]

  dt[, height := 0.5]

  dt
}

# -----------------------------------------------------------------------------
# 6. Circos plotting function
# -----------------------------------------------------------------------------

plot_circos_pairs <- function(dt_links, chromosome_info, plot_title) {

  circos.clear()

  circos.par(
    track.height = 0.05,
    cell.padding = c(0.04, 0.04, 0.04, 0.04),
    gap.degree = 1
  )

  circos.initialize(
    factors = chromosome_info$Chr,
    x = chromosome_info$V4
  )

  chr_colors <- viridis::viridis(
    22,
    option = "magma"
  )

  circos.track(
    ylim = c(0, 1),
    bg.col = chr_colors,
    bg.border = NA,
    track.height = 0.1,
    panel.fun = function(x, y) {

      circos.text(
        CELL_META$xcenter,
        CELL_META$cell.ylim[2] + mm_y(3),
        CELL_META$sector.index,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0.4, 0.5),
        font = 2.5
      )

      circos.axis(
        h = "top",
        major.tick = FALSE,
        labels.cex = 0.02
      )
    }
  )

  # Prepare genomic links
  bed1 <- dt_links[
    ,
    .(
      chr = chr1,
      start = pos1,
      end = pos1,
      probe_id = probe1,
      height
    )
  ]

  bed2 <- dt_links[
    ,
    .(
      chr = chr2,
      start = pos2,
      end = pos2,
      probe_id = probe2,
      height
    )
  ]

  bed1[, chr := factor(chr, levels = as.character(1:22))]
  bed2[, chr := factor(chr, levels = as.character(1:22))]

  if (nrow(bed1) > 0) {
    circos.genomicLink(
      bed1,
      bed2
    )
  }

  title(
    main = plot_title,
    line = -1
  )

  circos.clear()
}

# -----------------------------------------------------------------------------
# 7. Generate two-panel figure
# -----------------------------------------------------------------------------

output_file <- file.path(
  output_dir,
  "Fig_S4c_stable_changing_circos.pdf"
)

pdf(
  output_file,
  width = 9,
  height = 14
)

par(
  mfrow = c(2, 1),
  mar = c(1.5, 1.5, 1.5, 2),
  xpd = NA
)

for (sig_set in c("Stable", "Changing")) {

  dt_sub <- dt_pairs[
    sigGroup == sig_set
  ]

  dt_links <- prepare_links(
    dt_sub,
    opi
  )

  message(
    sig_set,
    ": ",
    nrow(dt_links),
    " genomic links"
  )

  plot_circos_pairs(
    dt_links,
    geneLoc,
    sig_set
  )
}

dev.off()


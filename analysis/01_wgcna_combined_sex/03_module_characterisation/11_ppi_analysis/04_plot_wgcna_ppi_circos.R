
# =============================================================================
# Purpose: Plot genomic Circos networks of WGCNA-associated PPI interactions.
# Inputs:  WGCNA module summary, rGREAT annotations, GRCh37 gene locations.
# Outputs: Supplementary Figures 11a and 11b.
# =============================================================================

library(biomaRt)
library(rbioapi)
library(circlize)
library(data.table)
library(readxl)
library(viridis)

# -------------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------------

module_file <- "/results/01_wgcna_combined_sex/01a_selected_mods_summary_stats.xlsx"

gene_location_file <- "data/gsDNAm/NCBI37.3.gene.loc"

output_dir <- "/figures"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Load module data
# -------------------------------------------------------------------------

all_mods_stats <- as.data.table(read_xlsx(module_file))

all_mods_stats[, chr := sub("^chr", "", as.character(chr))]
all_mods_stats[, moduleNumID := as.character(moduleNumID)]

geneLoc <- fread(gene_location_file, header = FALSE)

setnames(geneLoc, "V2", "Chr")

geneLoc[, Chr := as.character(Chr)]
geneLoc <- geneLoc[order(as.numeric(Chr))]

# -------------------------------------------------------------------------
# Ensembl gene coordinates
# -------------------------------------------------------------------------

# GRCh37 host is used to match the original NCBI37.3 reference.
# Verify this against the Ensembl release used in the original analysis.

mart <- useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = "https://grch37.ensembl.org"
)

# -------------------------------------------------------------------------
# Plot settings
# -------------------------------------------------------------------------

module_colours <- list(
  Changing = c(
    "purple", "red", "pink", "tan", "orange", "blue",
    "green", "grey70", "brown4", "deeppink2", "cyan"
  ),
  Stable = c(
    "purple", "red", "pink", "tan", "orange", "blue",
    "green", "brown4", "grey90", "deeppink2",
    "cyan", "darkgreen", "goldenrod", "turquoise4",
    "orchid", "dodgerblue", "firebrick", "gold",
    "mediumpurple"
  )
)

figure_files <- c(
  Stable = "Figure_S11a_WGCNA_PPI_Circos.pdf",
  Changing = "Figure_S11b_WGCNA_PPI_Circos.pdf"
)

# -------------------------------------------------------------------------
# Generate Circos plots
# -------------------------------------------------------------------------

for (sigGrp in c("Stable", "Changing")) {

  cat("Processing:", sigGrp, "\n")

  dt_sub <- copy(all_mods_stats[Category == sigGrp])

  module_lst <- as.character(unique(dt_sub$moduleNumID))

  if (length(module_lst) == 0L) {
    warning("No modules found for category: ", sigGrp)
    next
  }

  # Original module colour assignment
  sel_col <- module_colours[[sigGrp]]

  if (length(module_lst) > length(sel_col)) {
    stop("Insufficient colours for ", sigGrp, " modules.")
  }

  mod_col <- setNames(
    sel_col[seq_along(module_lst)],
    module_lst
  )

  # -----------------------------------------------------------------------
  # Retrieve STRING interactions
  # -----------------------------------------------------------------------

  dt_interact_list <- list()

  for (mod in module_lst) {

    cat("Retrieving STRING interactions for module:", mod, "\n")

    sub <- dt_sub[moduleNumID == mod]

    gene_list <- unique(na.omit(
      unlist(strsplit(sub$gene_rgreat, ";", fixed = TRUE))
    ))

    gene_list <- trimws(gene_list)
    gene_list <- gene_list[nzchar(gene_list)]

    if (length(gene_list) < 2L) {
      next
    }

    gene_interactions <- rba_string_interactions_network(
      gene_list,
      species = 9606,
      required_score = 0
    )

    if (!is.null(gene_interactions) && nrow(gene_interactions) > 0L) {

      g_int_sub <- as.data.table(
        gene_interactions[
          , c("preferredName_A", "preferredName_B", "score")
        ]
      )

    } else {

      g_int_sub <- data.table(
        preferredName_A = gene_list,
        preferredName_B = gene_list,
        score = NA_real_
      )
    }

    g_int_sub[, moduleNumID := mod]

    dt_interact_list[[mod]] <- g_int_sub
  }

  if (length(dt_interact_list) == 0L) {
    warning("No STRING results for ", sigGrp)
    next
  }

  dt_interactions <- rbindlist(
    dt_interact_list,
    use.names = TRUE
  )

  # -----------------------------------------------------------------------
  # Retrieve genomic coordinates
  # -----------------------------------------------------------------------

  gene_symbols <- unique(c(
    dt_interactions$preferredName_A,
    dt_interactions$preferredName_B
  ))

  gene_symbols <- gene_symbols[!is.na(gene_symbols)]

  gene_coords <- as.data.table(
    getBM(
      attributes = c(
        "hgnc_symbol",
        "chromosome_name",
        "start_position",
        "end_position"
      ),
      filters = "hgnc_symbol",
      values = gene_symbols,
      mart = mart
    )
  )

  # Exclude alternative assembly loci, retaining autosomes.
  gene_coords_sub <- gene_coords[
    chromosome_name %in% as.character(1:22)
  ]

  # -----------------------------------------------------------------------
  # Join PPI interactions with gene coordinates
  # -----------------------------------------------------------------------

  network_with_coords <- merge(
    dt_interactions[
      , .(preferredName_A, preferredName_B, moduleNumID)
    ],
    gene_coords_sub,
    by.x = "preferredName_A",
    by.y = "hgnc_symbol",
    all.x = TRUE
  )

  network_with_coords <- merge(
    gene_coords_sub,
    network_with_coords,
    by.x = "hgnc_symbol",
    by.y = "preferredName_B",
    all.y = TRUE,
    suffixes = c("_from", "_to")
  )

  network_with_coords <- na.omit(
    as.data.table(network_with_coords)
  )

  # -----------------------------------------------------------------------
  # Prepare chromosome tracks
  # -----------------------------------------------------------------------

  chr_levels <- as.character(1:22)

  geneLoc_sub <- copy(
    geneLoc[Chr %in% all_mods_stats$chr]
  )

  geneLoc_sub[
    , Chr := factor(Chr, levels = chr_levels)
  ]

  geneLoc_sub <- geneLoc_sub[!is.na(Chr)]

  # -----------------------------------------------------------------------
  # Save Circos plot
  # -----------------------------------------------------------------------

  pdf(
    file = file.path(output_dir, figure_files[[sigGrp]]),
    width = 11,
    height = 9
  )

  par(mar = c(1.5, 5, 1, 5), xpd = TRUE)

  circos.clear()
  circos.par(start.degree = 90)

  circos.initialize(
    factors = geneLoc_sub$Chr,
    x = geneLoc_sub$V4
  )

  circos.track(
    ylim = c(0, 1),
    bg.col = viridis(22, option = "magma"),
    bg.border = NA,
    track.height = 0.1,
    panel.fun = function(x, y) {

      circos.text(
        CELL_META$xcenter,
        CELL_META$cell.ylim[2] + mm_y(3.5),
        CELL_META$sector.index,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0.4, 0.5),
        font = 2,
        cex = 1.2
      )

      circos.axis(
        h = "top",
        major.tick = FALSE,
        labels.cex = 0.001,
        col = "black"
      )
    }
  )

  # -----------------------------------------------------------------------
  # Draw module-specific genomic PPI links
  # -----------------------------------------------------------------------

  for (mod in module_lst) {

    network_mod <- network_with_coords[
      moduleNumID == mod
    ]

    if (nrow(network_mod) == 0L) {
      next
    }

    bed1 <- network_mod[
      , .(
        chr = chromosome_name_from,
        start = start_position_from,
        end = end_position_from
      )
    ]

    bed2 <- network_mod[
      , .(
        chr = chromosome_name_to,
        start = start_position_to,
        end = end_position_to
      )
    ]

    bed1[, chr := as.character(chr)]
    bed2[, chr := as.character(chr)]

    circos.genomicLink(
      bed1,
      bed2,
      col = mod_col[[mod]]
    )
  }

  # -----------------------------------------------------------------------
  # Draw CpG positions
  # -----------------------------------------------------------------------

  for (mod in module_lst) {

    sub <- copy(dt_sub[moduleNumID == mod])

    if (nrow(sub) == 0L) {
      next
    }

    sub[, height := 0.5]

    sub <- unique(
      sub[, .(probe_id, chr, position, height)]
    )

    circos.trackLines(
      sub$chr,
      sub$position,
      sub$height,
      lwd = 2,
      type = "h"
    )
  }

  # -----------------------------------------------------------------------
  # Module legend
  # -----------------------------------------------------------------------

  sorted_ids <- sort(as.numeric(names(mod_col)))
  sorted_cols <- mod_col[as.character(sorted_ids)]

  legend(
    "right",
    legend = sorted_ids,
    col = sorted_cols,
    lty = 1,
    lwd = 2.5,
    title = "Module ID",
    inset = c(-0.20, 0)
  )

  circos.clear()
  dev.off()

  cat("Saved:", figure_files[[sigGrp]], "\n")
}



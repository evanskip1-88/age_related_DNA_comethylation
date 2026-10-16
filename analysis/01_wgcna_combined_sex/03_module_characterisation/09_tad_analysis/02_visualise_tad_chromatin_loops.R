
#!/usr/bin/env Rscript

# =============================================================================
# TAD boundary and chromatin loop enrichment heatmap
#
# Input:
#   01a_chrom_loop_enrich.tsv
#   01b_tad_boundery_enrich.tsv
#
# Output:
#   TAD_chromatin_loop_enrichment_WGCNA.pdf
#   TAD_chromatin_loop_enrichment_Heteroskid.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# -----------------------------------------------------------------------------
# 1. Paths and settings
# -----------------------------------------------------------------------------

input_dir <- "results/tad_chromatin_loop"
output_dir <- input_dir

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Choose "WGCNA" or "Heteroskid"
analName <- "WGCNA"

# -----------------------------------------------------------------------------
# 2. Load enrichment results
# -----------------------------------------------------------------------------

chrom_loop <- fread(
  file.path(input_dir, "01a_chrom_loop_enrich.tsv")
)

tad <- fread(
  file.path(input_dir, "01b_tad_boundery_enrich.tsv")
)

# Chromatin loops: retain ENCODE datasets for selected analysis
chrom_loop_sub <- copy(
  chrom_loop[
    grepl("ENCODE", tissue) &
    startsWith(analysis, paste0(analName, "_"))
  ]
)

chrom_loop_sub[, tissue_clean := gsub(
  "^ENCODE3\\.|\\.hg19\\.peakachu-merged$",
  "",
  tissue
)]

chrom_loop_sub[
  ,
  tissue_clean := sub("LNCAP", "LNCaP", tissue_clean)
]

# TAD boundaries: retain original raw datasets for selected analysis
tad_sub <- copy(
  tad[
    grepl("_raw-", tissue) &
    startsWith(analysis, paste0(analName, "_"))
  ]
)

tad_sub[, tissue_clean := sub("_raw.*", "", tissue)]
tad_sub[, tissue_clean := sub("T470", "T47D", tissue_clean)]

# -----------------------------------------------------------------------------
# 3. Combine chromatin loops and TAD boundaries
# -----------------------------------------------------------------------------

chrom_loop_sub2 <- chrom_loop_sub[
  ,
  .(analysis, tissue_clean, pvalue, OR, CI_lower, CI_upper)
]

chrom_loop_sub2[, Group := "Chromatin_Loops"]

tad_sub2 <- tad_sub[
  ,
  .(analysis, tissue_clean, pvalue, OR, CI_lower, CI_upper)
]

tad_sub2[, Group := "TAD_boundary"]

dt_genom_org_res <- rbindlist(
  list(chrom_loop_sub2, tad_sub2),
  use.names = TRUE
)

dt_genom_org_res[
  ,
  col_id := paste(Group, analysis, sep = "_")
]

# Preserve original cell-line exclusions
dt_genom_org_res <- dt_genom_org_res[
  !tissue_clean %in% c("SJCRH30", "SKNMC", "SKNDZ")
]

# -----------------------------------------------------------------------------
# 4. Transform and annotate results
# -----------------------------------------------------------------------------

dt_genom_org_res[, `:=`(
  OR = round(log(OR), 2),
  CI_lower = log(CI_lower),
  CI_upper = log(CI_upper)
)]

dt_genom_org_res[, OR_CI_label := fifelse(
  is.na(OR) | is.na(CI_lower) | is.na(CI_upper),
  NA_character_,
  sprintf("%.2f (%.2f,%.2f)", OR, CI_lower, CI_upper)
)]

dt_genom_org_res[, p_stars := fifelse(
  pvalue < 0.001, "***",
  fifelse(
    pvalue < 0.01, "**",
    fifelse(pvalue < 0.05, "*", NA_character_)
  )
)]

dt_genom_org_res[, pval_fmt := fifelse(
  is.na(pvalue),
  NA_character_,
  gsub(
    "e([+-]?[0-9]+)",
    " %*% 10^{\\1}",
    formatC(pvalue, format = "e", digits = 1)
  )
)]

# -----------------------------------------------------------------------------
# 5. Define column order
# -----------------------------------------------------------------------------

if (analName == "Heteroskid") {

  col_names_order <- c(
    "Chromatin_Loops_Heteroskid_Stable",
    "TAD_boundary_Heteroskid_Stable",
    "Chromatin_Loops_Heteroskid_Changing",
    "TAD_boundary_Heteroskid_Changing"
  )

  col_split_factor <- factor(
    rep(c("aSCPs", "aVCPs"), each = 2),
    levels = c("aSCPs", "aVCPs")
  )

} else {

  col_names_order <- c(
    "Chromatin_Loops_WGCNA_Stable",
    "TAD_boundary_WGCNA_Stable",
    "Chromatin_Loops_WGCNA_Changing",
    "TAD_boundary_WGCNA_Changing"
  )

  col_split_factor <- factor(
    rep(c("aSCMs", "aVCMs"), each = 2),
    levels = c("aSCMs", "aVCMs")
  )
}

col_names_labs <- rep(
  c("Chromatin_Loops", "TAD_Boundary"),
  2
)

# -----------------------------------------------------------------------------
# 6. Construct aligned matrices
# -----------------------------------------------------------------------------

make_matrix <- function(dt, value_col) {

  wide <- dcast(
    dt,
    tissue_clean ~ col_id,
    value.var = value_col
  )

  missing_cols <- setdiff(col_names_order, names(wide))

  for (nm in missing_cols) {
    wide[, (nm) := NA]
  }

  mat <- as.matrix(wide[, ..col_names_order])
  rownames(mat) <- wide$tissue_clean

  mat
}

or_matrix <- make_matrix(dt_genom_org_res, "OR")
storage.mode(or_matrix) <- "numeric"

ci_matrix <- make_matrix(dt_genom_org_res, "OR_CI_label")
pval_matrix <- make_matrix(dt_genom_org_res, "pval_fmt")
pstar_matrix <- make_matrix(dt_genom_org_res, "p_stars")

# -----------------------------------------------------------------------------
# 7. Heatmap colours
# -----------------------------------------------------------------------------

log_min <- min(or_matrix, na.rm = TRUE)
log_max <- max(or_matrix, na.rm = TRUE)

mat_color <- colorRamp2(
  c(log_min, 0, log_max),
  c("skyblue", "white", "red")
)

# -----------------------------------------------------------------------------
# 8. Construct heatmap
# -----------------------------------------------------------------------------

ht <- Heatmap(
  or_matrix,
  name = "OR",
  col = mat_color,
  na_col = "grey90",

  cluster_rows = FALSE,
  cluster_columns = FALSE,

  row_names_side = "left",
  row_labels = rownames(or_matrix),

  column_split = col_split_factor,
  column_names_rot = 45,
  column_labels = col_names_labs,

  column_title_gp = gpar(fontface = "italic", fontsize = 14),
  column_names_gp = gpar(fontsize = 16, fontface = "bold"),
  row_names_gp = gpar(fontsize = 14, fontface = "bold"),

  border_gp = gpar(col = "black", lty = 1, lwd = 1.5),
  rect_gp = gpar(col = "black", lwd = 0.2),

  cell_fun = function(j, i, x, y, width, height, fill) {

    ci_val <- ci_matrix[i, j]
    pval_val <- pval_matrix[i, j]
    star_val <- pstar_matrix[i, j]

    if (!is.na(ci_val)) {
      grid.text(
        ci_val,
        x,
        y + unit(4, "mm"),
        gp = gpar(fontsize = 14)
      )
    }

    if (!is.na(pval_val)) {
      grid.text(
        parse(text = pval_val),
        x,
        y - unit(1, "mm"),
        gp = gpar(fontsize = 14)
      )
    }

    if (!is.na(star_val) && nzchar(star_val)) {
      grid.text(
        star_val,
        x,
        y - unit(6, "mm"),
        gp = gpar(fontsize = 14)
      )
    }
  },

  column_gap = unit(2, "mm"),
  row_gap = unit(1, "mm"),

  row_title = "Tissue cell lines",
  row_title_gp = gpar(fontsize = 15, fontface = "bold")
)

# -----------------------------------------------------------------------------
# 9. Save heatmap
# -----------------------------------------------------------------------------

output_file <- file.path(
  output_dir,
  paste0("TAD_chromatin_loop_enrichment_", analName, ".pdf")
)

pdf(output_file, width = 14, height = 12)
draw(ht)
dev.off()

message("Heatmap saved: ", output_file)


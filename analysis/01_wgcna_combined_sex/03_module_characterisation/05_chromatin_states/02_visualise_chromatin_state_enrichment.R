#==============================================================================
# Visualise chromatin-state enrichment
#
# Heatmap of log odds ratios, 95% confidence intervals and P-values
# for age-stable and age-variable co-methylation CpG sets.
#
# Published ChromHMM chromatin states (hg19).
#==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

rm(list = ls())

#------------------------------------------------------------------------------
# Paths
#------------------------------------------------------------------------------

input_path <- "/results/01_wgcna_combined_sex"

#------------------------------------------------------------------------------
# Load enrichment results
#------------------------------------------------------------------------------

res_enrich <- fread(
  file.path(
    input_path,
    "01a_chrom_state_enrich_combined_mods_n_het.tsv"
  )
)

#------------------------------------------------------------------------------
# Select top 10 chromatin states per analysis
#------------------------------------------------------------------------------

chrom_res <- res_enrich[
  overlap > 10 & !is.na(OR) & OR > 0 & !is.na(chrom_state)
]

chrom_res <- chrom_res[
  order(-OR),
  head(.SD, 10),
  by = analysis
]

#------------------------------------------------------------------------------
# Group chromatin states into functional categories
#------------------------------------------------------------------------------

chrom_res[
  ,
  chrom_group := fcase(
    grepl("Enh", chrom_state),     "Enhancers",
    grepl("BivProm", chrom_state), "Bivalent promoters",
    grepl("TSS", chrom_state),     "TSS",
    grepl("PromF", chrom_state),   "Core promoters",
    grepl("ReprPC", chrom_state),  "ReprPC",
    grepl("HET", chrom_state),     "HET",
    default = "Other"
  )
]

group_levels <- c(
  "TSS",
  "Core promoters",
  "Bivalent promoters",
  "Enhancers",
  "ReprPC",
  "HET",
  "Other"
)

chrom_res[
  ,
  chrom_group := factor(
    chrom_group,
    levels = group_levels
  )
]

#------------------------------------------------------------------------------
# Transform odds ratios and confidence intervals
#------------------------------------------------------------------------------

chrom_res[
  ,
  `:=`(
    log_OR = log(OR),
    log_CI_lower = log(CI_lower),
    log_CI_upper = log(CI_upper)
  )
]

chrom_res[
  ,
  OR_CI_label := fifelse(
    is.na(log_OR) | is.na(log_CI_lower) | is.na(log_CI_upper),
    NA_character_,
    sprintf(
      "%.2f (%.2f–%.2f)",
      log_OR,
      log_CI_lower,
      log_CI_upper
    )
  )
]

chrom_res[
  ,
  pval_label := fifelse(
    is.na(pvalue),
    NA_character_,
    formatC(pvalue, format = "e", digits = 1)
  )
]

#------------------------------------------------------------------------------
# Define row and column order
#------------------------------------------------------------------------------

analysis_order <- c(
  "WGCNA_Stable",
  "Heteroskid_Stable",
  "WGCNA_Changing",
  "Heteroskid_Changing"
)

analysis_labels <- c(
  "aSCMs",
  "aSCPs",
  "aVCMs",
  "aVCPs"
)

row_info <- unique(
  chrom_res[, .(chrom_state, chrom_group)]
)

setorder(row_info, chrom_group, chrom_state)

state_order <- row_info$chrom_state

#------------------------------------------------------------------------------
# Construct matrices
#------------------------------------------------------------------------------

make_matrix <- function(value_column) {

  mat <- dcast(
    chrom_res,
    chrom_state ~ analysis,
    value.var = value_column
  )

  mat <- as.data.frame(mat)
  rownames(mat) <- mat$chrom_state
  mat$chrom_state <- NULL

  # Preserve all four analysis groups
  missing_cols <- setdiff(analysis_order, colnames(mat))

  for (col in missing_cols) {
    mat[[col]] <- NA
  }

  as.matrix(
    mat[state_order, analysis_order, drop = FALSE]
  )
}

or_matrix <- make_matrix("log_OR")
ci_matrix <- make_matrix("OR_CI_label")
pval_matrix <- make_matrix("pval_label")

storage.mode(or_matrix) <- "numeric"

#------------------------------------------------------------------------------
# Heatmap annotations
#------------------------------------------------------------------------------

row_groups <- factor(
  row_info$chrom_group,
  levels = group_levels
)

row_labels <- sub(
  "^[0-9]+_",
  "",
  rownames(or_matrix)
)

column_groups <- factor(
  c("Stable", "Stable", "Changing", "Changing"),
  levels = c("Stable", "Changing")
)

#------------------------------------------------------------------------------
# Heatmap colours
#------------------------------------------------------------------------------

log_min <- min(-1, min(or_matrix, na.rm = TRUE))
log_max <- max(1, max(or_matrix, na.rm = TRUE))

mat_col <- colorRamp2(
  c(log_min, 0, log_max),
  c("skyblue", "white", "red")
)

#------------------------------------------------------------------------------
# Plot chromatin-state enrichment heatmap
#------------------------------------------------------------------------------

ht <- Heatmap(
  or_matrix,
  name = "log(OR)",
  col = mat_col,
  na_col = "grey90",

  cluster_rows = FALSE,
  cluster_columns = FALSE,

  row_labels = row_labels,
  row_names_side = "left",

  column_labels = analysis_labels,
  column_split = column_groups,
  column_names_rot = 45,

  row_split = row_groups,
  row_title = NULL,
  column_title = NULL,

  border_gp = gpar(
    col = "black",
    lwd = 1.5
  ),

  rect_gp = gpar(
    col = "black",
    lwd = 0.2
  ),

  row_names_gp = gpar(
    fontsize = 14,
    fontface = "bold"
  ),

  column_names_gp = gpar(
    fontsize = 14,
    fontface = "bold"
  ),

  cell_fun = function(j, i, x, y, width, height, fill) {

    ci_val <- ci_matrix[i, j]
    pval_val <- pval_matrix[i, j]

    if (!is.na(ci_val)) {
      grid.text(
        ci_val,
        x,
        y + unit(2, "mm"),
        gp = gpar(fontsize = 11)
      )
    }

    if (!is.na(pval_val)) {
      grid.text(
        pval_val,
        x,
        y - unit(2.5, "mm"),
        gp = gpar(fontsize = 10)
      )
    }
  },

  column_gap = unit(2, "mm"),
  row_gap = unit(1, "mm")
)

draw(ht)

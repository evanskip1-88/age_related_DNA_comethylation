
#!/usr/bin/env Rscript

# =============================================================================
# Chromatin accessibility enrichment heatmap
# Patrick et al. (2024) age-associated DARs
#
# Input:
#   02a_DAR_chrom_access_patrick_et_al.tsv
#
# Output:
#   Patrick_DAR_enrichment_heatmap.pdf
#
# Values displayed are natural-log-transformed odds ratios and 95% CIs.
# P-values are from the original enrichment analysis.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# -----------------------------------------------------------------------------
# 1. Paths
# -----------------------------------------------------------------------------

input_file <- "results/01_wgcna_combined_sex/02a_DAR_chrom_access_patrick_et_al.tsv"
output_dir <- "/figures"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_file <- file.path(
  output_dir,
  "Patrick_DAR_enrichment_heatmap.pdf"
)

# -----------------------------------------------------------------------------
# 2. Load and prepare enrichment results
# -----------------------------------------------------------------------------

res_dars <- fread(input_file)

# Preserve the original overlap threshold
res_dars_sub <- copy(res_dars[overlap > 5])

res_dars_sub[, tissue := sub(" DARs", "", tissue)]

# Natural-log transform ORs and confidence intervals
res_dars_sub[, `:=`(
  OR = round(log(OR), 2),
  CI_lower = log(CI_lower),
  CI_upper = log(CI_upper)
)]

res_dars_sub[, col_id := paste(DAR_group, analysis, sep = "_")]

# -----------------------------------------------------------------------------
# 3. Define published column order and labels
# -----------------------------------------------------------------------------

col_mat_order <- c(
  "Opening_Heteroskid_Stable",
  "Closing_Heteroskid_Stable",
  "Opening_Heteroskid_Changing",
  "Closing_Heteroskid_Changing",
  "Opening_WGCNA_Stable",
  "Closing_WGCNA_Stable",
  "Opening_WGCNA_Changing",
  "Closing_WGCNA_Changing"
)

col_names_labs <- rep(
  c("DAR_opening", "DAR_closing"),
  4
)

col_split_group <- factor(
  rep(c("aSCPs", "aVCPs", "aSCMs", "aVCMs"), each = 2),
  levels = c("aSCPs", "aVCPs", "aSCMs", "aVCMs")
)

# -----------------------------------------------------------------------------
# 4. Construct aligned heatmap matrices
# -----------------------------------------------------------------------------

make_matrix <- function(dt, value_col) {

  wide <- dcast(
    dt,
    tissue ~ col_id,
    value.var = value_col
  )

  missing_cols <- setdiff(col_mat_order, names(wide))

  for (nm in missing_cols) {
    wide[, (nm) := NA]
  }

  mat <- as.matrix(wide[, ..col_mat_order])
  rownames(mat) <- wide$tissue

  mat
}

# Log odds ratios
or_matrix <- make_matrix(res_dars_sub, "OR")
storage.mode(or_matrix) <- "numeric"

# Log OR and 95% CI labels
res_dars_sub[, OR_CI_label := fifelse(
  is.na(OR) | is.na(CI_lower) | is.na(CI_upper),
  NA_character_,
  sprintf(
    "%.2f (%.2f-%.2f)",
    OR, CI_lower, CI_upper
  )
)]

ci_matrix <- make_matrix(res_dars_sub, "OR_CI_label")

# Significance stars
res_dars_sub[, p_stars := fifelse(
  pvalue < 0.001, "***",
  fifelse(
    pvalue < 0.01, "**",
    fifelse(pvalue < 0.05, "*", NA_character_)
  )
)]

pstar_matrix <- make_matrix(res_dars_sub, "p_stars")

# Scientific-notation P-value labels
res_dars_sub[, pval_fmt := fifelse(
  is.na(pvalue),
  NA_character_,
  gsub(
    "e([+-]?[0-9]+)",
    " %*% 10^{\\1}",
    formatC(pvalue, format = "e", digits = 1)
  )
)]

pval_matrix <- make_matrix(res_dars_sub, "pval_fmt")

# -----------------------------------------------------------------------------
# 5. Heatmap colour scale
# -----------------------------------------------------------------------------

log_min <- min(or_matrix, na.rm = TRUE)
log_max <- max(or_matrix, na.rm = TRUE)

mat_color <- colorRamp2(
  c(log_min, 0, log_max),
  c("skyblue", "white", "red")
)

# -----------------------------------------------------------------------------
# 6. Construct heatmap
# -----------------------------------------------------------------------------

ht <- Heatmap(
  or_matrix,
  name = "OR",
  col = mat_color,
  na_col = "white",

  cluster_rows = FALSE,
  cluster_columns = FALSE,

  row_names_side = "left",
  row_names_max_width = unit(12, "cm"),

  column_labels = col_names_labs,
  column_split = col_split_group,
  column_names_rot = 45,

  column_title_gp = gpar(fontface = "italic"),
  column_names_gp = gpar(fontsize = 14, fontface = "bold"),
  row_names_gp = gpar(fontsize = 14, fontface = "bold"),

  border_gp = gpar(col = "black", lty = 1, lwd = 1.5),
  rect_gp = gpar(col = "black", lwd = 0.2),

  cell_fun = function(j, i, x, y, width, height, fill) {

    ci_val <- ci_matrix[i, j]
    pval_val <- pval_matrix[i, j]
    star_val <- pstar_matrix[i, j]

    # Log OR and confidence interval
    if (!is.na(ci_val)) {
      grid.text(
        ci_val,
        x,
        y + unit(2.5, "mm"),
        gp = gpar(fontsize = 12)
      )
    }

    # P-value
    if (!is.na(pval_val)) {
      grid.text(
        parse(text = pval_val),
        x,
        y - unit(2.0, "mm"),
        gp = gpar(fontsize = 12)
      )
    }

    # Significance stars
    if (!is.na(star_val) && nzchar(star_val)) {
      grid.text(
        star_val,
        x,
        y - unit(6.0, "mm"),
        gp = gpar(fontsize = 12)
      )
    }
  },

  column_gap = unit(2, "mm"),
  row_gap = unit(1, "mm"),
  row_title = NULL,
  show_column_names = TRUE
)

# -----------------------------------------------------------------------------
# 7. Save figure
# -----------------------------------------------------------------------------

pdf(output_file, width = 18, height = 12)

draw(ht)

dev.off()

message("Heatmap saved: ", output_file)


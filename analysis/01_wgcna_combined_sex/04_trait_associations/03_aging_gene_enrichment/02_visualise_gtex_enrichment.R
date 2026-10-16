# 02_visualise_gtex_enrichment.R
# Reproduce the two GTEx DEG enrichment heatmaps for Supplementary Fig. S17.
# Pairwise (Heteroskid: aSCPs/aVCPs) and module-level (WGCNA: aSCMs/aVCMs)
# panels are exported separately; the manuscript panels were assembled in PPT.
# Run from the repository root after 01_run_aging_gene_enrichment.R.
# Plotting and statistical display follow the original 17b script.

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(ComplexHeatmap)
  library(circlize)
  library(reshape2)
  library(grid)
})

results_dir <- "/results/01_wgcna_combined_sex"
input_file <- file.path(results_dir, "10a_GTEx_age_DEGs_enrichment.xlsx")
output_dir <- "/figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) stop("Missing GTEx enrichment results: ", input_file)
res_twas_wang <- as.data.table(read_xlsx(input_file))

col_mat_order <- c(
  "sex-biased (male) Stable", "sex-biased (female) Stable", "sex-common Stable",
  "sex-biased (male) Changing", "sex-biased (female) Changing", "sex-common Changing"
)
col_names_labs <- c(
  "sex-biased (male)", "sex-biased (female)", "sex-common",
  "sex-biased (male)", "sex-biased (female)", "sex-common"
)
col_split_group <- factor(rep(c("Stable", "Changing"), each = 3),
                          levels = c("Stable", "Changing"))

plot_gtex_heatmap <- function(sigGroup, col_title_lab, filename) {
  # Original 17b used %like% to select Heteroskid or WGCNA and Group == DEGs.
  dt <- copy(res_twas_wang[grepl(sigGroup, analysis, fixed = TRUE) & Group == "DEGs"])
  if (!nrow(dt)) stop("No GTEx DEG rows for analysis: ", sigGroup)

  # Original figure displays log(OR) and log(CI), not raw OR.
  dt[, `:=`(OR = log(OR),
            CI_lower = log(CI_lower),
            CI_upper = log(CI_upper))]
  dt[, OR := round(OR, 2)]

  # Original label convention: analysis names split into prefix and Stable/Changing.
  parts <- tstrsplit(dt$analysis, "_", fixed = TRUE)
  if (length(parts) < 2L) stop("Expected analysis labels like Heteroskid_Stable: ", sigGroup)
  dt[, sigGroup_label := parts[[2]]]
  dt[, col_id := paste(sex_class, sigGroup_label)]

  missing_cols <- setdiff(col_mat_order, unique(dt$col_id))
  if (length(missing_cols)) {
    stop("Missing expected GTEx categories for ", sigGroup, ": ",
         paste(missing_cols, collapse = ", "),
         ". Check the original enrichment labels before changing the plot.")
  }
  if (anyDuplicated(dt[, .(tissue, col_id)])) {
    stop("Duplicate tissue/category rows in GTEx results for ", sigGroup)
  }

  dt[, OR_CI_label := fifelse(
    is.na(OR) | is.na(CI_lower) | is.na(CI_upper),
    NA_character_,
    paste0(round(OR, 2), "\n(", round(CI_lower, 2), "-", round(CI_upper, 2), ")")
  )]
  dt[, p_stars := fifelse(pvalue < 0.001, "***",
                         fifelse(pvalue < 0.01, "**",
                                 fifelse(pvalue < 0.05, "*", NA_character_)))]

  or_matrix <- reshape2::acast(dt, tissue ~ col_id, value.var = "OR")
  ci_matrix <- reshape2::acast(dt, tissue ~ col_id, value.var = "OR_CI_label")
  pstar_matrix <- reshape2::acast(dt, tissue ~ col_id, value.var = "p_stars")
  or_matrix <- or_matrix[, col_mat_order, drop = FALSE]
  ci_matrix <- ci_matrix[rownames(or_matrix), col_mat_order, drop = FALSE]
  pstar_matrix <- pstar_matrix[rownames(or_matrix), col_mat_order, drop = FALSE]
  ci_matrix[ci_matrix == "NA"] <- NA_character_

  log_min <- min(or_matrix, na.rm = TRUE)
  log_max <- max(or_matrix, na.rm = TRUE)
  if (!is.finite(log_min) || !is.finite(log_max)) {
    stop("No finite log odds ratios available for ", sigGroup)
  }
  # Original blue-white-red scale, centred at zero.
  mat_color <- circlize::colorRamp2(
    c(log_min, 0, log_max), c("skyblue", "white", "red")
  )

  ht <- ComplexHeatmap::Heatmap(
    or_matrix,
    name = "OR",  # Retain original legend label; underlying values are log(OR).
    col = mat_color,
    na_col = "grey90",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    row_names_max_width = grid::unit(12, "cm"),
    column_labels = col_names_labs,
    column_split = col_split_group,
    column_names_rot = 45,
    column_title = col_title_lab,
    column_title_gp = grid::gpar(fill = c("white", "white"), font = 2),
    border_gp = grid::gpar(col = "black", lty = 1, lwd = 1.5),
    rect_gp = grid::gpar(col = "black", lwd = 0.2),
    row_names_gp = grid::gpar(fontsize = 14, fontface = "bold"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      ci_val <- ci_matrix[i, j]
      if (!is.na(ci_val)) {
        grid::grid.text(ci_val, x, y + grid::unit(1, "mm"),
                        gp = grid::gpar(fontsize = 10))
      }
      star_val <- pstar_matrix[i, j]
      if (!is.na(star_val)) {
        grid::grid.text(star_val, x, y - grid::unit(4.5, "mm"),
                        gp = grid::gpar(fontsize = 10))
      }
    },
    column_gap = grid::unit(2, "mm"),
    column_names_gp = grid::gpar(fontsize = 14, fontface = "bold"),
    row_gap = grid::unit(1, "mm"),
    row_title = NULL,
    show_column_names = TRUE
  )

  pdf(filename, width = 15, height = max(8, 0.45 * nrow(or_matrix) + 3),
      useDingbats = FALSE)
  tryCatch(ComplexHeatmap::draw(ht), finally = dev.off())
  message("Saved: ", filename)
  invisible(ht)
}

plot_gtex_heatmap(
  sigGroup = "Heteroskid",
  col_title_lab = c("aSCPs", "aVCPs"),
  filename = file.path(output_dir, "Fig_S17a_GTEx_pairwise.pdf")
)
plot_gtex_heatmap(
  sigGroup = "WGCNA",
  col_title_lab = c("aSCMs", "aVCMs"),
  filename = file.path(output_dir, "Fig_S17b_GTEx_modules.pdf")
)



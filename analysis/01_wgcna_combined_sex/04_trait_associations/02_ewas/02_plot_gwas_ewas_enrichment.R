
# =============================================================================
# Purpose: Plot Figure 6a showing GWAS and EWAS Catalog enrichment
#          across age-stable and age-variable methylation groups.
#
# Inputs:  MAGMA GWAS enrichment results, GWAS trait names, EWAS Catalog
#          enrichment results, and manually curated EWAS trait themes.
# Outputs: Figure 6a GWAS/EWAS enrichment heatmap.
# =============================================================================

library(data.table)
library(readxl)
library(circlize)
library(ComplexHeatmap)
library(grid)

# -------------------------------------------------------------------------
# 1. Input and output paths
# -------------------------------------------------------------------------

res_dir <- "/results/01_wgcna_combined_sex"
trait_dir <- "data/gwas"
figure_dir <- "figures"

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# 2. Load GWAS enrichment results and trait names
# -------------------------------------------------------------------------

gwas_trait_names <- fread(
  file.path(trait_dir, "MAGMA_GWAS_traits_v2.txt")
)[, .(Trait, DatasetID)]

res_gwas <- fread(
  file.path(res_dir, "comb_gwas_gene_overlap_2.5percent.tsv")
)

res_gwas <- merge(
  gwas_trait_names,
  res_gwas,
  by.x = "DatasetID",
  by.y = "trait"
)

# Preserve original selection:
# Top five associations per analysis, using MHC == "yes".
res_gwas_sub_top <- unique(
  res_gwas[MHC == "yes"][
    ,
    .SD[order(pvalue)][1:5],
    by = analysis
  ]
)

# -------------------------------------------------------------------------
# 3. Load EWAS Catalog enrichment and curated trait themes
# -------------------------------------------------------------------------

res_ewas <- as.data.table(
  read_xlsx(
    file.path(res_dir, "12a_ewas_catalog_top_overlaps.xlsx")
  )
)

# Retain original top-five EWAS selection for reference
res_ewas_topN <- res_ewas[
  ,
  .SD[order(pvalue)][1:5],
  by = analysis
][
  ,
  !c("FDR", "CI_range"),
  with = FALSE
]

# Curated trait themes used in the original Figure 6a
res_ewas_curated <- as.data.table(
  read_xlsx(
    file.path(res_dir, "12b_ewas_catalog_top_10_sub_traits.xlsx")
  )
)

res_ewas_curated[
  ,
  study_id := sub("_.*", "", ewas_study)
]

res_ewas_curated[
  ,
  trait_theme := paste0(study_id, "_", trait_theme)
]

# Preserve original exclusions
res_ewas_curated <- res_ewas_curated[
  !trait_theme %in% c(
    "26864933_Smoking",
    "28603561_Cleft_lip_palate"
  )
]

# Select the lowest-p-value study per analysis and trait theme
res_ewas_sub_top <- res_ewas_curated[
  order(pvalue),
  .SD[1],
  by = .(analysis, trait_theme)
][order(analysis)]

# -------------------------------------------------------------------------
# 4. Combine GWAS and EWAS enrichment results
# -------------------------------------------------------------------------

gwas_plot <- res_gwas_sub_top[
  ,
  .(
    analysis,
    DatasetID,
    Trait,
    pvalue,
    OR,
    CI_lower,
    CI_upper
  )
][
  ,
  group := "GWAS"
]

ewas_plot <- res_ewas_sub_top[
  ,
  .(
    analysis,
    ewas_study,
    trait_theme,
    pvalue,
    OR,
    CI_lower,
    CI_upper
  )
][
  ,
  group := "MWAS"
]

res_combined <- rbind(
  gwas_plot,
  ewas_plot,
  use.names = FALSE
)

# Preserve original top-five selection per analysis and group
res_combined <- res_combined[
  ,
  .SD[order(pvalue)][1:5],
  by = .(analysis, group)
]

dt_plot <- res_combined[
  group %in% c("GWAS", "MWAS") & !is.na(OR)
]

dt_plot[
  ,
  row_id := paste0(group, ":", Trait)
]

# -------------------------------------------------------------------------
# 5. Prepare heatmap matrices
# -------------------------------------------------------------------------

# Log-transform odds ratios and confidence intervals
dt_plot[
  ,
  `:=`(
    OR = log(OR),
    CI_lower = log(CI_lower),
    CI_upper = log(CI_upper)
  )
]

dt_plot[, OR := round(OR, 2)]

column_labs <- c("aSCMs", "aSCPs", "aVCMs", "aVCPs")

col_order <- c(
  "WGCNA_Stable",
  "Heteroskid_Stable",
  "WGCNA_Changing",
  "Heteroskid_Changing"
)

# Odds-ratio matrix
or_matrix <- reshape2::acast(
  dt_plot,
  row_id ~ analysis,
  value.var = "OR"
)

# Preserve original row grouping and ordering
group_order <- sub(":.*", "", rownames(or_matrix))

group_factor <- factor(
  group_order,
  levels = c("GWAS", "MWAS", "Longevity exp", "DNAm Clocks")
)

row_reorder <- order(group_factor)

or_matrix_ordered <- or_matrix[row_reorder, col_order]

# Confidence-interval labels
dt_plot[
  ,
  OR_CI_label := ifelse(
    is.na(OR) | is.na(CI_lower) | is.na(CI_upper),
    "NA",
    paste0(
      round(OR, 2),
      " (",
      round(CI_lower, 2),
      "-",
      round(CI_upper, 2),
      ")"
    )
  )
]

ci_matrix <- reshape2::acast(
  dt_plot,
  row_id ~ analysis,
  value.var = "OR_CI_label"
)

ci_matrix[ci_matrix == "NA"] <- NA

ci_matrix_ordered <- ci_matrix[row_reorder, col_order]

# Scientific notation for p-values
dt_plot[
  ,
  pval_fmt := ifelse(
    is.na(pvalue),
    NA_character_,
    gsub(
      "e([+-]?[0-9]+)",
      " %*% 10^{\\1}",
      formatC(pvalue, format = "e", digits = 1)
    )
  )
]

pval_matrix <- reshape2::acast(
  dt_plot,
  row_id ~ analysis,
  value.var = "pval_fmt"
)

pval_matrix_ordered <- pval_matrix[row_reorder, col_order]

# Significance-star matrix (retained for consistency)
dt_plot[
  ,
  p_stars := cut(
    pvalue,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", ""),
    right = FALSE
  )
]

pstar_matrix <- reshape2::acast(
  dt_plot,
  row_id ~ analysis,
  value.var = "p_stars"
)

pstar_matrix_ordered <- pstar_matrix[row_reorder, col_order]

# -------------------------------------------------------------------------
# 6. Heatmap annotations and colour scale
# -------------------------------------------------------------------------

row_names_labs <- sub(
  ".*:",
  "",
  rownames(or_matrix_ordered)
)

row_split_labs <- sub(
  ":.*",
  "",
  rownames(or_matrix_ordered)
)

row_split_labs <- factor(
  row_split_labs,
  levels = unique(row_split_labs)
)

col_split_labs <- factor(
  c("Stable", "Stable", "Changing", "Changing"),
  levels = c("Stable", "Changing")
)

left_annot_labs <- c("GWAS", "MWAS")

# Preserve original skyblue-white-red colour scale
log_max <- max(or_matrix_ordered, na.rm = TRUE)

mat_col <- colorRamp2(
  c(-1, 0, log_max),
  c("skyblue", "white", "red")
)

# -------------------------------------------------------------------------
# 7. Draw Figure 6a
# -------------------------------------------------------------------------

heatmap_fun <- function(or_mat, ci_mat, pval_mat, pstar_mat, mat_cols) {

  Heatmap(
    or_mat,
    name = "OR",
    col = mat_cols,
    na_col = "grey90",

    cluster_rows = FALSE,
    cluster_columns = FALSE,

    row_labels = row_names_labs,
    row_names_side = "left",
    column_labels = column_labs,

    column_split = col_split_labs,
    row_split = row_split_labs,

    column_names_rot = 45,
    row_names_max_width = unit(25, "cm"),

    column_title_gp = gpar(
      fill = c("white", "white"),
      font = 2,
      fontsize = 15
    ),

    border_gp = gpar(
      col = "black",
      lty = 1,
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

    column_gap = unit(2, "mm"),

    column_names_gp = gpar(
      fontsize = 14,
      fontface = "bold"
    ),

    row_title = NULL,
    column_title = NULL,

    cell_fun = function(j, i, x, y, width, height, fill) {

      ci_val <- ci_mat[i, j]
      pval_val <- pval_mat[i, j]

      # Confidence interval above cell centre
      if (!is.na(ci_val)) {
        grid.text(
          ci_val,
          x,
          y + unit(2.5, "mm"),
          gp = gpar(fontsize = 11)
        )
      }

      # P-value below cell centre
      if (!is.na(pval_val)) {
        grid.text(
          parse(text = pval_val),
          x,
          y - unit(1.0, "mm"),
          gp = gpar(fontsize = 11)
        )
      }
    },

    left_annotation = rowAnnotation(
      cpgAnnot = anno_block(
        labels = left_annot_labs,
        labels_gp = gpar(fontsize = 14)
      )
    )
  )
}

# Render heatmap
ht <- heatmap_fun(
  or_mat = or_matrix_ordered,
  ci_mat = ci_matrix_ordered,
  pval_mat = pval_matrix_ordered,
  pstar_mat = pstar_matrix_ordered,
  mat_cols = mat_col
)

draw(ht)

# -------------------------------------------------------------------------
# 8. Export figure
# -------------------------------------------------------------------------

pdf(
  file.path(figure_dir, "Figure_6a_GWAS_EWAS_enrichment.pdf"),
  width = 14,
  height = 12
)

draw(ht)
dev.off()


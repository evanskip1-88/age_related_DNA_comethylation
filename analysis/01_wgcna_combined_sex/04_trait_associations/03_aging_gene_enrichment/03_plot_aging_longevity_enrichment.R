
# =============================================================================
# Purpose: Visualise enrichment of longevity experiments and DNAm clocks.
# Inputs:  Longevity intervention and epigenetic clock enrichment results.
# Outputs: Figure 6b enrichment heatmap (PDF).
# =============================================================================

library(data.table)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(grid)

# Run from the repository root
results_dir <- "/results/01_wgcna_combined_sex"
figure_dir <- "/figures"

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# 1. Load enrichment results
# -------------------------------------------------------------------------

res_all_intervent <- as.data.table(read_xlsx(
  file.path(results_dir, "10b_intervention_enrich.xlsx")
))
res_all_intervent[, group := "Longevity exp"]

res_clock <- as.data.table(read_xlsx(
  file.path(results_dir, "01a_all_clock_enrichements.xlsx")
))
res_clock[, group := "DNAm clocks"]

# -------------------------------------------------------------------------
# 2. Combine longevity experiments and DNAm clocks
# -------------------------------------------------------------------------

# Retain the original clock exclusions
res_clock_sub <- res_clock[
  !Clock %in% c("S3.1_Clock 1", "S3.2_Clock 2"),
  .(analysis, group, Study, Clock, pvalue, OR, CI_lower, CI_upper)
]

# Preserve the original complete-case filtering for interventions
res_intervent_sub <- na.omit(
  res_all_intervent[
    ,
    .(analysis, group, Dataset, Group, pvalue, OR, CI_lower, CI_upper)
  ]
)

# Positional binding is intentional: Study/Dataset and Clock/Group
# represent the corresponding identifier and display-label columns.
res_all_pheno_overlaps <- rbind(
  res_clock_sub,
  res_intervent_sub,
  use.names = FALSE
)

# Preserve the original clock label
res_all_pheno_overlaps[
  Clock == "S3.3_Clock 3",
  Clock := "Pan_mammalian_clock3"
]

# -------------------------------------------------------------------------
# 3. Prepare heatmap data
# -------------------------------------------------------------------------

dt_merge_final_sub <- copy(
  res_all_pheno_overlaps[
    !is.na(OR)
  ]
)

# In the original combined table, Clock holds the display labels
# for both DNAm clocks and longevity experiments.
dt_merge_final_sub[
  ,
  row_id := paste0(group, ":", Clock)
]

# Log-transform odds ratios and confidence intervals
dt_merge_final_sub[
  ,
  `:=`(
    OR = round(log(OR), 2),
    CI_lower = log(CI_lower),
    CI_upper = log(CI_upper)
  )
]

column_labs <- c("aSCMs", "aSCPs", "aVCMs", "aVCPs")

col_order <- c(
  "WGCNA_Stable",
  "Heteroskid_Stable",
  "WGCNA_Changing",
  "Heteroskid_Changing"
)

# Odds-ratio matrix
or_matrix <- reshape2::acast(
  dt_merge_final_sub,
  row_id ~ analysis,
  value.var = "OR"
)

# Preserve the original row-ordering procedure
group_order <- sub(":.*", "", rownames(or_matrix))

new_order <- c(
  "GWAS",
  "MWAS",
  "Longevity exp",
  "DNAm Clocks"
)

group_factor <- factor(group_order, levels = new_order)
row_reorder <- order(group_factor)

or_matrix_ordered <- or_matrix[
  row_reorder,
  col_order,
  drop = FALSE
]

# -------------------------------------------------------------------------
# 4. Confidence intervals and p-value matrices
# -------------------------------------------------------------------------

dt_merge_final_sub[
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
  dt_merge_final_sub,
  row_id ~ analysis,
  value.var = "OR_CI_label"
)

ci_matrix[ci_matrix == "NA"] <- NA

ci_matrix_ordered <- ci_matrix[
  row_reorder,
  col_order,
  drop = FALSE
]

# Format p-values for scientific notation
dt_merge_final_sub[
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
  dt_merge_final_sub,
  row_id ~ analysis,
  value.var = "pval_fmt"
)

pval_matrix_ordered <- pval_matrix[
  row_reorder,
  col_order,
  drop = FALSE
]

# -------------------------------------------------------------------------
# 5. Heatmap labels and colours
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

left_annot_labs <- c(
  "Longevity experiments",
  "DNAm clocks"
)

# Original colour scale
log_max <- max(or_matrix_ordered, na.rm = TRUE)

mat_col <- colorRamp2(
  c(-1, 0, log_max),
  c("skyblue", "white", "red")
)

# -------------------------------------------------------------------------
# 6. Plot enrichment heatmap
# -------------------------------------------------------------------------

heatmap_fun <- function(or_mat, ci_mat, pval_mat, mat_cols) {

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

      if (!is.na(ci_val)) {
        grid.text(
          ci_val,
          x,
          y + unit(2.5, "mm"),
          gp = gpar(fontsize = 11)
        )
      }

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

ht <- heatmap_fun(
  or_mat = or_matrix_ordered,
  ci_mat = ci_matrix_ordered,
  pval_mat = pval_matrix_ordered,
  mat_cols = mat_col
)

# -------------------------------------------------------------------------
# 7. Export heatmap
# -------------------------------------------------------------------------

output_file <- file.path(
  figure_dir,
  "Figure_6b_aging_longevity_enrichment.pdf"
)

pdf(output_file, width = 14, height = 12)
draw(ht)
dev.off()

message("Saved: ", output_file)


#==============================================================================
# Visualise transcription factor binding-site enrichment
#
# Input: ENCODE TFBS enrichment results (no permutations)
# Outputs:
#   1. Dot plot for all four CpG groups
#   2. Heatmap for aSCPs and aVCPs
#==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

#------------------------------------------------------------------------------
# Paths
#------------------------------------------------------------------------------

results_dir <- "/results/01_wgcna_combined_sex"
figure_dir <- "/figures"

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

#------------------------------------------------------------------------------
# Load enrichment results
#------------------------------------------------------------------------------

tfbs <- fread(
  file.path(results_dir, "01a_tfbs_enrich_no_perm.tsv")
)

# Retain significant TFBS enrichment results
# Preserve the original unadjusted p-value threshold.

tfbs_sub <- copy(tfbs)[
  !is.na(pval) & pval < 0.05
]

tfbs_sub[
  ,
  SigGroup := factor(
    SigGroup,
    levels = c("Stable", "Changing")
  )
]

tfbs_sub[
  ,
  SigGroup2 := fcase(
    Dataset == "wgcna_Stable", "aSCMs",
    Dataset == "wgcna_Changing", "aVCMs",
    Dataset == "het_Stable", "aSCPs",
    Dataset == "het_Changing", "aVCPs",
    default = NA_character_
  )
]

#------------------------------------------------------------------------------
# Figure 1: TFBS enrichment dot plot
#------------------------------------------------------------------------------

dot_data <- copy(tfbs_sub)[
  !is.na(OR) & !is.na(CI_lower) & !is.na(CI_upper)
]

dot_data[
  ,
  TF_group := paste(SigGroup, TF, sep = "_")
]

dot_data[
  ,
  TF_plot := factor(
    TF_group,
    levels = unique(TF_group[order(OR)])
  ),
  by = SigGroup
]

p_dot <- ggplot(
  dot_data,
  aes(x = OR, y = TF_plot)
) +
  geom_point(aes(size = overlap)) +
  geom_errorbar(
    aes(xmin = CI_lower, xmax = CI_upper),
    orientation = "y",
    width = 0.15,
    linewidth = 0.6
  ) +
  scale_y_discrete(
    labels = function(x) {
      sub("^(Stable|Changing)_", "", x)
    }
  ) +
  scale_size_continuous(range = c(2, 8)) +
  facet_wrap(~SigGroup2, scales = "free_y", nrow = 2) +
  labs(
    x = "Odds Ratio",
    y = "Transcription factor",
    size = "Overlap"
  ) +
  theme_bw(base_size = 16) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(
      size = 14, colour = "black", face = "bold"
    ),
    axis.text.y = element_text(
      size = 14, colour = "black", face = "bold"
    ),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(linewidth = 1)
  )

ggsave(
  file.path(figure_dir, "tfbs_enrichment_dotplot.pdf"),
  plot = p_dot,
  width = 14,
  height = 12,
  limitsize = FALSE
)

#------------------------------------------------------------------------------
# Figure 2: TFBS enrichment heatmap
#
# Original visualization focused on Heteroskedasticity:
# aSCPs (Stable) versus aVCPs (Changing).
#------------------------------------------------------------------------------

analysis_type <- "Heteroskedasticity"

heatmap_data <- copy(tfbs_sub)[
  Analysis == analysis_type
]

# Convert odds ratios and confidence intervals to log scale
for (col in c("OR", "CI_lower", "CI_upper")) {
  set(
    heatmap_data,
    j = col,
    value = log(heatmap_data[[col]])
  )
}

heatmap_data[
  ,
  OR := round(OR, 2)
]

# Original labels use log(OR) and log confidence intervals
heatmap_data[
  ,
  OR_CI_label := fifelse(
    is.na(OR) | is.na(CI_lower) | is.na(CI_upper),
    NA_character_,
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

heatmap_data[
  ,
  p_stars := fcase(
    pval < 0.001, "***",
    pval < 0.01, "**",
    pval < 0.05, "*",
    default = NA_character_
  )
]

heatmap_data[
  ,
  pval_fmt := fifelse(
    is.na(pval),
    NA_character_,
    gsub(
      "e([+-]?[0-9]+)",
      " %*% 10^{\\1}",
      formatC(pval, format = "e", digits = 1)
    )
  )
]

#------------------------------------------------------------------------------
# Construct heatmap matrices
#------------------------------------------------------------------------------

column_order <- c("het_Stable", "het_Changing")
column_labels <- c("aSCPs", "aVCPs")

make_matrix <- function(value_col) {

  wide <- dcast(
    heatmap_data,
    TF ~ Dataset,
    value.var = value_col,
    fill = NA
  )

  for (column in column_order) {
    if (!column %in% names(wide)) {
      wide[, (column) := NA]
    }
  }

  matrix_data <- as.matrix(
    wide[, ..column_order]
  )

  rownames(matrix_data) <- wide$TF

  matrix_data
}

or_matrix <- make_matrix("OR")
ci_matrix <- make_matrix("OR_CI_label")
pstar_matrix <- make_matrix("p_stars")
pval_matrix <- make_matrix("pval_fmt")

#------------------------------------------------------------------------------
# Heatmap colour scale
#------------------------------------------------------------------------------

log_max <- max(or_matrix, na.rm = TRUE)

mat_color <- colorRamp2(
  c(-1, 0, log_max),
  c("skyblue", "white", "red")
)

#------------------------------------------------------------------------------
# Draw heatmap
#------------------------------------------------------------------------------

tfbs_heatmap <- Heatmap(
  or_matrix,
  name = "OR",
  col = mat_color,
  na_col = "white",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  row_names_max_width = unit(12, "cm"),
  column_labels = column_labels,
  column_split = factor(
    c("Stable", "Changing"),
    levels = c("Stable", "Changing")
  ),
  column_names_rot = 45,
  column_title = c("aSCPs", "aVCPs"),
  column_title_gp = gpar(font = 2),
  border_gp = gpar(col = "black", lty = 1, lwd = 1.5),
  rect_gp = gpar(col = "black", lwd = 0.2),
  row_names_gp = gpar(fontsize = 14, fontface = "bold"),
  column_gap = unit(2, "mm"),
  row_gap = unit(1, "mm"),
  row_title = NULL,
  show_column_names = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {

    ci_val <- ci_matrix[i, j]
    pval_val <- pval_matrix[i, j]
    star_val <- pstar_matrix[i, j]

    if (!is.na(ci_val)) {
      grid.text(
        ci_val,
        x,
        y + unit(4.5, "mm"),
        gp = gpar(fontsize = 12)
      )
    }

    if (!is.na(pval_val)) {
      grid.text(
        parse(text = pval_val),
        x,
        y - unit(1, "mm"),
        gp = gpar(fontsize = 12)
      )
    }

    if (!is.na(star_val) && nzchar(star_val)) {
      grid.text(
        star_val,
        x,
        y - unit(6, "mm"),
        gp = gpar(fontsize = 12)
      )
    }
  }
)

pdf(
  file.path(figure_dir, "tfbs_enrichment_heatmap.pdf"),
  width = 12,
  height = max(8, nrow(or_matrix) * 0.5)
)

draw(tfbs_heatmap)
dev.off()


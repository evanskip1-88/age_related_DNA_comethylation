# ============================================================
# Midlife sex-biased module enrichment heatmap
# Comparison: Age2 vs Age4
#
# Selected modules: 59, 43, 63
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(reshape2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# Directories
results_path <- "/results/02_wgcna_within_sex"
output_path <- "/results/02_wgcna_within_sex"

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Load manually annotated GREAT results
# ------------------------------------------------------------

go <- setDT(
  read_xlsx(
    file.path(
      results_path,
      "01a_age2v4_topN_BP_per_mod_by_sex_modref_M.xlsx"
    )
  )
)

pheno <- setDT(
  read_xlsx(
    file.path(
      results_path,
      "01a_age2v4_topN_pheno_per_mod_by_sex_modref_M.xlsx"
    )
  )
)

go_pheno <- rbindlist(
  list(go, pheno),
  use.names = TRUE,
  fill = TRUE
)

# Select the strongest enrichment for each biological theme
go_pheno_grp <- go_pheno[
  ,
  .SD[order(Hyper_Adjp_BH)][1],
  by = .(Group, moduleColor, Biol_theme)
]

# ------------------------------------------------------------
# 2. Attach module IDs and sex-bias classification
# ------------------------------------------------------------
dt_age2v4 <- read_xlsx(
  file.path(
    results_path,
    "01a_sex_specific_mods_modref_M_age2v4.xlsx"
  )
)
setDT(dt_age2v4)

dt_age2v4 <- unique(
  dt_age2v4[
    , .(
      moduleColor,
      moduleNumID,
      sexChangeGroup2v4
    )
  ]
)

go_pheno_grp <- merge(
  dt_age2v4,
  go_pheno_grp,
  by = "moduleColor"
)

go_pheno_grp[
  , moduleNumID := paste0("module_", moduleNumID)
]

go_pheno_grp[
  , col_id := paste0(
    moduleNumID,
    "_",
    sexChangeGroup2v4
  )
]

go_pheno_grp[
  , Hyper_Adjp_BH_log10 := -log10(Hyper_Adjp_BH)
]

# ------------------------------------------------------------
# 3. Construct matrices
# ------------------------------------------------------------

col_order <- c(
  "module_59_change_Male",
  "module_43_change_Female",
  "module_63_change_Female"
)

categories <- c(
  "BP",
  "Mouse pheno sKO",
  "Human pheno"
)

make_matrix <- function(category, value_col) {

  dt_sub <- go_pheno_grp[Group == category]

  mat <- reshape2::acast(
    dt_sub,
    Biol_theme ~ col_id,
    value.var = value_col
  )

  # Ensure all three modules are represented
  missing_cols <- setdiff(col_order, colnames(mat))

  if (length(missing_cols) > 0) {
    missing <- matrix(
      NA,
      nrow = nrow(mat),
      ncol = length(missing_cols),
      dimnames = list(rownames(mat), missing_cols)
    )

    mat <- cbind(mat, missing)
  }

  mat[, col_order, drop = FALSE]
}

# Numerical enrichment matrices
p_matrix_BP <- make_matrix(
  "BP",
  "Hyper_Adjp_BH_log10"
)

p_matrix_mpko <- make_matrix(
  "Mouse pheno sKO",
  "Hyper_Adjp_BH_log10"
)

p_matrix_hp <- make_matrix(
  "Human pheno",
  "Hyper_Adjp_BH_log10"
)

combine_mat_log10 <- rbind(
  p_matrix_BP,
  p_matrix_mpko,
  p_matrix_hp
)

# ------------------------------------------------------------
# 4. Format adjusted p-values for cell labels
# ------------------------------------------------------------

go_pheno_grp[
  , pval_fmt := fifelse(
    is.na(Hyper_Adjp_BH),
    NA_character_,
    gsub(
      "e([+-]?[0-9]+)",
      " %*% 10^{\\1}",
      formatC(
        Hyper_Adjp_BH,
        format = "e",
        digits = 1
      )
    )
  )
]

pval_matrix_BP <- make_matrix("BP", "pval_fmt")
pval_matrix_mpko <- make_matrix("Mouse pheno sKO", "pval_fmt")
pval_matrix_hp <- make_matrix("Human pheno", "pval_fmt")

pval_combine_mat <- rbind(
  pval_matrix_BP,
  pval_matrix_mpko,
  pval_matrix_hp
)

# ------------------------------------------------------------
# 5. Define row and column annotations
# ------------------------------------------------------------

row_split_labs <- factor(
  c(
    rep("Biol. process", nrow(p_matrix_BP)),
    rep("Mouse pheno sKO", nrow(p_matrix_mpko)),
    rep("Human pheno", nrow(p_matrix_hp))
  ),
  levels = c(
    "Biol. process",
    "Mouse pheno sKO",
    "Human pheno"
  )
)

col_split_labs <- factor(
  c(
    "Male biased",
    "Female biased",
    "Female biased"
  ),
  levels = c(
    "Male biased",
    "Female biased"
  )
)

col_names_labs <- c(
  "module_59",
  "module_43",
  "module_63"
)

mat_col <- circlize::colorRamp2(
  c(5, 15, 40, 60, 125),
  c(
    "lightyellow",
    "gold",
    "orange",
    "orangered",
    "tomato"
  )
)

# ------------------------------------------------------------
# 6. Generate heatmap
# ------------------------------------------------------------

ht <- Heatmap(
  combine_mat_log10,
  name = "-log10(Adj_Pval)",
  col = mat_col,
  na_col = "grey90",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  column_labels = col_names_labs,
  column_split = col_split_labs,
  row_split = row_split_labs,
  column_names_rot = 45,
  row_names_max_width = unit(25, "cm"),
  column_title_gp = gpar(font = 2, fontsize = 15),
  border_gp = gpar(col = "black", lty = 1, lwd = 1.5),
  rect_gp = gpar(col = "black", lwd = 0.2),
  row_names_gp = gpar(fontsize = 12, fontface = "bold"),
  column_gap = unit(2, "mm"),
  column_names_gp = gpar(fontsize = 14, fontface = "bold"),
  row_title = NULL,
  show_column_names = TRUE,
  show_heatmap_legend = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {

    pval_val <- pval_combine_mat[i, j]

    if (!is.na(pval_val)) {
      grid.text(
        parse(text = pval_val),
        x,
        y - unit(0.1, "mm"),
        gp = gpar(fontsize = 14)
      )
    }
  },
  left_annotation = rowAnnotation(
    cpgAnnot = anno_block(
      gp = gpar(
        fill = c(
          "skyblue",
          "lightgreen",
          "slateblue"
        )
      ),
      labels = c(
        "Biological process",
        "Mouse pheno sKO",
        "H.pheno"
      ),
      labels_gp = gpar(fontsize = 14)
    )
  )
)

pdf(
  file.path(
    output_path,
    "midlife_rGREAT_biological_themes.pdf"
  ),
  width = 12,
  height = 10
)

draw(ht)
dev.off()

cat("Midlife GREAT heatmap completed.\n")

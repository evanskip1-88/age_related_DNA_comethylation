# ============================================================
# Biological themes of sex-biased co-methylation modules
# Supplementary Figure S21 / S22
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(reshape2)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# Directories
results_path <- "/results/02_wgcna_within_sex"
output_path <- "/results/02_wgcna_within_sex"

mod_ref <- "M"

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Load reference module labels
# ------------------------------------------------------------

dt_mod_labs <- fread(
  file.path(
    results_path,
    paste0(
      "precorr_gs_qc_chr_all_refage1_MinMod20_",
      mod_ref, ".modLabs"
    )
  )
)

dt_mod_map <- unique(
  dt_mod_labs[, .(moduleNumID, moduleColors)]
)

# ------------------------------------------------------------
# 2. Load manually curated biological themes
# ------------------------------------------------------------

dt_themes <- setDT(
  read_xlsx(
    file.path(
      results_path,
      paste0(
        "03a_top5_BP_per_mod_by_sex_modref_",
        mod_ref, ".xlsx"
      )
    ),
    sheet = 2
  )
)

dt_themes[
  , Hyper_Adjp_BH_log10 := -log10(Hyper_Adjp_BH)
]

# Attach module IDs
dt_themes <- merge(
  dt_mod_map,
  dt_themes,
  by.x = "moduleColors",
  by.y = "moduleColor",
  all.y = TRUE
)

dt_themes[
  , moduleNumID := paste0("module_", moduleNumID)
]

# ------------------------------------------------------------
# 3. Select the top enrichment per theme and module
# ------------------------------------------------------------

dt_themes_sub <- dt_themes[
  moduleColors != "All",
  .SD[order(Hyper_Adjp_BH)][1],
  by = .(
    Group,
    sigGroup,
    moduleColors,
    Biological_Theme,
    ModRef
  )
]

dt_themes_sub[
  , col_id := paste0(sigGroup, "_", moduleNumID)
]

setorder(
  dt_themes_sub,
  -sigGroup,
  moduleNumID
)

# ------------------------------------------------------------
# 4. Construct enrichment matrix
# ------------------------------------------------------------

p_matrix <- reshape2::acast(
  dt_themes_sub,
  Biological_Theme ~ col_id,
  value.var = "Hyper_Adjp_BH_log10"
)

# Column annotations
col_split_labs <- ifelse(
  grepl("change_Female", colnames(p_matrix)),
  "Female-biased change",
  "Male-biased change"
)

col_names_labs <- sub(
  "^change_(Male|Female)_",
  "",
  colnames(p_matrix)
)

# Original colour scale
mat_col <- circlize::colorRamp2(
  c(5, 15, 40, 60, 125),
  c(
    "lightyellow",
    "gold",
    "orange",
    "orangered",
    "darkred"
  )
)

# ------------------------------------------------------------
# 5. Generate heatmap
# ------------------------------------------------------------

ht <- Heatmap(
  p_matrix,
  name = "-log10(Adj_P)",
  col = mat_col,
  na_col = "grey90",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  column_labels = col_names_labs,
  column_split = col_split_labs,
  column_names_rot = 45,
  row_names_max_width = unit(25, "cm"),
  column_title_gp = gpar(font = 2),
  border_gp = gpar(col = "black", lty = 1, lwd = 1.5),
  rect_gp = gpar(col = "black", lwd = 0.2),
  row_names_gp = gpar(fontsize = 14, fontface = "bold"),
  column_gap = unit(2, "mm"),
  column_names_gp = gpar(fontsize = 14, fontface = "bold"),
  row_title = NULL,
  show_column_names = TRUE
)

pdf(
  file.path(
    output_path,
    paste0(
      "rGREAT_biological_themes_modref_",
      mod_ref, ".pdf"
    )
  ),
  width = 12,
  height = 10
)

draw(ht)
dev.off()

cat("Biological-theme heatmap completed.\n")

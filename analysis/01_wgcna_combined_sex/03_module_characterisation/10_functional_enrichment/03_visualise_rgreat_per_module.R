
#!/usr/bin/env Rscript

# =============================================================================
# Visualise rGREAT GO enrichment across co-methylation modules
#
# Inputs:
#   - all_res_rGREAT_out.tsv
#   - res_top_go_terms.xlsx (manually curated)
#   - 01a_selected_mods_summary_stats.xlsx
#
# Output:
#   rGREAT_GO_enrichment_per_module.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# -----------------------------------------------------------------------------
# 1. Paths
# -----------------------------------------------------------------------------

input_dir <- "/results/01_wgcna_combined_sex"
output_dir <- input_dir

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

rgreat_file <- file.path(input_dir, "all_res_rGREAT_out.tsv")

module_file <- file.path(
  input_dir,
  "01a_selected_mods_summary_stats_wgcna.xlsx"
)

curated_file <- file.path(
  input_dir,
  "res_top_go_terms.xlsx"
)

# -----------------------------------------------------------------------------
# 2. Load module metadata
# -----------------------------------------------------------------------------

sel_mods <- as.data.table(read_xlsx(module_file))

sel_mods[, moduleNumID := paste0("module_", moduleNumID)]

sel_mods_id <- unique(
  sel_mods[, .(moduleColor, moduleNumID)]
)

# -----------------------------------------------------------------------------
# 3. Load consolidated rGREAT results
# -----------------------------------------------------------------------------

res_out <- fread(rgreat_file)

# Combined module results
res_out_sub_all <- res_out[
  moduleColor == "All"
]

# Select GO Molecular Function and Biological Process
res_sub <- res_out_sub_all[
  Group %in% c("MF", "BP"),
  .(
    ID,
    name,
    Hyper_Adjp_BH,
    Group,
    moduleColor,
    statsPres,
    sigGroup,
    Hyper_Adjp_BH_log10,
    analysis
  )
]

# -----------------------------------------------------------------------------
# 4. Reproduce original top-term selection
# -----------------------------------------------------------------------------

top_go_terms <- res_sub[
  ,
  .SD[which.max(Hyper_Adjp_BH_log10)],
  by = .(Group, moduleColor, statsPres, sigGroup)
]

# The original workflow exported top_go_terms for manual curation.
# The curated workbook is used below rather than overwriting it.

# -----------------------------------------------------------------------------
# 5. Load manually curated GO terms
# -----------------------------------------------------------------------------

defined_go <- as.data.table(
  read_xlsx(curated_file)
)

required_cols <- c(
  "moduleColor",
  "sigGroup",
  "Hyper_Adjp_BH_log10",
  "Biological_theme"
)

missing_cols <- setdiff(required_cols, names(defined_go))

if (length(missing_cols) > 0L) {
  stop(
    "Missing columns in curated workbook: ",
    paste(missing_cols, collapse = ", ")
  )
}

# Match module colours to numeric module identifiers
defined_go <- merge(
  defined_go,
  sel_mods_id,
  by = "moduleColor",
  all.x = TRUE
)

defined_go[
  ,
  col_id := paste0(moduleNumID, "_", sigGroup)
]

# Preserve original selection of one term per module
defined_go_sub <- defined_go[
  ,
  .SD[which.max(Hyper_Adjp_BH_log10)],
  by = moduleNumID
]

defined_go_sub <- defined_go_sub[
  order(-sigGroup, moduleNumID)
]

col_order <- defined_go_sub$col_id

# -----------------------------------------------------------------------------
# 6. Construct GO enrichment matrix
# -----------------------------------------------------------------------------

go_matrix <- dcast(
  defined_go_sub,
  Biological_theme ~ col_id,
  value.var = "Hyper_Adjp_BH_log10"
)

row_labels <- go_matrix$Biological_theme

go_matrix <- as.matrix(
  go_matrix[, ..col_order]
)

storage.mode(go_matrix) <- "numeric"

rownames(go_matrix) <- row_labels

# Module labels
col_names_labs <- sub(
  "_Stable|_Changing",
  "",
  colnames(go_matrix)
)

# Column groups derived from actual data
col_groups <- ifelse(
  grepl("_Stable$", colnames(go_matrix)),
  "Stable",
  "Changing"
)

col_split_labs <- factor(
  col_groups,
  levels = c("Stable", "Changing")
)

# -----------------------------------------------------------------------------
# 7. Heatmap colour scale
# -----------------------------------------------------------------------------

# Preserve original colour thresholds
mat_col <- colorRamp2(
  c(5, 20, 50, 100, 245),
  c(
    "lightyellow",
    "gold",
    "orange",
    "orangered",
    "darkred"
  )
)

# -----------------------------------------------------------------------------
# 8. Plot heatmap
# -----------------------------------------------------------------------------

ht <- Heatmap(
  go_matrix,

  name = "-log10(Adj_Pvalue)",
  col = mat_col,
  na_col = "grey",

  cluster_rows = FALSE,
  cluster_columns = FALSE,

  row_names_side = "left",

  column_labels = col_names_labs,
  column_split = col_split_labs,
  column_title = c("aSCMs", "aVCMs"),

  column_names_rot = 45,

  column_title_gp = gpar(fontface = "italic"),

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
  show_column_names = TRUE
)

# -----------------------------------------------------------------------------
# 9. Save figure
# -----------------------------------------------------------------------------

output_file <- file.path(
  output_dir,
  "rGREAT_GO_enrichment_per_module.pdf"
)

pdf(output_file, width = 18, height = 10)
draw(ht)
dev.off()

message("Heatmap saved: ", output_file)


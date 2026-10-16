#==============================================================================
# Visualise histone-mark enrichment in stable and changing WGCNA modules
# Input: 01a_histone_enrichment_res_combined_mods_n_het.tsv
# Output: heatmap of log(OR), 95% CIs, and nominal significance stars
#==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# Configure input/output paths
results_dir <- "/results/01_wgcna_combined_sex"
figure_dir <- "/figures"

# Plot heatmao for aSCMs/aVCMs or aSCPs/aVCPs
analysis <- "WGCNA"  # NB: Use "Heteroskedasticity" for aSCPs/aVCPs

res <- fread(file.path(results_dir, "01a_histone_enrichment_res_combined_mods_n_het.tsv"))
res <- res[Analysis == analysis]
if (nrow(res) == 0L) stop("No results found for analysis: ", analysis)

# Replace zero/non-positive OR or CI values in their own columns only.
for (column in c("OR", "CI_lower", "CI_upper")) {
  set(res, which(!is.na(res[[column]]) & res[[column]] <= 0), column, NA_real_)
  set(res, j = column, value = log(res[[column]]))
}

shorten_names <- function(x) {
  substitutions <- c(
    "A549.*Lung_Carcinoma" = "A549_lung_cancer",
    "HepG2.*Carcinoma" = "HepG2_liver_cancer",
    "Brain_Hippocampus_Middle" = "Brain_HPC",
    "Brain_Mid_Frontal_Lobe" = "Brain_Frontal",
    "Penis_Foreskin_Fibroblast_Primary_Cells_skin01" = "Skin_Fibroblast",
    "Peripheral_Blood_Mononuclear_Primary_Cells" = "PBMC",
    "CD34_Primary_Cells" = "CD34",
    "HMEC_Mammary_Epithelial" = "HMEC_mammary",
    "Skeletal_Muscle_Male" = "Muscle",
    "Adult_Liver" = "Liver",
    "Adipose_Nuclei" = "Adipose",
    "Left_Ventricle" = "Heart"
  )
  for (pattern in names(substitutions)) {
    x <- gsub(pattern, substitutions[[pattern]], x)
  }
  x
}

active_marks <- c("H3K27ac", "H3K4me3", "H3K9ac", "H3K4me1", "H3K36me3")
repressive_marks <- c("H3K27me3", "H3K9me3")
marks <- c(active_marks, repressive_marks)
groups <- c("Stable", "Changing")
ordered_rows <- unlist(lapply(groups, function(g) paste(marks, g, sep = "_")))
ordered_cols <- c(
  "PBMC", "CD34", "A549_lung_cancer", "HepG2_liver_cancer",
  "HMEC_mammary", "Skin_Fibroblast", "Adipose", "Liver",
  "Aorta", "Heart", "Brain_HPC", "Brain_Frontal", "Esophagus",
  "Lung", "Spleen", "Muscle"
)

res[, row_id := paste(histone_mark, sigGroup, sep = "_")]
res[, short_standard_epigenome := shorten_names(edacc_epigenome_name)]

# Require unique row/column combinations rather than silently aggregating.
if (res[, anyDuplicated(paste(row_id, short_standard_epigenome))] > 0L) {
  stop("Duplicate histone-mark/tissue combinations; check input metadata.")
}
missing_rows <- setdiff(ordered_rows, res$row_id)
missing_cols <- setdiff(ordered_cols, res$short_standard_epigenome)
if (length(missing_rows) || length(missing_cols)) {
  stop("Missing expected rows: ", paste(missing_rows, collapse = ", "),
       "; missing tissues: ", paste(missing_cols, collapse = ", "))
}

# Match historical display: round log odds ratios before plotting.
res[, OR := round(OR, 2)]
res[, OR_CI_label := fifelse(
  is.na(OR) | is.na(CI_lower) | is.na(CI_upper),
  NA_character_,
  paste0(round(OR, 2), "\n(", round(CI_lower, 2), ",",
         round(CI_upper, 2), ")")
)]
res[, p_stars := fifelse(p_value < 0.001, "***",
                        fifelse(p_value < 0.01, "**",
                                fifelse(p_value < 0.05, "*", NA_character_)))]

make_matrix <- function(column, fill = NA) {
  m <- data.table::dcast(res, row_id ~ short_standard_epigenome,
                         value.var = column, fill = fill)
  rownames(m) <- m$row_id
  m$row_id <- NULL
  as.matrix(m)[ordered_rows, ordered_cols, drop = FALSE]
}

or_matrix <- make_matrix("OR", NA_real_)
ci_matrix <- make_matrix("OR_CI_label", NA_character_)
star_matrix <- make_matrix("p_stars", NA_character_)

finite_or <- or_matrix[is.finite(or_matrix)]
if (!length(finite_or)) stop("No finite odds ratios to plot.")
log_min <- min(floor(min(finite_or)), -0.01)
log_max <- max(ceiling(max(finite_or)), 0.01)
colour_fun <- circlize::colorRamp2(
  c(log_min, 0, log_max), c("skyblue", "white", "red")
)

row_labels <- sub("_.*", "", ordered_rows)
row_colours <- ifelse(row_labels %in% active_marks, "forestgreen", "purple")
row_split <- factor(rep(groups, each = length(marks)), levels = groups)
block_labels <- if (analysis == "WGCNA") c("aSCMs", "aVCMs") else c("aSCPs", "aVCPs")

ht <- ComplexHeatmap::Heatmap(
  or_matrix,
  name = "log(OR)",
  col = colour_fun,
  na_col = "grey90",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  row_labels = row_labels,
  column_names_rot = 45,
  border_gp = grid::gpar(col = "black", lty = 1, lwd = 1.5),
  rect_gp = grid::gpar(col = "black", lwd = 0.2),
  row_names_gp = grid::gpar(fontsize = 14, fontface = "bold", col = row_colours),
  column_names_gp = grid::gpar(fontsize = 14, fontface = "bold"),
  row_gap = grid::unit(1, "mm"),
  column_gap = grid::unit(2, "mm"),
  row_split = row_split,
  row_title = NULL,
  show_column_names = TRUE,
  left_annotation = ComplexHeatmap::rowAnnotation(
    cpgAnnot = ComplexHeatmap::anno_block(
      gp = grid::gpar(fill = c("white", "white")),
      labels = block_labels,
      labels_gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(ci_matrix[i, j])) {
      grid::grid.text(ci_matrix[i, j], x, y + grid::unit(2.5, "mm"),
                      gp = grid::gpar(fontsize = 12))
    }
    if (!is.na(star_matrix[i, j]) && nzchar(star_matrix[i, j])) {
      grid::grid.text(star_matrix[i, j], x, y - grid::unit(6, "mm"),
                      gp = grid::gpar(fontsize = 10))
    }
  }
)

# Save figure. Large dimensions accommodate text in the original cell layout.
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
pdf(file.path(figure_dir, "histone_mark_enrichment_heatmap.pdf"),
    width = 22, height = 13, useDingbats = FALSE)
ComplexHeatmap::draw(ht)
dev.off()


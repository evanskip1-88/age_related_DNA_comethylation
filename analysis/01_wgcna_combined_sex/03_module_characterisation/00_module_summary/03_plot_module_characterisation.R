
# =============================================================================
# Purpose: Plot module characteristics for aSCMs and aVCMs (Figure 2c).
# Inputs:  Selected modules, annotations, covariates and raw DNAm beta values.
# Outputs: Figure 2c.
# Run from the repository root.
# =============================================================================

library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(wesanderson)
library(grid)

# -------------------------------------------------------------------------
# Input and output paths
# -------------------------------------------------------------------------

input_dir <- "results/01_wgcna_combined_sex/"
beta_dir <- "..." # Individual level restricted covariate data
heteroskid_dir <- "results/02_wgcna_within_sex"

output_file <- "figures/Figure_2c_module_characterisation.pdf"

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Load selected modules and annotations
# -------------------------------------------------------------------------

dt_WGCNA_sel_mods <- as.data.table(
  readxl::read_xlsx(
    file.path(input_dir, "01a_selected_mods_summary_stats_wgcna.xlsx")
  )
)

# Number of unique genes per module
gene_counts <- dt_WGCNA_sel_mods[
  ,
  .(
    num_genes = length(
      unique(na.omit(unlist(strsplit(gene_rgreat, ";"))))
    )
  ),
  by = moduleNumID
]

dt_WGCNA_sel_mods <- merge(
  dt_WGCNA_sel_mods,
  gene_counts,
  by = "moduleNumID"
)

dt_WGCNA_sel_mods[
  ,
  moduleNumID := paste0("module_", moduleNumID)
]

genic_annotr <- fread(
  file.path(input_dir, "GS_sub_qc1_h19_genic.annotr")
)

cpgs_annotr <- fread(
  file.path(input_dir, "GS_sub_qc1_h19_cpgs.annotr")
)

# -------------------------------------------------------------------------
# Prepare module annotation data
# -------------------------------------------------------------------------

dt_WGCNA_sel_mods_m <- merge(
  genic_annotr[, .(probe_id, annot.type)],
  dt_WGCNA_sel_mods,
  by = "probe_id"
)

dt_WGCNA_sel_mods_m <- merge(
  cpgs_annotr[, .(probe_id, annot.type)],
  dt_WGCNA_sel_mods_m,
  by = "probe_id"
)

dt_WGCNA_sel_mods_m_sub <- unique(
  dt_WGCNA_sel_mods_m[
    Category != "Others",
    .(
      probe_id,
      annot_cpg = annot.type.x,
      annot_genic = annot.type.y,
      moduleNumID,
      Category,
      change_presStats,
      moduleSize,
      dist_to_TSS,
      gene_rgreat,
      num_genes
    )
  ]
)

# Minimum absolute distance to TSS
dt_WGCNA_sel_mods_m_sub[
  ,
  abs_distTSS := vapply(
    strsplit(as.character(dist_to_TSS), ";"),
    function(x) min(abs(as.numeric(x))),
    numeric(1)
  )
]

dt_WGCNA_sel_mods_m_sub[
  ,
  dist_TSS_grp := fcase(
    abs_distTSS / 1000 < 1, "<1kb",
    abs_distTSS / 1000 > 1 &
      abs_distTSS / 1000 <= 5, "1kb~5kb",
    abs_distTSS / 1000 > 5 &
      abs_distTSS / 1000 <= 10, "5kb~10kb",
    abs_distTSS / 1000 > 10, ">10kb",
    default = NA_character_
  )
]

# -------------------------------------------------------------------------
# Summarise annotation percentages
# -------------------------------------------------------------------------

annotation_summary <- function(dt, annot_col) {

  dt_modsize <- unique(
    dt[, .(moduleNumID, Category, moduleSize)]
  )

  dt_total <- dt[
    ,
    .(Total_Probes = .N),
    by = .(moduleNumID, Category)
  ]

  dt_annot_counts <- dt[
    ,
    .N,
    by = .(
      moduleNumID,
      Category,
      annotation = get(annot_col)
    )
  ]

  dt_annot_counts <- merge(
    dt_annot_counts,
    dt_total,
    by = c("moduleNumID", "Category")
  )

  dt_annot_counts[
    ,
    Proportion := (N / Total_Probes) * 100
  ]

  dt_wide <- dcast(
    dt_annot_counts,
    moduleNumID + Category + Total_Probes ~ annotation,
    value.var = "Proportion",
    fill = 0
  )

  dt_wide <- merge(
    dt_wide,
    dt_modsize,
    by = c("moduleNumID", "Category")
  )

  dt_wide[, Total_Probes := NULL]

  dt_wide
}

dt_cpg <- annotation_summary(
  dt_WGCNA_sel_mods_m_sub,
  "annot_cpg"
)

setnames(dt_cpg, sub("hg19_", "", names(dt_cpg)))

dt_genic <- annotation_summary(
  dt_WGCNA_sel_mods_m_sub,
  "annot_genic"
)

setnames(dt_genic, sub("hg19_genes_", "", names(dt_genic)))

# -------------------------------------------------------------------------
# Prepare CpG and genic annotation matrices
# -------------------------------------------------------------------------

cpg_cols <- c(
  "cpg_islands",
  "cpg_shores",
  "cpg_shelves",
  "cpg_inter"
)

genic_cols <- c(
  "1to5kb",
  "promoters",
  "5UTRs",
  "exons",
  "introns",
  "3UTRs",
  "intergenic"
)

dt_cpg_order <- dt_cpg[
  order(-Category, moduleNumID)
]

cpg_mat <- as.matrix(
  dt_cpg_order[, ..cpg_cols]
)

rownames(cpg_mat) <- dt_cpg_order$moduleNumID

dt_genic_order <- dt_genic[
  order(-Category, moduleNumID)
]

genic_mat <- as.matrix(
  dt_genic_order[, ..genic_cols]
)

rownames(genic_mat) <- dt_genic_order$moduleNumID

# -------------------------------------------------------------------------
# Module CpG counts, gene counts and preservation changes
# -------------------------------------------------------------------------

dt_counts_summary <- unique(
  dt_WGCNA_sel_mods_m_sub[
    ,
    .(
      moduleNumID,
      Category,
      moduleSize,
      num_genes,
      change_presStats
    )
  ]
)

# Preserve original selection: maximum change per module
dt_counts_summary <- dt_counts_summary[
  ,
  .SD[which.max(change_presStats)],
  by = moduleNumID
]

dt_counts_summary_order <- dt_counts_summary[
  order(-Category, moduleNumID)
]

dt_counts_summary_mat <- as.matrix(
  dt_counts_summary_order[
    ,
    .(moduleSize, num_genes, change_presStats)
  ]
)

rownames(dt_counts_summary_mat) <-
  dt_counts_summary_order$moduleNumID

colnames(dt_counts_summary_mat) <- c(
  "num_cpgs",
  "num_genes",
  "ChangePresStats"
)

# -------------------------------------------------------------------------
# Distance-to-TSS annotation matrix
# -------------------------------------------------------------------------

dt_dist_TSS_grp_counts <- dt_WGCNA_sel_mods_m_sub[
  ,
  .N,
  by = .(moduleNumID, Category, dist_TSS_grp)
]

dt_total_probes <- dt_WGCNA_sel_mods_m_sub[
  ,
  .(Total_Probes = .N),
  by = .(moduleNumID, Category)
]

dt_dist_TSS_grp_counts <- merge(
  dt_dist_TSS_grp_counts,
  dt_total_probes,
  by = c("moduleNumID", "Category")
)

dt_dist_TSS_grp_counts[
  ,
  Proportion := (N / Total_Probes) * 100
]

dt_dist_TSS_grp <- dcast(
  dt_dist_TSS_grp_counts,
  moduleNumID + Category + Total_Probes ~ dist_TSS_grp,
  value.var = "Proportion",
  fill = 0
)

tss_cols <- c(
  "<1kb",
  "1kb~5kb",
  "5kb~10kb",
  ">10kb"
)

dt_dist_TSS_grp_order <- dt_dist_TSS_grp[
  order(-Category, moduleNumID)
]

dist_TSS_mat <- as.matrix(
  dt_dist_TSS_grp_order[, ..tss_cols]
)

rownames(dist_TSS_mat) <-
  dt_dist_TSS_grp_order$moduleNumID

# -------------------------------------------------------------------------
# Mean DNA methylation across age groups
# -------------------------------------------------------------------------

dt_cov <- fread(
  file.path(beta_dir, "GS_18869_qcov.cov"),
  fill = TRUE
)

dt_cov[
  ,
  age_grp := fcase(
    age >= 17.9 & age < 36.0, "Age1",
    age >= 36.0 & age < 44.0, "Age2",
    age >= 44.0 & age < 49.5, "Age3",
    age >= 49.5 & age < 54.2, "Age4",
    age >= 54.2 & age < 57.8, "Age5",
    age >= 57.8 & age < 61.5, "Age6",
    age >= 61.5 & age < 66.0, "Age7",
    age >= 66.0 & age < 94.0, "Age8",
    default = "NA"
  )
]

calculate_mean_dnam <- function(category) {

  beta_file <- file.path(
    beta_dir,
    paste0("res_raw_beta_WGCNA_", category, ".txt")
  )

  dt_beta <- fread(beta_file)

  dt_beta <- merge(
    dt_cov[, .(IID_DNAm, age, age_grp)],
    dt_beta[, -1, with = FALSE],
    by.x = "IID_DNAm",
    by.y = "IID"
  )

  dt_beta <- dt_beta[order(age)]

  dt_long <- melt(
    dt_beta,
    id.vars = c("IID_DNAm", "age", "age_grp"),
    variable.name = "probe_id",
    value.name = "DNAm"
  )

  dt_long <- merge(
    unique(
      dt_WGCNA_sel_mods[
        Category == category,
        .(probe_id, moduleNumID)
      ]
    ),
    dt_long,
    by = "probe_id"
  )

  dt_means <- dt_long[
    ,
    .(Mean_DNAm = mean(DNAm, na.rm = TRUE)),
    by = .(moduleNumID, age_grp)
  ]

  dt_means_wide <- dcast(
    dt_means,
    moduleNumID ~ age_grp,
    value.var = "Mean_DNAm"
  )

  dt_means_wide[
    ,
    Category := category
  ]

  dt_means_wide
}

dt_means_c <- rbindlist(
  list(
    calculate_mean_dnam("Stable"),
    calculate_mean_dnam("Changing")
  ),
  use.names = TRUE,
  fill = TRUE
)

dt_means_c_order <- dt_means_c[
  order(-Category, moduleNumID)
]

age_cols <- paste0("Age", 1:8)

means_c_mat <- as.matrix(
  dt_means_c_order[, ..age_cols]
)

rownames(means_c_mat) <-
  dt_means_c_order$moduleNumID

# -------------------------------------------------------------------------
# Verify matrix alignment
# -------------------------------------------------------------------------

module_order <- dt_counts_summary_order$moduleNumID

stopifnot(
  setequal(module_order, rownames(cpg_mat)),
  setequal(module_order, rownames(genic_mat)),
  setequal(module_order, rownames(dist_TSS_mat)),
  setequal(module_order, rownames(means_c_mat))
)

cpg_mat <- cpg_mat[module_order, , drop = FALSE]
genic_mat <- genic_mat[module_order, , drop = FALSE]
dist_TSS_mat <- dist_TSS_mat[module_order, , drop = FALSE]
means_c_mat <- means_c_mat[module_order, , drop = FALSE]

group_order <- dt_counts_summary_order$Category

probe_counts_mat <- dt_counts_summary_mat[, 1, drop = FALSE]
gene_counts_mat <- dt_counts_summary_mat[, 2, drop = FALSE]
presstat_change_mat <- dt_counts_summary_mat[, 3, drop = FALSE]

# -------------------------------------------------------------------------
# Figure 2c: Colours
# -------------------------------------------------------------------------

cpgs_col <- c(
  "#008000",
  "#C2E105",
  "#8A91D0",
  "#CD5C5C"
)

genic_col <- wes_palette(
  "Zissou1",
  10,
  type = "continuous"
)

tss_col <- wes_palette(
  "Moonrise3",
  10,
  type = "continuous"
)

# -------------------------------------------------------------------------
# Figure 2c: Module labels and annotation barplots
# -------------------------------------------------------------------------

mod_labels <- sub("module_", "", module_order)

ht_label <- rowAnnotation(
  text = anno_text(
    mod_labels,
    location = unit(1, "npc"),
    just = "right",
    gp = gpar(fontsize = 10)
  )
)

ht_probe_counts <- rowAnnotation(
  probe_counts = anno_barplot(
    probe_counts_mat,
    bar_width = 0.5,
    gp = gpar(fill = "grey30"),
    width = unit(ncol(dt_counts_summary_mat) * 7, "mm"),
    axis_param = list(gp = gpar(fontsize = 10))
  ),
  show_annotation_name = FALSE
)

ht_gene_counts <- rowAnnotation(
  gene_counts = anno_barplot(
    gene_counts_mat,
    bar_width = 0.5,
    gp = gpar(fill = "grey"),
    width = unit(ncol(dt_counts_summary_mat) * 7, "mm"),
    axis_param = list(gp = gpar(fontsize = 10))
  ),
  show_annotation_name = FALSE
)

ht_absdiffChange <- rowAnnotation(
  absdiffChange = anno_barplot(
    presstat_change_mat,
    bar_width = 1,
    gp = gpar(fill = "orange"),
    width = unit(ncol(dt_counts_summary_mat) * 7, "mm"),
    axis_param = list(gp = gpar(fontsize = 10))
  ),
  show_annotation_name = FALSE
)

ht_cpgs_annot <- rowAnnotation(
  cpgs_annot = anno_barplot(
    cpg_mat,
    bar_width = 1,
    gp = gpar(fill = cpgs_col),
    width = unit(ncol(cpg_mat) * 8, "mm"),
    axis_param = list(gp = gpar(fontsize = 10))
  ),
  show_annotation_name = FALSE
)

ht_genic_annot <- rowAnnotation(
  genic_annot = anno_barplot(
    genic_mat,
    bar_width = 1,
    gp = gpar(fill = genic_col),
    width = unit(ncol(genic_mat) * 4, "mm"),
    axis_param = list(gp = gpar(fontsize = 10))
  ),
  show_annotation_name = FALSE
)

ht_tss <- rowAnnotation(
  dist_tss = anno_barplot(
    dist_TSS_mat,
    bar_width = 1,
    gp = gpar(fill = tss_col),
    width = unit(ncol(dist_TSS_mat) * 6, "mm"),
    axis_param = list(gp = gpar(fontsize = 10))
  ),
  show_annotation_name = FALSE
)

# -------------------------------------------------------------------------
# Invisible dummy heatmap for group labels and row separation
# -------------------------------------------------------------------------

dummy_mat <- matrix(
  0,
  nrow = length(group_order),
  ncol = 1,
  dimnames = list(module_order, "dummy")
)

ht_split <- Heatmap(
  dummy_mat,
  name = "dummy",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  width = unit(1, "mm"),
  col = c("0" = "white"),
  border = FALSE,
  rect_gp = gpar(col = NA, fill = NA),
  row_split = factor(
    group_order,
    levels = c("Stable", "Changing")
  ),
  row_gap = unit(2, "mm"),
  row_title = NULL,
  show_heatmap_legend = FALSE,
  left_annotation = rowAnnotation(
    foo = anno_block(
      gp = gpar(
        alpha = 0.5,
        fill = c("gray100", "gray100")
      ),
      labels = c("aSCMs", "aVCMs"),
      labels_gp = gpar(
        fontsize = 12,
        fontface = "bold"
      )
    ),
    width = unit(6, "mm")
  )
)

# -------------------------------------------------------------------------
# Combine Figure 2c panels
# -------------------------------------------------------------------------

ht_list <-
  ht_label +
  ht_probe_counts +
  ht_gene_counts +
  ht_absdiffChange +
  ht_cpgs_annot +
  ht_genic_annot +
  ht_tss +
  ht_split

# -------------------------------------------------------------------------
# Legends
# -------------------------------------------------------------------------

legend_title_gp <- gpar(
  fontsize = 12,
  fontface = "bold"
)

legend_labels_gp <- gpar(fontsize = 12)

cpg_lgd <- Legend(
  labels = c(
    "cpg_islands",
    "cpg_shores",
    "cpg_shelves",
    "cpg_inter"
  ),
  title = "CpG annotation",
  title_gp = legend_title_gp,
  labels_gp = legend_labels_gp,
  legend_gp = gpar(fill = cpgs_col)
)

genic_lgd <- Legend(
  labels = c(
    "1to5kb",
    "promoters",
    "5UTRs",
    "exons",
    "introns",
    "intergenic",
    "3UTRs"
  ),
  title = "Gene annotation",
  title_gp = legend_title_gp,
  labels_gp = legend_labels_gp,
  legend_gp = gpar(fill = genic_col)
)

tss_lgd <- Legend(
  labels = c(
    "<1kb",
    "1kb~5kb",
    "5kb~10kb",
    ">10kb"
  ),
  title = "Distance to TSS",
  title_gp = legend_title_gp,
  labels_gp = legend_labels_gp,
  legend_gp = gpar(fill = tss_col)
)

# -------------------------------------------------------------------------
# Save Figure 2c
# -------------------------------------------------------------------------

pdf(output_file, width = 16, height = 10)

draw(
  ht_list,
  padding = unit(c(3, 3, 10, 3), "mm"),
  heatmap_legend_list = list(
    cpg_lgd,
    genic_lgd,
    tss_lgd
  ),
  heatmap_legend_side = "bottom"
)

# Annotation headings
anno_title <- c(
  probe_counts = "Number of \nCpGs",
  gene_counts = "Number of \ngenes",
  absdiffChange = "Change in \npres. stats",
  cpgs_annot = "CpG \nannotation",
  genic_annot = "Gene \nannotation",
  dist_tss = "Distance to \nTSS"
)

for (an in names(anno_title)) {
  decorate_annotation(an, {
    grid.text(
      anno_title[an],
      y = unit(1.1, "npc") + unit(3, "mm"),
      just = "top",
      gp = gpar(
        fontsize = 10,
        fontface = "bold"
      )
    )
  })
}

grid.text(
  "Module ID",
  x = unit(-1, "mm"),
  y = unit(0.5, "npc"),
  rot = 90,
  gp = gpar(
    fontsize = 12,
    fontface = "bold"
  )
)

dev.off()


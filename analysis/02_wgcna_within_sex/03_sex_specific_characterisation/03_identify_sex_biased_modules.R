# ============================================================
# Identify sex-biased changes in WGCNA module preservation
#
# Reference network: Male or female Age1
# Main comparison: Age2 vs Age8
# Additional comparison: Age2 vs Age4
#
# Sex bias is defined using the difference in absolute
# preservation change between males and females.
#
# Outputs:
#   1. Sex-biased module summary
#   2. GREAT-associated gene annotations
#   3. Male-female preservation-change scatterplot
#   4. Age2–Age4 exploratory comparison
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
})

# ------------------------------------------------------------
# 1. Configure directories and parameters
# ------------------------------------------------------------

data_path <- "results/02_wgcna_within_sex/"
opi_path <- "data/gsDNAm/"
output_path <- "..."

mod_ref <- "M"  # Reference network: "M" or "F"

sexDifferenceCutoff <- 0.05
minGeneCpGs <- 5

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Load paired preservation statistics from Part B
# ------------------------------------------------------------

dt_paired <- fread(
  file.path(
    data_path,
    paste0(
      "paired_preservation_by_sex_modref_",
      mod_ref, ".tsv"
    )
  )
)

# Load reference module labels
dt_mod_labs <- fread(
  file.path(
    data_path,
    paste0(
      "precorr_gs_qc_chr_all_refage1_MinMod20_",
      mod_ref, ".modLabs"
    )
  )
)

# Load GREAT CpG-gene annotations
opi_rgreat <- fread(
  file.path(
    opi_path,
    "GS_sub_qc1_gene_cpgs_assoc_rGREAT_h19_all_CpGs.opi"
  )
)

setnames(
  opi_rgreat,
  old = names(opi_rgreat)[c(2, 4, 7)],
  new = c("chr", "gene_rgreat", "orig_gene")
)

opi_rgreat <- opi_rgreat[position != "NA"]

# ------------------------------------------------------------
# 3. Identify sex-biased modules: Age2 vs Age8
# ------------------------------------------------------------

# Retain one row per module at Age8
# Part A already calculated absolute Age2–Age8 changes
dt_sub <- copy(
  dt_paired[testAge == "Age8" & presStats == "corkME"]
)

# Difference in magnitude of preservation change
dt_sub[
  , abs_diff_change := abs(
    change_presStats_M - change_presStats_F
  )
]

# Identify the sex with the larger absolute change
dt_sub[
  , sexChangeGroup := fifelse(
    abs_diff_change < sexDifferenceCutoff,
    "Others",
    fifelse(
      change_presStats_M > change_presStats_F,
      "change_Male",
      "change_Female"
    )
  )
]

# Check whether either sex was classified as age-variable
dt_sub[
  , consistChangeTrend := fifelse(
    sigGroup_M == "aVCMs" | sigGroup_F == "aVCMs",
    "yes",
    "no"
  )
]

# Rank modules by sex difference
setorder(dt_sub, -abs_diff_change)

# ------------------------------------------------------------
# 4. Attach module IDs
# ------------------------------------------------------------

# One row per reference module
dt_mod_map <- unique(
  dt_mod_labs[
    , .(moduleNumID, moduleColors)
  ]
)

dt_annotated <- merge(
  dt_sub,
  dt_mod_map,
  by.x = "moduleColor",
  by.y = "moduleColors",
  all.x = TRUE
)

# ------------------------------------------------------------
# 5. Annotate sex-biased modules with GREAT genes
# ------------------------------------------------------------

# Restrict gene annotation to sex-biased modules
dt_sex_biased <- dt_annotated[
  sexChangeGroup != "Others"
]

# Match CpGs to reference modules
dt_mod_cpgs <- merge(
  dt_mod_labs,
  dt_sex_biased[
    , .(moduleColor, moduleNumID)
  ],
  by.x = c("moduleColors", "moduleNumID"),
  by.y = c("moduleColor", "moduleNumID")
)

# Attach GREAT gene associations
dt_mod_genes <- merge(
  dt_mod_cpgs,
  opi_rgreat[
    , .(probe_id, gene_rgreat)
  ],
  by = "probe_id"
)

# Split multiple gene assignments
dt_gene_counts <- dt_mod_genes[
  !is.na(gene_rgreat) & gene_rgreat != "",
  .(
    gene = unlist(
      strsplit(gene_rgreat, ";", fixed = TRUE)
    )
  ),
  by = .(moduleColors, moduleNumID, probe_id)
][
  ,
  .(count = uniqueN(probe_id)),
  by = .(moduleColors, moduleNumID, gene)
]

# Retain genes associated with more than five CpGs
dt_gene_counts <- dt_gene_counts[
  count > minGeneCpGs
]

setorder(
  dt_gene_counts,
  moduleColors,
  -count,
  gene
)

# Select up to two representative genes per module
top2_genes_per_module <- dt_gene_counts[
  ,
  .(
    geneCluster = paste(
      head(gene, 2),
      collapse = "/"
    )
  ),
  by = .(moduleColors, moduleNumID)
]

# Add gene labels to module summary
dt_annotated <- merge(
  dt_annotated,
  top2_genes_per_module,
  by.x = c("moduleColor", "moduleNumID"),
  by.y = c("moduleColors", "moduleNumID"),
  all.x = TRUE
)

setorder(dt_annotated, -abs_diff_change)

# ------------------------------------------------------------
# 6. Plot sex-biased changes in preservation
# ------------------------------------------------------------

group_colors <- c(
  "change_Male" = "red",
  "change_Female" = "blue"
)

if (mod_ref == "M") {

  x_labs <- "Δ preservation (males - reference network)"
  y_labs <- "Δ preservation (females)"

} else {

  x_labs <- "Δ preservation (males)"
  y_labs <- "Δ preservation (females - reference network)"
}

p_sex_bias <- ggplot(
  dt_annotated,
  aes(
    x = change_presStats_M,
    y = change_presStats_F
  )
) +
  geom_text_repel(
    data = dt_annotated[sexChangeGroup == "Others"],
    aes(label = moduleNumID),
    size = 4,
    color = "grey40",
    segment.color = NA,
    max.overlaps = Inf
  ) +
  geom_text_repel(
    data = dt_annotated[sexChangeGroup != "Others"],
    aes(label = moduleNumID),
    size = 4,
    fontface = "bold",
    color = "grey20",
    segment.color = NA,
    max.overlaps = Inf
  ) +
  geom_text_repel(
    data = dt_annotated[
      sexChangeGroup != "Others" & !is.na(geneCluster)
    ],
    aes(
      label = geneCluster,
      color = sexChangeGroup
    ),
    size = 4,
    fontface = "bold",
    segment.color = "grey40",
    max.overlaps = Inf
  ) +
  stat_cor(
    aes(label = after_stat(r.label)),
    method = "pearson",
    label.x.npc = "left",
    label.y.npc = "bottom",
    size = 6
  ) +
  scale_color_manual(
    values = group_colors,
    name = NULL
  ) +
  labs(
    x = x_labs,
    y = y_labs
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  )

print(p_sex_bias)

# ------------------------------------------------------------
# 7. Explore midlife preservation changes: Age2 vs Age4
# ------------------------------------------------------------

dt_age2 <- dt_paired[
  testAge == "Age2" & presStats == "corkME",
  .(
    moduleColor,
    pres_Age2_M = pres_obs_value_M,
    pres_Age2_F = pres_obs_value_F
  )
]

dt_age4 <- dt_paired[
  testAge == "Age4" & presStats == "corkME",
  .(
    moduleColor,
    pres_Age4_M = pres_obs_value_M,
    pres_Age4_F = pres_obs_value_F
  )
]

changeAge2vs4 <- merge(
  dt_age2,
  dt_age4,
  by = "moduleColor"
)

# Absolute preservation change in each sex
changeAge2vs4[
  , change_presStats_Age2vs4_M :=
    abs(pres_Age4_M - pres_Age2_M)
]

changeAge2vs4[
  , change_presStats_Age2vs4_F :=
    abs(pres_Age4_F - pres_Age2_F)
]

# Fix missing calculation in original script
changeAge2vs4[
  , abs_diff_change2v4 := abs(
    change_presStats_Age2vs4_M -
      change_presStats_Age2vs4_F
  )
]

changeAge2vs4[
  , sexChangeGroup2v4 := fifelse(
    abs_diff_change2v4 < sexDifferenceCutoff,
    "Others",
    fifelse(
      change_presStats_Age2vs4_M >
        change_presStats_Age2vs4_F,
      "change_Male",
      "change_Female"
    )
  )
]

# Attach Age2–Age4 comparison to module summary
dt_annotated <- merge(
  dt_annotated,
  changeAge2vs4[
    , .(
      moduleColor,
      change_presStats_Age2vs4_M,
      change_presStats_Age2vs4_F,
      abs_diff_change2v4,
      sexChangeGroup2v4
    )
  ],
  by = "moduleColor",
  all.x = TRUE
)

setorder(dt_annotated, -abs_diff_change)

# ------------------------------------------------------------
# 8. Save results
# ------------------------------------------------------------

fwrite(
  dt_annotated,
  file.path(
    output_path,
    paste0(
      "01a_sex_specific_mods_modref_F",
      mod_ref, ".tsv"
    )
  ),
  sep = "\t"
)

fwrite(
  dt_gene_counts,
  file.path(
    output_path,
    paste0(
      "01a_sex_specific_mods_modref_F",
      mod_ref, ".tsv"
    )
  ),
  sep = "\t"
)

fwrite(
  changeAge2vs4,
  file.path(
    output_path,
    paste0(
      "01a_sex_specific_mods_modref_M_age2v4",
      mod_ref, ".tsv"
    )
  ),
  sep = "\t"
)

ggsave(
  filename = file.path(
    output_path,
    paste0(
      "sex_biased_preservation_changes_modref_",
      mod_ref, ".png"
    )
  ),
  plot = p_sex_bias,
  width = 10,
  height = 8,
  dpi = 300
)

cat("Sex-biased module analysis completed.\n")

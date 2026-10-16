
# =============================================================================
# Purpose: Plot the WGCNA module eigengene dendrogram (Figure 2a).
# Inputs:  Age1 module eigengenes, module labels, selected module summary.
# Outputs: Figure 2a module eigengene dendrogram.
# =============================================================================

library(data.table)
library(WGCNA)
library(ggplot2)
library(readxl)

# -------------------------------------------------------------------------
# Input and output paths
# -------------------------------------------------------------------------

module_summary_file <- "results/01_wgcna_combined_sex/01a_selected_mods_summary_stats.xlsx"

module_labels_file <- "results/01_wgcna_combined_sex/precorr_gs_qc_chr_all_refage1_wg.modLabs"

eigengene_file <- "..." # Individual level data identifiers

figure_file <- "figures/Figure_2a_module_eigengene_dendrogram.pdf"

dir.create(dirname(figure_file), recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Load selected age-stable and age-variable modules
# -------------------------------------------------------------------------

dt_sel_mods <- as.data.table(read_xlsx(module_summary_file))

sel_mods_uniq <- unique(
  dt_sel_mods[
    Category %in% c("Stable", "Changing"),
    .(moduleNumID, Category)
  ]
)

sel_mods_uniq[, moduleNumID := as.character(moduleNumID)]

aSCMs <- sel_mods_uniq[
  Category == "Stable",
  moduleNumID
]

aVCMs <- sel_mods_uniq[
  Category == "Changing",
  moduleNumID
]

# -------------------------------------------------------------------------
# Load module labels and eigengenes
# -------------------------------------------------------------------------

dt_mod_labs <- fread(module_labels_file)

dt_mod_labs[
  , ME_mods_labs := paste0("ME", moduleColors)
]

dt_mod_labs[
  , moduleNumID := as.character(moduleNumID)
]

dt_uniq_mods <- unique(
  dt_mod_labs[, .(ME_mods_labs, moduleNumID)]
)

MEs <- fread(eigengene_file)

MEs_mat <- as.matrix(MEs[, -1, with = FALSE])
rownames(MEs_mat) <- MEs[[1]]

# Map eigengene colour labels to numeric module IDs
name_map <- setNames(
  dt_uniq_mods$moduleNumID,
  dt_uniq_mods$ME_mods_labs
)

mapped_names <- name_map[colnames(MEs_mat)]

if (anyNA(mapped_names)) {
  stop("Some eigengene columns could not be mapped to module IDs.")
}

colnames(MEs_mat) <- unname(mapped_names)

MEs_df <- as.data.frame(MEs_mat)

# -------------------------------------------------------------------------
# Construct module eigengene dendrogram
# -------------------------------------------------------------------------

ME_cor <- cor(MEs_df, use = "p")

ME_diss <- 1 - ME_cor

ME_tree <- hclust(
  as.dist(ME_diss),
  method = "average"
)

# -------------------------------------------------------------------------
# Classify modules for annotation
# -------------------------------------------------------------------------

mod_class <- ifelse(
  colnames(MEs_df) %in% aSCMs,
  "Stable",
  ifelse(
    colnames(MEs_df) %in% aVCMs,
    "Changing",
    "Other"
  )
)

mod_colors <- ifelse(
  mod_class == "Stable",
  "blue",
  ifelse(
    mod_class == "Changing",
    "red",
    "grey"
  )
)

# -------------------------------------------------------------------------
# Plot Figure 2a
# -------------------------------------------------------------------------

pdf(
  figure_file,
  width = 10,
  height = 7
)

plotDendroAndColors(
  ME_tree,
  colors = mod_colors,
  groupLabels = "Module class",
  main = ""
)

dev.off()

cat("Figure 2a saved to:", figure_file, "\n")

# -------------------------------------------------------------------------
# Figure 2b: Module preservation trajectories
# -------------------------------------------------------------------------

# Load module preservation summary
preservation_file <- "results/01_wgcna_combined_sex/01a_selected_mods_summary_stats_wgcna.xlsx"

res <- as.data.table(read_xlsx(preservation_file))

# Explicit group order
res[, sigGroup := factor(
  sigGroup,
  levels = c("Others", "aSCMs", "aVCMs")
)]

mod_colors_grps <- c(
  "aSCMs" = "blue",
  "aVCMs" = "red"
)

p <- ggplot() +
  geom_line(
    data = res[sigGroup == "Others"],
    aes(
      x = testAgeNum,
      y = pres_obs_value,
      group = moduleColor
    ),
    color = "grey",
    linewidth = 1.1,
    alpha = 1.0
  ) +
  geom_point(
    data = res[sigGroup == "Others"],
    aes(
      x = testAgeNum,
      y = pres_obs_value,
      group = moduleColor
    ),
    color = "grey",
    size = 1.4,
    alpha = 1.0
  ) +
  geom_line(
    data = res[sigGroup != "Others"],
    aes(
      x = testAgeNum,
      y = pres_obs_value,
      color = sigGroup,
      group = moduleColor
    ),
    linewidth = 1.1
  ) +
  geom_point(
    data = res[sigGroup != "Others"],
    aes(
      x = testAgeNum,
      y = pres_obs_value,
      color = sigGroup,
      group = moduleColor
    ),
    size = 2
  ) +
  ylim(0, 1.0) +
  scale_color_manual(
    values = mod_colors_grps,
    breaks = c("aSCMs", "aVCMs"),
    name = NULL
  ) +
  facet_wrap(~ presStats, scales = "free_y", nrow = 1) +
  labs(
    x = "Test age group",
    y = "Module preservation"
  ) +
  theme_bw(base_size = 16, base_family = "Arial") +
  theme(
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.3
    ),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      color = "black",
      linewidth = 0.6
    ),
    axis.title.x = element_text(
      face = "bold",
      size = 16,
      margin = margin(t = 8)
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 16,
      margin = margin(r = 10)
    ),
    axis.text = element_text(
      size = 14,
      color = "black"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "top",
    legend.text = element_text(
      size = 14,
      color = "black"
    ),
    legend.key.width = unit(1, "cm"),
    strip.text = element_text(
      size = 14,
      face = "bold",
      color = "black"
    ),
    strip.background = element_rect(
      fill = "grey95",
      color = "white"
    )
  )

print(p)

ggsave(
  filename = "figures/Figure_2b_module_preservation.pdf",
  plot = p,
  width = 11,
  height = 5
)



#!/usr/bin/env Rscript

# =============================================================================
# Supplementary Figure S4d
# Selected GO Biological Process enrichments for changing CpG pairs (aVCPs)
#
# Input:
#   res_Heteroskid_Changing_GO.out
#
# Output:
#   Fig_S4d_changing_pair_GO_BP.pdf
#
# Five GO terms were selected manually for biological interpretation,
# consistent with the original figure.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# 1. Paths
# -----------------------------------------------------------------------------

input_file <- "/results/03_pairwise_correlations/res_Heteroskid_Changing_GO.out"
output_dir <- "/figures"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 2. Load GO enrichment results
# -----------------------------------------------------------------------------

dt <- fread(input_file)

required_cols <- c(
  "ID",
  "name",
  "Hyper_Adjp_BH",
  "Group"
)

if (!all(required_cols %in% names(dt))) {
  stop("Missing required columns in GO enrichment results.")
}

# Retain GO Biological Process results
dt <- dt[Group == "BP"]

# -----------------------------------------------------------------------------
# 3. Manually selected biological themes
# -----------------------------------------------------------------------------

selected_go <- data.table(
  ID = c(
    "GO:0060480",
    "GO:0009790",
    "GO:0060574",
    "GO:0035603",
    "GO:0060615"
  ),
  biolTheme = c(
    "Respiratory/lung development",
    "Embryonic morphogenesis",
    "Epithelial differentiation",
    "FGF signalling",
    "Mammary development"
  )
)

# Match labels by GO identifier, not row order
dt_plot <- merge(
  selected_go,
  dt,
  by = "ID",
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(dt_plot$Hyper_Adjp_BH)) {
  missing_ids <- dt_plot[is.na(Hyper_Adjp_BH), ID]

  stop(
    "Selected GO terms missing enrichment results: ",
    paste(missing_ids, collapse = ", ")
  )
}

# Check for duplicate GO identifiers
if (anyDuplicated(dt_plot$ID)) {
  stop("Duplicate selected GO terms detected in enrichment results.")
}

# -----------------------------------------------------------------------------
# 4. Transform adjusted p-values
# -----------------------------------------------------------------------------

if (any(
  !is.finite(dt_plot$Hyper_Adjp_BH) |
  dt_plot$Hyper_Adjp_BH <= 0 |
  dt_plot$Hyper_Adjp_BH > 1
)) {
  stop("Adjusted p-values must be finite and between 0 and 1.")
}

dt_plot[
  ,
  Hyper_Adjp_BH_log10 := -log10(Hyper_Adjp_BH)
]

# Order biological themes by enrichment significance
setorder(dt_plot, Hyper_Adjp_BH_log10)

dt_plot[
  ,
  biolTheme := factor(
    biolTheme,
    levels = biolTheme
  )
]

# -----------------------------------------------------------------------------
# 5. Generate Supplementary Figure S4d
# -----------------------------------------------------------------------------

p_s4d <- ggplot(
  dt_plot,
  aes(
    y = biolTheme,
    x = Hyper_Adjp_BH_log10
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = Hyper_Adjp_BH_log10,
      yend = biolTheme
    ),
    linewidth = 1.1,
    colour = "grey55"
  ) +
  geom_point(
    size = 4.2,
    shape = 21,
    fill = "black",
    colour = "black"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(
    x = expression(-log[10]~"(adj.Pvalue)"),
    y = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.y = element_text(
      size = 13,
      face = "bold",
      colour = "black"
    ),
    axis.text.x = element_text(
      size = 12,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 13,
      face = "bold"
    ),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(5.5, 40, 5.5, 5.5)
  ) +
  coord_cartesian(clip = "off")

# -----------------------------------------------------------------------------
# 6. Save figure and selected enrichment data
# -----------------------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "Fig_S4d_changing_pair_GO_BP.pdf"
  ),
  plot = p_s4d,
  width = 9,
  height = 4.5
)

fwrite(
  dt_plot[
    ,
    .(
      ID,
      name,
      biolTheme,
      Hyper_Adjp_BH,
      Hyper_Adjp_BH_log10
    )
  ],
  file.path(
    output_dir,
    "Fig_S4d_selected_GO_terms.tsv"
  ),
  sep = "\t"
)


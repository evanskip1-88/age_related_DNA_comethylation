#==============================================================================
# Visualise smoking-related variance explained
#
# Visualises the proportion of variance in smoking phenotypes explained by
# DNAm variation across age-stable and age-variable co-methylation CpG sets.
#
# R2 and its standard error are obtained from OSCA OREML analyses.
#==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

rm(list = ls())

#------------------------------------------------------------------------------
# Paths
#------------------------------------------------------------------------------

input_path <- "results/01_wgcna_combined_sex/"

#------------------------------------------------------------------------------
# Load OREML results
#------------------------------------------------------------------------------

res_osca <- fread(
  paste0(
    input_path,
    "/01a_oreml_variance_summary.tsv"
  )
)

#------------------------------------------------------------------------------
# Assign analysis labels
#------------------------------------------------------------------------------

res_osca[
  ,
  anal_labs := fcase(
    grepl("WGCNA_Changing", analysis),       "aVCMs",
    grepl("WGCNA_Stable", analysis),         "aSCMs",
    grepl("Heteroskid_Changing", analysis),  "aVCPs",
    grepl("Heteroskid_Stable", analysis),    "aSCPs",
    default = NA_character_
  )
]

#------------------------------------------------------------------------------
# Calculate 95% confidence intervals
#------------------------------------------------------------------------------

res_osca[
  ,
  `:=`(
    R2_lower = R2 - 1.96 * SE,
    R2_upper = R2 + 1.96 * SE
  )
]

res_osca[
  ,
  `:=`(
    R2 = round(R2, 4),
    R2_lower = round(R2_lower, 4),
    R2_upper = round(R2_upper, 4)
  )
]

#------------------------------------------------------------------------------
# Parse analysis names
#------------------------------------------------------------------------------

res_osca[
  ,
  c("analysis_type", "sigGroup", "trait", "Model") :=
    tstrsplit(analysis, "_")
]

res_osca[
  ,
  analysis_type := ifelse(
    analysis_type == "WGCNA",
    "WGCNA",
    "Heteroskidasticity"
  )
]

res_osca[
  ,
  analysis_type := factor(
    analysis_type,
    levels = c(
      "Heteroskidasticity",
      "WGCNA"
    )
  )
]

res_osca[
  ,
  sigGroup := factor(
    sigGroup,
    levels = c(
      "Stable",
      "Changing"
    )
  )
]

# Order CpG groups by variance explained
res_osca[
  ,
  anal_labs := factor(
    anal_labs,
    levels = unique(
      anal_labs[
        order(R2)
      ]
    )
  )
]

#------------------------------------------------------------------------------
# Plot
#------------------------------------------------------------------------------

ggplot(
  res_osca,
  aes(
    x = anal_labs,
    y = R2,
    color = sigGroup
  )
) +
  geom_point(
    size = 4
  ) +
  geom_errorbar(
    aes(
      ymin = R2_lower,
      ymax = R2_upper
    ),
    width = 0.15
  ) +
  facet_wrap(
    ~trait,
    nrow = 1
  ) +
  scale_color_manual(
    values = c(
      "Stable" = "deepskyblue",
      "Changing" = "tomato"
    )
  ) +
  labs(
    x = NULL,
    y = expression(italic(R)^2)
  ) +
  theme_bw(
    base_size = 16
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold",
      color = "black",
      size = 14
    ),
    axis.text.y = element_text(
      face = "bold",
      color = "black",
      size = 14
    ),
    strip.text = element_text(
      size = 16,
      face = "bold"
    ),
    legend.position = "none"
  )

#==============================================================================
# Visualise enrichment in PMD and solo-WCGW chromatin domains
#
# Displays enrichment odds ratios and 95% confidence intervals for
# age-stable and age-variable co-methylation CpG sets across published
# human PMD/HMD and solo-WCGW genomic regions.
#==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

rm(list = ls())

#------------------------------------------------------------------------------
# Paths
#------------------------------------------------------------------------------

input_path <- "/results/01_wgcna_combined_sex"

#------------------------------------------------------------------------------
# Load enrichment results
#------------------------------------------------------------------------------

res_enrich <- fread(
  paste0(input_path, "/01a_PMD_enrichments_human.tsv")
)

# Exclude genome-wide solo-WCGW results from this figure
res_sub <- copy(
  res_enrich[pmd_group != "soloGenomewide"]
)

#------------------------------------------------------------------------------
# Assign CpG-set labels
#------------------------------------------------------------------------------

res_sub[
  ,
  anal_group := fcase(
    grepl("WGCNA_Changing", analysis),      "aVCMs",
    grepl("WGCNA_Stable", analysis),        "aSCMs",
    grepl("Heteroskid_Changing", analysis), "aVCPs",
    grepl("Heteroskid_Stable", analysis),   "aSCPs",
    default = NA_character_
  )
]

#------------------------------------------------------------------------------
# Transform odds ratios and confidence intervals
#------------------------------------------------------------------------------

res_sub[
  ,
  `:=`(
    OR = log(OR),
    CI_lower = log(CI_lower),
    CI_upper = log(CI_upper)
  )
]

#------------------------------------------------------------------------------
# Order CpG sets and genomic domains
#------------------------------------------------------------------------------

res_sub[
  ,
  anal_group := factor(
    anal_group,
    levels = c("aSCMs", "aSCPs", "aVCMs", "aVCPs")
  )
]

res_sub[
  ,
  pmd_group := factor(
    pmd_group,
    levels = c(
      "neither",
      "commonHMD",
      "soloinCommonPMDs",
      "commonPMD"
    )
  )
]

#------------------------------------------------------------------------------
# Plot enrichment
#------------------------------------------------------------------------------

ggplot(
  res_sub,
  aes(x = OR, y = pmd_group)
) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(xmin = CI_lower, xmax = CI_upper),
    orientation = "y",
    width = 0.15,
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "firebrick",
    linewidth = 0.8
  ) +
  facet_wrap(~anal_group) +
  labs(
    x = "log(OR)",
    y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 16
    ),
    axis.text.x = element_text(
      size = 14,
      color = "black",
      face = "bold",
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(
      size = 14,
      color = "black",
      face = "bold"
    ),
    axis.title.x = element_text(
      size = 16,
      face = "bold"
    ),
    axis.title.y = element_text(
      size = 16,
      face = "bold"
    ),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(linewidth = 1)
  )

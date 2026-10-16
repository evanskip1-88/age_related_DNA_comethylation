# ============================================================
# Supplementary Figure S22a-b
#
# Sex differences in module preservation during midlife
#
# a) Male vs female absolute preservation change, Age2-Age4
# b) Preservation trajectories of sex-biased modules
#
# Reference network: Male Age1
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(readxl)
})

# ------------------------------------------------------------
# 1. Configure directories
# ------------------------------------------------------------

results_path <- "/results/02_wgcna_within_sex"
output_path <- "/results/02_wgcna_within_sex"

mod_ref <- "M"

dir.create(
  output_path,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Load module classification and midlife changes
# ------------------------------------------------------------

# Output from Part C

dt_midlife <- read_xlsx(
  file.path(
    results_path,
    "01a_sex_specific_mods_modref_M_age2v4.xlsx"
  )
)
setDT(dt_midlife)

# One observation per module
dt_midlife <- unique(
  dt_midlife[
    ,
    .(
      moduleColor,
      moduleNumID,
      change_presStats_Age2vs4_M,
      change_presStats_Age2vs4_F,
      sexChangeGroup2v4
    )
  ]
)

# ------------------------------------------------------------
# 3. Figure S22a: Midlife preservation differences
# ------------------------------------------------------------

p_a <- ggplot(
  dt_midlife,
  aes(
    x = change_presStats_Age2vs4_M,
    y = change_presStats_Age2vs4_F
  )
) +
  geom_text_repel(
    data = dt_midlife[
      sexChangeGroup2v4 == "Others"
    ],
    aes(label = moduleNumID),
    size = 4,
    fontface = "plain",
    color = "grey40",
    box.padding = 0.3,
    point.padding = 0.3,
    segment.color = NA,
    max.overlaps = Inf
  ) +
  geom_text_repel(
    data = dt_midlife[
      sexChangeGroup2v4 != "Others"
    ],
    aes(label = moduleNumID),
    size = 4,
    fontface = "bold",
    color = "grey20",
    box.padding = 0.3,
    point.padding = 0.3,
    segment.color = NA,
    max.overlaps = Inf
  ) +
  stat_cor(
    aes(label = after_stat(r.label)),
    method = "pearson",
    label.x = 0,
    label.y = 0.19,
    label.y.npc = "bottom",
    size = 6,
    color = "black"
  ) +
  labs(
    x = "\u0394 preservation (males - reference network)",
    y = "\u0394 preservation (females)"
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    plot.margin = margin(
      t = 10,
      r = 20,
      b = 40,
      l = 10
    )
  )

# ------------------------------------------------------------
# 4. Figure S22b: Preservation trajectories
# ------------------------------------------------------------

# Output from Part A
dt_pres <- fread(
  file.path(
    results_path,
    paste0(
      "age_related_module_classification_modref_",
      mod_ref, ".tsv"
    )
  )
)

# Retain preservation statistic used in original figure
dt_pres <- dt_pres[
  presStats == "corkME"
]

# Attach midlife sex-bias classification
dt_pres <- merge(
  dt_pres,
  dt_midlife[
    ,
    .(
      moduleColor,
      moduleNumID,
      sexChangeGroup2v4
    )
  ],
  by = "moduleColor"
)

# Select sex-biased modules
dt_pres <- dt_pres[
  sexChangeGroup2v4 != "Others"
]

# Format labels
dt_pres[
  , moduleLabel := paste0(
    "module_",
    moduleNumID
  )
]

dt_pres[
  , sex := fifelse(
    Sex == "Male",
    "M",
    "F"
  )
]

dt_pres[
  , testAge := factor(
    testAge,
    levels = paste0("Age", 2:8)
  )
]

# Original male/female colours
sex_colors <- c(
  "M" = "red",
  "F" = "blue"
)

p_b <- ggplot(
  dt_pres,
  aes(
    x = testAge,
    y = pres_obs_value,
    color = sex,
    group = sex
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(
    values = sex_colors
  ) +
  facet_wrap(
    ~moduleLabel,
    scales = "free_y"
  ) +
  coord_cartesian(
    ylim = c(0.50, 1.0)
  ) +
  labs(
    x = "Age group",
    y = "Network preservation",
    color = NULL
  ) +
  theme_bw() +
  theme(
    text = element_text(
      size = 18,
      family = "Arial"
    ),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1
    ),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

# ------------------------------------------------------------
# 5. Save figures
# ---------------  library(readxl)
---------------------------------------------

ggsave(
  filename = file.path(
    output_path,
    "Fig_S22a_midlife_preservation_scatter.png"
  ),
  plot = p_a,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  filename = file.path(
    output_path,
    "Fig_S22b_midlife_preservation_trajectories.png"
  ),
  plot = p_b,
  width = 12,
  height = 8,
  dpi = 300
)

cat("Supplementary Figure S22a-b completed.\n")

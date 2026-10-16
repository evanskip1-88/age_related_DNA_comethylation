# ============================================================
# Supplementary Figure S19
#
# Age-related module preservation trajectories in males
# and females, using male or female Age1 reference modules.
#
# Input: Outputs from 01_classify_age_related_modules.R
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ------------------------------------------------------------
# 1. Configure paths
# ------------------------------------------------------------

results_path <- "/results/02_wgcna_within_sex"
output_path  <- "/figures"

mod_ref <- "M"  # Choose "M" or "F"

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Load previously classified modules
# ------------------------------------------------------------

dt <- fread(
  file.path(
    results_path,
    paste0(
      "age_related_module_classification_modref_",
      mod_ref, ".tsv"
    )
  )
)

# Preserve both statistics used in the original figure
dt <- dt[
  presStats %in% c("propVarExplained", "corkME")
]

# ------------------------------------------------------------
# 3. Format classification and facet labels
# ------------------------------------------------------------

dt[, sigGroup := factor(
  as.character(sigGroup),
  levels = c("aSCMs", "aVCMs", "Others")
)]

if (mod_ref == "F") {

  dt[, Sex_lab := fifelse(
    Sex == "Male",
    "Males - female reference",
    "Females - female reference"
  )]

  sex_levels <- c(
    "Males - female reference",
    "Females - female reference"
  )

} else {

  dt[, Sex_lab := fifelse(
    Sex == "Male",
    "Males - male reference",
    "Females - male reference"
  )]

  sex_levels <- c(
    "Males - male reference",
    "Females - male reference"
  )
}

dt[, Sex_lab := factor(Sex_lab, levels = sex_levels)]

dt[, presStats := factor(
  presStats,
  levels = c("propVarExplained", "corkME")
)]

# ------------------------------------------------------------
# 4. Plot preservation trajectories
# ------------------------------------------------------------

mod_colors_grps <- c(
  "aSCMs" = "blue",
  "aVCMs" = "red",
  "Others" = "grey"
)

p <- ggplot(
  dt,
  aes(
    x = testAgeNum,
    y = pres_obs_value,
    color = sigGroup,
    group = moduleColor
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = mod_colors_grps) +
  facet_wrap(Sex_lab ~ presStats, scales = "free_y") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Test age group",
    y = "",
    color = NULL
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 18, family = "Arial"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.position = "top"
  )

# ------------------------------------------------------------
# 5. Save Figure S19
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_path,
    paste0("Fig_S19_modref_", mod_ref, ".png")
  ),
  plot = p,
  width = 12,
  height = 9,
  dpi = 300
)

cat("Supplementary Figure S19 completed:", mod_ref, "\n")

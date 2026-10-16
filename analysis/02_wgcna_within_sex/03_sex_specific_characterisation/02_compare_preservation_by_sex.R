# ============================================================
# Compare module preservation between males and females
#
# Input: Module classifications from Part A
# Reference modules: Male or female Age1 network
# Test age groups: Age2–Age8
#
# Outputs:
#   1. Paired male/female preservation statistics
#   2. Age-specific preservation correlations
#   3. Male vs female preservation scatterplots
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggpubr)
})

# ------------------------------------------------------------
# 1. Configure directories
# ------------------------------------------------------------

data_path <- "results/02_wgcna_within_sex/"
output_path <- "/figures"

mod_ref <- "M"  # Reference network: "M" or "F"

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Load classified preservation statistics from Part A
# ------------------------------------------------------------

input_file <- file.path(
  data_path,
  paste0(
    "age_related_module_trajectories_modref_",
    mod_ref, ".tsv"
  )
)

dt <- fread(input_file)

# Select the corkME statistic used in the original
# male-female correlation analysis
dt <- dt[presStats == "corkME"]

# ------------------------------------------------------------
# 3. Match male and female observations
# ------------------------------------------------------------

join_cols <- c(
  "moduleColor",
  "testAge",
  "presStats",
  "module_ref"
)

# Retain relevant module-level information
keep_cols <- c(
  join_cols,
  "testAgeNum",
  "moduleSize",
  "pres_obs_value",
  "change_presStats",
  "sigGroup"
)

dt_male <- dt[
  Sex == "Male",
  ..keep_cols
]

dt_female <- dt[
  Sex == "Female",
  ..keep_cols
]

# Rename sex-specific columns
setnames(
  dt_male,
  setdiff(names(dt_male), join_cols),
  paste0(setdiff(names(dt_male), join_cols), "_M")
)

setnames(
  dt_female,
  setdiff(names(dt_female), join_cols),
  paste0(setdiff(names(dt_female), join_cols), "_F")
)

# Match observations by module and age group
dt_paired <- merge(
  dt_male,
  dt_female,
  by = join_cols,
  all = FALSE
)

setorder(dt_paired, testAgeNum_M, moduleColor)

# ------------------------------------------------------------
# 4. Calculate age-specific male-female correlations
# ------------------------------------------------------------

dt_correlations <- dt_paired[
  ,
  .(
    n_modules = sum(
      complete.cases(pres_obs_value_M, pres_obs_value_F)
    ),
    correlation = cor(
      pres_obs_value_M,
      pres_obs_value_F,
      use = "complete.obs",
      method = "pearson"
    )
  ),
  by = .(testAge, testAgeNum_M, module_ref)
]

setorder(dt_correlations, testAgeNum_M)

print(dt_correlations)

# ------------------------------------------------------------
# 5. Plot male versus female preservation
# ------------------------------------------------------------

if (mod_ref == "M") {

  x_labs <- "cor.kME (males - reference network)"
  y_labs <- "cor.kME (females)"

} else {

  x_labs <- "cor.kME (males)"
  y_labs <- "cor.kME (females - reference network)"
}

corr_plot <- ggplot(
  dt_paired,
  aes(
    x = pres_obs_value_M,
    y = pres_obs_value_F
  )
) +
  geom_point(
    aes(color = testAge),
    size = 2
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "black",
    linetype = "dashed"
  ) +
  stat_cor(
    method = "pearson",
    aes(label = after_stat(r.label)),
    label.x.npc = "left",
    label.y.npc = "top",
    size = 4
  ) +
  facet_wrap(
    ~testAge,
    scales = "free",
    ncol = 7
  ) +
  coord_cartesian(
    xlim = c(0.4, 1.0),
    ylim = c(0.4, 1.0)
  ) +
  labs(
    x = x_labs,
    y = y_labs
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 16),
    legend.position = "none",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )

print(corr_plot)

# ------------------------------------------------------------
# 6. Save outputs
# ------------------------------------------------------------

fwrite(
  dt_paired,
  file.path(
    output_path,
    paste0(
      "paired_preservation_by_sex_modref_",
      mod_ref, ".tsv"
    )
  ),
  sep = "\t"
)

fwrite(
  dt_correlations,
  file.path(
    output_path,
    paste0(
      "preservation_correlations_by_age_modref_",
      mod_ref, ".tsv"
    )
  ),
  sep = "\t"
)

ggsave(
  filename = file.path(
    output_path,
    paste0(
      "preservation_correlation_plot_modref_",
      mod_ref, ".png"
    )
  ),
  plot = corr_plot,
  width = 18,
  height = 5,
  dpi = 300
)

cat("Male-female preservation comparison completed.\n")

# ============================================================
# Classify age-related module preservation trajectories by sex
#
# Reference modules: male Age1 network
# Test datasets: male and female Age1–Age8
#
# Classification:
#   Stable   -> aSCMs
#   Changing -> aVCMs
#   Others   -> modules not meeting selection criteria
#
# Preservation statistics:
#   propVarExplained and corkME
# ============================================================

suppressPackageStartupMessages({
  library(mgcv)
  library(data.table)
  library(ggplot2)
})

# ------------------------------------------------------------
# 1. Configure directories and analysis parameters
# ------------------------------------------------------------

data_path <- "/results/02_wgcna_within_sex"
output_path <- "/results/02_wgcna_within_sex"

mod_ref <- "M"  # Reference modules: "M" or "F"

connectionCutoff <- 0.40
gamPvalCutoff <- 0.05
stableSDcutoff <- 0.20
changeQuantCutoff <- 0.50
changePercentCutoff <- 0.10

preservation_stats <- c("propVarExplained", "corkME")

# ------------------------------------------------------------
# 2. Classify modules using preservation trajectories
# ------------------------------------------------------------

select_top_modules <- function(
    sex_id,
    mod_ref,
    pres_stats,
    connectionCutoff,
    gamPvalCutoff,
    stableSDcutoff,
    changeQuantCutoff,
    changePercentCutoff
) {

  # Load observed preservation statistics
  obs_file <- file.path(
    data_path,
    paste0(
      "precorr_gs_qc_chr_all_refage1_MinMod20_",
      sex_id, "_modLabref_", mod_ref, ".obsStat"
    )
  )

  dt_obs_stats <- fread(obs_file)

  # Exclude randomly generated gold module
  dt_obs_stats <- dt_obs_stats[moduleColor != "gold"]

  # Convert test age group to numeric
  dt_obs_stats[
    , testAgeNum := as.numeric(sub("Age", "", testAge))
  ]

  # Original preservation statistic names
  setnames(
    dt_obs_stats,
    old = names(dt_obs_stats)[6:8],
    new = c("corkIM", "corkME", "corCor")
  )

  # Fit GAM to each module's age-related trajectory
  gam_results <- lapply(
    unique(dt_obs_stats$moduleColor),
    function(mod) {

      module_data <- dt_obs_stats[moduleColor == mod]

      gam_model <- gam(
        as.formula(
          paste0(pres_stats, " ~ s(testAgeNum, k = 4)")
        ),
        data = module_data
      )

      gam_summary <- summary(gam_model)

      data.table(
        moduleColor = mod,
        edf = gam_summary$s.table[1, "edf"],
        F_value = gam_summary$s.table[1, "F"],
        gam_p_value = gam_summary$s.table[1, "p-value"],
        R2_adj = gam_summary$r.sq
      )
    }
  )

  dt_gam_stats <- rbindlist(gam_results)

  # Initial classification based on GAM p-value
  dt_gam_stats[
    , sigGroup := fifelse(
      gam_p_value < gamPvalCutoff,
      "Changing",
      "Stable"
    )
  ]

  dt_obs_stats <- merge(
    dt_obs_stats,
    dt_gam_stats[
      , .(moduleColor, edf, F_value, gam_p_value, R2_adj, sigGroup)
    ],
    by = "moduleColor",
    all.x = TRUE
  )

  # ----------------------------------------------------------
  # 3. Apply secondary selection criteria
  # ----------------------------------------------------------

  # Standard deviation across Age2–Age8
  sd_values <- dt_obs_stats[
    , .(
      sd_presStats = sd(get(pres_stats), na.rm = TRUE)
    ),
    by = moduleColor
  ]

  # Absolute preservation change between Age2 and Age8
  extreme_changes <- dt_obs_stats[
    testAge %in% c("Age2", "Age8")
  ][
    order(moduleColor, testAgeNum)
  ][
    , .(
      change_presStats = abs(diff(get(pres_stats)))
    ),
    by = moduleColor
  ]

  dt_obs_stats <- merge(
    dt_obs_stats,
    sd_values,
    by = "moduleColor",
    all.x = TRUE
  )

  dt_obs_stats <- merge(
    dt_obs_stats,
    extreme_changes,
    by = "moduleColor",
    all.x = TRUE
  )

  # Original quantile-based thresholds
  sd_cutoff <- quantile(
    dt_obs_stats$sd_presStats,
    stableSDcutoff,
    na.rm = TRUE
  )

  change_cutoff <- max(
    quantile(
      dt_obs_stats$change_presStats,
      changeQuantCutoff,
      na.rm = TRUE
    ),
    changePercentCutoff
  )

  # Modules with sufficiently strong preservation at Age2
  connected_mods <- dt_obs_stats[
    testAge == "Age2" & get(pres_stats) > connectionCutoff,
    unique(moduleColor)
  ]

  # Secondary classification filters
  dt_obs_stats[
    sigGroup == "Stable" &
      sd_presStats >= sd_cutoff,
    sigGroup := "Others"
  ]

  dt_obs_stats[
    sigGroup == "Stable" &
      !(moduleColor %in% connected_mods),
    sigGroup := "Others"
  ]

  dt_obs_stats[
    sigGroup == "Changing" &
      change_presStats < change_cutoff,
    sigGroup := "Others"
  ]

  dt_obs_stats[
    sigGroup == "Changing" &
      !(moduleColor %in% connected_mods),
    sigGroup := "Others"
  ]

  # Rename selected groups
  dt_obs_stats[
    sigGroup == "Stable",
    sigGroup := "aSCMs"
  ]

  dt_obs_stats[
    sigGroup == "Changing",
    sigGroup := "aVCMs"
  ]

  # Identify preservation statistic
  dt_obs_stats[, presStats := pres_stats]

  # Retain original long-format structure for downstream work
  dt_obs_stats[
    , pres_obs_value := get(pres_stats)
  ]

  return(dt_obs_stats)
}

# ------------------------------------------------------------
# 4. Run classification for males and females
# ------------------------------------------------------------

results <- list()

for (sex_id in c("M", "F")) {

  for (pres_stats in preservation_stats) {

    cat(
      "Processing sex:", sex_id,
      "| Statistic:", pres_stats, "\n"
    )

    res <- select_top_modules(
      sex_id = sex_id,
      mod_ref = mod_ref,
      pres_stats = pres_stats,
      connectionCutoff = connectionCutoff,
      gamPvalCutoff = gamPvalCutoff,
      stableSDcutoff = stableSDcutoff,
      changeQuantCutoff = changeQuantCutoff,
      changePercentCutoff = changePercentCutoff
    )

    res[, Sex := ifelse(sex_id == "M", "Male", "Female")]

    results[[paste(sex_id, pres_stats, sep = "_")]] <- res
  }
}

# Combine results
sel_mods_c_all_sex <- rbindlist(
  results,
  use.names = TRUE,
  fill = TRUE
)

sel_mods_c_all_sex[, module_ref := mod_ref]

# Reference network labels
sel_mods_c_all_sex[
  , Sex_lab := paste0(
    ifelse(Sex == "Male", "Males", "Females"),
    " - ",
    ifelse(mod_ref == "M", "male", "female"),
    " reference"
  )
]

# ------------------------------------------------------------
# 5. Plot preservation trajectories
# ------------------------------------------------------------

mod_colors_grps <- c(
  "aSCMs" = "blue",
  "aVCMs" = "red",
  "Others" = "grey"
)

sel_mods_c_all_sex[
  , sigGroup := factor(
    sigGroup,
    levels = c("aSCMs", "aVCMs", "Others")
  )
]

p <- ggplot(
  sel_mods_c_all_sex,
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
  facet_grid(Sex_lab ~ presStats, scales = "free_y") +
  labs(
    x = "Test age group",
    y = "Preservation value",
    color = NULL
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

print(p)

# ------------------------------------------------------------
# 6. Save results
# ------------------------------------------------------------

fwrite(
  sel_mods_c_all_sex,
  file.path(
    output_path,
    paste0(
      "age_related_module_trajectories_modref_",
      mod_ref, ".tsv"
    )
  ),
  sep = "\t"
)

ggsave(
  filename = file.path(
    output_path,
    paste0(
      "age_related_module_trajectories_modref_",
      mod_ref, ".png"
    )
  ),
  plot = p,
  width = 12,
  height = 7,
  dpi = 300
)

cat("Module classification completed.\n")

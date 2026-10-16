#!/usr/bin/env Rscript

# =============================================================================
# Supplementary Figure S1
# Sex composition and age distribution across eight age groups
#
# Inputs:
#   N_samples_male_female_witinSex.txt
#   male_female_N7532.txt
#
# Output:
#   Fig_S1_age_sex_distribution.pdf
#
# Additional output:
#   KS_test_age_distributions.txt
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggpubr)
})

# -----------------------------------------------------------------------------
# 1. Configure paths
# -----------------------------------------------------------------------------

input_dir <- "..."
output_dir <- "..."

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

count_file <- file.path(
  input_dir,
  "N_samples_male_female_witinSex.txt"
)

individual_file <- file.path(
  input_dir,
  "male_female_N7532.txt"
)

age_groups <- paste0("Age", 1:8)

sex_colours <- c(
  Male = "#1f78b4",
  Female = "#ff7f00"
)

# -----------------------------------------------------------------------------
# 2. Load sample counts by sex and age group
# -----------------------------------------------------------------------------

dt_counts <- fread(count_file)

stopifnot(
  all(c("Var1", "sex", "Freq") %in% names(dt_counts))
)

dt_counts[
  ,
  sex := fifelse(
    tolower(sex) == "male",
    "Male",
    "Female"
  )
]

dt_counts[
  ,
  Var1 := factor(Var1, levels = age_groups)
]

dt_counts[
  ,
  sex := factor(sex, levels = c("Male", "Female"))
]

# -----------------------------------------------------------------------------
# 3. Panel A: Number of participants by sex
# -----------------------------------------------------------------------------

plt_bar <- ggplot(
  dt_counts,
  aes(
    x = Var1,
    y = Freq,
    fill = sex
  )
) +
  geom_col(
    position = position_dodge(width = 0.9)
  ) +
  scale_fill_manual(
    values = sex_colours
  ) +
  labs(
    x = "",
    y = "Number of samples",
    fill = "Sex"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    legend.title = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

# -----------------------------------------------------------------------------
# 4. Load individual-level age and sex data
# -----------------------------------------------------------------------------

dt_individuals <- fread(individual_file)

stopifnot(
  all(c("age", "sex") %in% names(dt_individuals))
)

# Reproduce original age-group boundaries

dt_individuals[
  ,
  age_grp := fcase(
    age >= 17.9 & age < 36.0, "Age1",
    age >= 36.0 & age < 44.0, "Age2",
    age >= 44.0 & age < 49.5, "Age3",
    age >= 49.5 & age < 54.2, "Age4",
    age >= 54.2 & age < 57.8, "Age5",
    age >= 57.8 & age < 61.5, "Age6",
    age >= 61.5 & age < 66.0, "Age7",
    age >= 66.0, "Age8",
    default = NA_character_
  )
]

dt_individuals[
  ,
  sex := fifelse(
    sex == "M",
    "Male",
    "Female"
  )
]

dt_individuals[
  ,
  age_grp := factor(age_grp, levels = age_groups)
]

dt_individuals[
  ,
  sex := factor(sex, levels = c("Male", "Female"))
]

# -----------------------------------------------------------------------------
# 5. Panel B: Age distributions by sex
# -----------------------------------------------------------------------------

plt_violin <- ggplot(
  dt_individuals,
  aes(
    x = age_grp,
    y = age,
    fill = sex
  )
) +
  geom_violin(
    trim = FALSE,
    position = position_dodge(width = 1)
  ) +
  geom_boxplot(
    width = 0.1,
    position = position_dodge(width = 1)
  ) +
  scale_fill_manual(
    values = sex_colours
  ) +
  labs(
    x = "",
    y = "Age (years)",
    fill = "Sex"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    legend.title = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

# -----------------------------------------------------------------------------
# 6. Combine panels and export Figure S1
# -----------------------------------------------------------------------------

fig_s1 <- ggarrange(
  plt_bar,
  plt_violin,
  ncol = 2,
  common.legend = TRUE,
  legend = "right"
)

ggsave(
  filename = file.path(
    output_dir,
    "Fig_S1_age_sex_distribution.pdf"
  ),
  plot = fig_s1,
  width = 13,
  height = 5.5
)

# -----------------------------------------------------------------------------
# 7. Kolmogorov-Smirnov tests
# -----------------------------------------------------------------------------

ks_results <- capture.output({

  cat("Kolmogorov-Smirnov tests: male versus female age distributions\n")
  cat("=============================================================\n\n")

  for (i in seq_along(age_groups)) {

    group <- age_groups[i]

    female_ages <- dt_individuals[
      sex == "Female" & age_grp == group,
      age
    ]

    male_ages <- dt_individuals[
      sex == "Male" & age_grp == group,
      age
    ]

    cat("\n", group, "\n", sep = "")

    if (length(male_ages) > 0 && length(female_ages) > 0) {
      print(ks.test(male_ages, female_ages))
    } else {
      cat("Insufficient observations.\n")
    }
  }

  cat("\nOverall male versus female comparison\n")

  print(
    ks.test(
      dt_individuals[sex == "Female", age],
      dt_individuals[sex == "Male", age]
    )
  )
})

writeLines(
  ks_results,
  file.path(
    output_dir,
    "KS_test_age_distributions.txt"
  )
)


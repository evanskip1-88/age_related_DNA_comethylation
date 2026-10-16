# ============================================================
# Supplementary Figure S20
#
# Distribution of absolute sex differences in the magnitude
# of module preservation changes:
#
# a) Age2 vs Age8, female reference modules
# b) Age2 vs Age8, male reference modules
# c) Age2 vs Age4, male reference modules
#
# Dashed red line: sex-difference threshold = 0.05
# ============================================================

library(data.table)

# Directories
data_path <- "/results/02_wgcna_within_sex"
output_path <- "/figures"

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# Parameters
sexDifferenceCutoff <- 0.05

# ------------------------------------------------------------
# Load module summaries from Part C
# ------------------------------------------------------------

dt_F <- fread(file.path(
  data_path,
  "01a_sex_specific_mods_modref_F.tsv"
))

dt_M <- fread(file.path(
  data_path,
  "01a_sex_specific_mods_modref_M.tsv"
))

# One observation per module
dt_F <- unique(dt_F[
  , .(moduleColor, moduleNumID, abs_diff_change)
])

dt_M <- unique(dt_M[
  , .(
    moduleColor,
    moduleNumID,
    abs_diff_change,
    abs_diff_change2v4
  )
])

# ------------------------------------------------------------
# Generate Supplementary Figure S20
# ------------------------------------------------------------

png(
  filename = file.path(output_path, "Fig_S20.png"),
  width = 2400,
  height = 1800,
  res = 300
)

par(mfrow = c(2, 2))

plot_histogram <- function(values, panel_label) {

  values <- values[is.finite(values)]

  hist(
    values,
    main = "",
    xlab = expression("|" * Delta * "Pres. Stats" * "|"),
    ylab = "Num. of modules",
    cex.lab = 1.5,
    cex.axis = 1.5
  )

01a_sex_specific_mods_modref_F  abline(
    v = sexDifferenceCutoff,
    col = "red",
    lwd = 2,
    lty = 2
  )

  mtext(
    panel_label,
    side = 3,
    line = 1,
    adj = 0,
    cex = 1.5,
    font = 2
  )

  # Print distribution summaries
  cat("\nPanel", panel_label, "\n")
  print(quantile(values))

  cat(
    "Proportion below 0.05:",
    ecdf(values)(sexDifferenceCutoff),
    "\n"
  )
}

# a) Female reference, Age2 vs Age8
plot_histogram(dt_F$abs_diff_change, "a")

# b) Male reference, Age2 vs Age8
plot_histogram(dt_M$abs_diff_change, "b")

# c) Male reference, Age2 vs Age4
plot_histogram(dt_M$abs_diff_change2v4, "c")

dev.off()

cat("\nSupplementary Figure S20 generated.\n")

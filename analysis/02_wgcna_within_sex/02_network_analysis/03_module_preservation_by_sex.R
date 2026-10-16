# ============================================================
# Assess WGCNA module preservation across age groups
# within and between sexes
#
# Reference modules: Age1, selected by mod_ref
# Test DNAm data: selected by sex_id
# ============================================================

suppressPackageStartupMessages({
  library(WGCNA)
  library(data.table)
})

# Configure directories
dnam_inpath <- "..."
output_path <- "/results/02_wgcna_within_sex"

# Analysis settings
sex_id <- "F"  # Sex of DNAm data: "M" or "F"
mod_ref <- "M" # Sex defining reference modules: "M" or "F"
modSize <- 20

# Reference DNAm data (Age1)
ref_fname <- file.path(
  dnam_inpath,
  paste0("precorr_gs_qc_chr_all_age_grp1_", sex_id, ".txt")
)

# Reference module assignments
mod_file <- file.path(
  output_path,
  paste0(
    "precorr_gs_qc_chr_all_refage1_MinMod",
    modSize, "_", mod_ref, ".modLabs"
  )
)

dt_modules <- fread(mod_file)

# Select modules containing 20–4,999 CpGs
dt_module_counts <- dt_modules[
  , .(probeCounts = .N),
  by = moduleNumID
]

dt_select_mods <- dt_modules[
  moduleNumID %in% dt_module_counts[
    probeCounts >= 20 & probeCounts < 5000,
    moduleNumID
  ]
]

# Load reference DNAm once
ref_dnam <- as.matrix(
  fread(ref_fname, select = dt_select_mods$probe_id)
)

# Assess module preservation across Age2–Age8
res_pres_lst <- vector("list", 7)

for (ageGrp in 2:8) {

  cat("Processing test age group:", ageGrp, "\n")

  test_fname <- file.path(
    dnam_inpath,
    paste0(
      "precorr_gs_qc_chr_all_age_grp",
      ageGrp, "_", sex_id, ".txt"
    )
  )

  # Load test DNAm for selected CpGs
  test_dnam <- as.matrix(
    fread(test_fname, select = dt_select_mods$probe_id)
  )

  # Calculate observed module preservation statistics
  preserv_res <- modulePreservation(
    multiData = list(
      ref = list(data = ref_dnam),
      test = list(data = test_dnam)
    ),
    multiColor = list(
      ref = dt_select_mods$moduleColors
    ),
    dataIsExpr = TRUE,
    nPermutations = 0,
    verbose = 5
  )

  # Extract observed preservation statistics
  res_pres_obs <- data.table(
    moduleColor = as.factor(
      rownames(
        preserv_res$observed$ref$inter$ref_vs_test
      )
    ),
    refAge = "Age1",
    testAge = paste0("Age", ageGrp),

    cbind(
      preserv_res$observed$ref$quality[[2]][, 1:2],
      preserv_res$observed$ref$intra$ref_vs_test[, c(2, 1)],
      preserv_res$observed$ref$inter$ref_vs_test[, c(1, 2, 4)]
    )
  )

  # Retain the original selected preservation statistics
  res_pres_obs[
    , c("propVarExplained.qual", "meanSignAwareKME.pres") := NULL
  ]

  setnames(
    res_pres_obs,
    old = names(res_pres_obs)[5],
    new = "propVarExplained"
  )

  res_pres_lst[[ageGrp - 1]] <- res_pres_obs

  cat("Completed age group:", ageGrp, "\n")
}

# Combine observed preservation statistics
dt_res_pres_obs <- rbindlist(res_pres_lst)

# Round preservation statistics to four decimal places
numeric_cols <- names(dt_res_pres_obs)[
  vapply(dt_res_pres_obs, is.numeric, logical(1))
]

dt_res_pres_obs[
  , (numeric_cols) := lapply(.SD, round, digits = 4),
  .SDcols = numeric_cols
]

# Save results
out_file <- file.path(
  output_path,
  paste0(
    "precorr_gs_qc_chr_all_refage1_MinMod",
    modSize, "_", sex_id,
    "_modLabref_", mod_ref, ".obsStat"
  )
)

fwrite(
  dt_res_pres_obs,
  out_file,
  sep = "\t",
  quote = FALSE
)

cat("Module preservation analysis completed.\n")

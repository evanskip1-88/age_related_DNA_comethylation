
#!/usr/bin/env Rscript

# =============================================================================
# Consolidate rGREAT functional enrichment results
#
# Combines GO, phenotype and enriched-gene results across Stable
# and Changing co-methylation modules.
#
# Output:
#   all_res_rGREAT_out.tsv
#
# Subsequent biological theme curation was performed manually in Excel.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# 1. Paths
# -----------------------------------------------------------------------------

input_dir <- "/results/01_wgcna_combined_sex/functional_enrichment"
output_dir <- "/results/01_wgcna_combined_sex"

output_file <- file.path(output_dir, "all_res_rGREAT_out.tsv")

# -----------------------------------------------------------------------------
# 2. Load rGREAT output files
# -----------------------------------------------------------------------------

files <- list.files(
  path = input_dir,
  pattern = "\\.out$",
  full.names = TRUE
)

if (length(files) == 0L) {
  stop("No rGREAT .out files found in: ", input_dir)
}

load_rgreat_file <- function(file) {

  dt <- fread(file)

  if (nrow(dt) == 0L) {
    return(NULL)
  }

  dt[, ID := as.character(ID)]

  # Identify Stable or Changing analysis from filename
  sig_group <- if (grepl("Stable", basename(file))) {
    "Stable"
  } else if (grepl("Changing", basename(file))) {
    "Changing"
  } else {
    NA_character_
  }

  dt[, sigGroup := sig_group]

  dt
}

res_out <- rbindlist(
  lapply(files, load_rgreat_file),
  use.names = TRUE,
  fill = TRUE
)

# -----------------------------------------------------------------------------
# 3. Prepare enrichment statistics
# -----------------------------------------------------------------------------

res_out[, Hyper_Adjp_BH := as.numeric(Hyper_Adjp_BH)]

# Rounded adjusted P-value, as in the original script
res_out[, Hyper_Adjp_BH_2 := as.numeric(
  sprintf("%.1e", Hyper_Adjp_BH)
)]

# Negative log10 adjusted P-value
res_out[, Hyper_Adjp_BH_log10 := -log10(Hyper_Adjp_BH)]

# Analysis label used by downstream scripts
res_out[, analysis := paste0("WGCNA_", sigGroup)]

# -----------------------------------------------------------------------------
# 4. Save consolidated results
# -----------------------------------------------------------------------------

fwrite(
  res_out,
  output_file,
  sep = "\t",
  quote = FALSE
)

message("Consolidated rGREAT results saved: ", output_file)


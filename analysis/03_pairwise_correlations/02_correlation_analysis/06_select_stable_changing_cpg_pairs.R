#!/usr/bin/env Rscript

# =============================================================================
# Select changing and stable CpG pairs for genomic annotation
#
# Inputs:
#   1. Combined pairwise correlation table
#   2. Combined heteroskedasticity results
#   3. rGREAT CpG-to-gene annotation
#
# Outputs:
#   Selected CpG pairs
#   Unique CpG lists for stable and changing groups
#
# IMPORTANT:
#   This version defines one result per unordered CpG pair using
#   the minimum p-value across its heteroskedasticity records.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# 1. Paths and parameters
# -----------------------------------------------------------------------------

corr_file <- "/results/03_pairwise_correlations/precorr_gs_qc_chr_all_Ages_vsOthers_long.cor"
heterosk_file <- "/results/03_pairwise_correlations/chr_all_AgeBaseAll_age_beta_est_resiq_sq.res"
rgreat_file <- "/data/gsDNAm/GS_sub_qc1_gene_cpgs_assoc_rGREAT_h19_all_CpGs.opi"
output_dir <- "/results/03_pairwise_correlations"

top_n <- 1000L
alpha <- 0.05

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 2. Load input tables
# -----------------------------------------------------------------------------

dt_corr <- fread(corr_file)
dt_het <- fread(heterosk_file)
opi_rgreat <- fread(rgreat_file)

stopifnot(
  all(c("probe1", "probe2", "AgeBase") %in% names(dt_corr)),
  all(c("probe1", "probe2", "p_val") %in% names(dt_het)),
  all(c("probe_id", "rgreat_genes", "position", "chr") %in%
        names(opi_rgreat))
)

# -----------------------------------------------------------------------------
# 3. Create direction-independent pair identifiers
# -----------------------------------------------------------------------------

make_pair_id <- function(x, y) {
  paste(pmin(x, y), pmax(x, y), sep = "_")
}

dt_corr[, probe_pair := make_pair_id(probe1, probe2)]
dt_het[, probe_pair := make_pair_id(probe1, probe2)]

# Retain only pairs represented in the correlation results
eligible_pairs <- unique(dt_corr$probe_pair)

dt_het <- dt_het[probe_pair %in% eligible_pairs]

# -----------------------------------------------------------------------------
# 4. Aggregate heteroskedasticity results per CpG pair
# -----------------------------------------------------------------------------

dt_het <- dt_het[is.finite(p_val) & p_val >= 0 & p_val <= 1]

# Keep the row with the smallest p-value for each pair.
# Ties are resolved deterministically using AgeBase.

if ("AgeBase" %in% names(dt_het)) {
  setorder(dt_het, probe_pair, p_val, AgeBase)
} else {
  setorder(dt_het, probe_pair, p_val)
}

dt_pairs <- dt_het[, .SD[1L], by = probe_pair]

n_tests <- uniqueN(dt_het$probe_pair)

if (n_tests == 0L) {
  stop("No eligible heteroskedasticity results found.")
}

bonferroni_threshold <- alpha / n_tests

message("Unique eligible CpG pairs: ", n_tests)
message("Bonferroni threshold: ", bonferroni_threshold)

# -----------------------------------------------------------------------------
# 5. Select changing pairs
# -----------------------------------------------------------------------------

dt_changing <- copy(
  dt_pairs[p_val < bonferroni_threshold]
)

setorder(dt_changing, p_val, probe_pair)

dt_changing <- head(dt_changing, top_n)

dt_changing[, sigGroup := "Changing"]

# -----------------------------------------------------------------------------
# 6. Select stable pairs
# -----------------------------------------------------------------------------

dt_stable <- copy(
  dt_pairs[p_val > bonferroni_threshold]
)

# Highest p-values first
setorder(dt_stable, -p_val, probe_pair)

dt_stable <- head(dt_stable, top_n)

dt_stable[, sigGroup := "Stable"]

# -----------------------------------------------------------------------------
# 7. Combine selected pairs
# -----------------------------------------------------------------------------

dt_selected <- rbindlist(
  list(dt_stable, dt_changing),
  fill = TRUE
)

dt_selected[, Analysis := "Heteroskedasticity"]

message("Selected stable pairs: ", nrow(dt_stable))
message("Selected changing pairs: ", nrow(dt_changing))

# -----------------------------------------------------------------------------
# 8. Annotate probe1 and probe2
# -----------------------------------------------------------------------------

probe_map <- unique(
  opi_rgreat[
    ,
    .(
      probe_id,
      gene = rgreat_genes,
      position,
      chr
    )
  ],
  by = "probe_id"
)

map1 <- copy(probe_map)

setnames(
  map1,
  c("probe_id", "gene", "position", "chr"),
  c("probe1", "probe1_gene", "probe1_pos", "probe1_chr")
)

map2 <- copy(probe_map)

setnames(
  map2,
  c("probe_id", "gene", "position", "chr"),
  c("probe2", "probe2_gene", "probe2_pos", "probe2_chr")
)

dt_selected <- merge(
  dt_selected,
  map1,
  by = "probe1",
  all.x = TRUE,
  sort = FALSE
)

dt_selected <- merge(
  dt_selected,
  map2,
  by = "probe2",
  all.x = TRUE,
  sort = FALSE
)

# -----------------------------------------------------------------------------
# 9. Export selected pair table
# -----------------------------------------------------------------------------

pair_columns <- c(
  "probe_pair",
  "probe1",
  "probe2",
  "probe1_chr",
  "probe2_chr",
  "probe1_pos",
  "probe2_pos",
  "probe1_gene",
  "probe2_gene",
  "beta",
  "beta.se",
  "t_val",
  "p_val",
  "sigGroup",
  "Analysis"
)

pair_columns <- intersect(pair_columns, names(dt_selected))

fwrite(
  dt_selected[, ..pair_columns],
  file.path(output_dir, "10a_top_1000_stable_changing.tsv"),
  sep = "\t"
)

# -----------------------------------------------------------------------------
# 10. Export unique CpGs for downstream annotation
# -----------------------------------------------------------------------------

for (group in c("Stable", "Changing")) {

  dt_group <- dt_selected[sigGroup == group]

  probe_ids <- sort(unique(c(
    dt_group$probe1,
    dt_group$probe2
  )))

  probe_ids <- probe_ids[!is.na(probe_ids)]

  dt_probes <- data.table(
    probe_id = probe_ids,
    sigGroup = group
  )

  output_file <- file.path(
    output_dir,
    paste0(
      "chr_all_AgeBaseAll_",
      tolower(group),
      "PosNeg_heteroskid_top1000_cpgs_list.tsv"
    )
  )

  fwrite(
    dt_probes,
    output_file,
    sep = "\t"
  )

  message(group, ": ", length(probe_ids), " unique CpGs")
}

message("Stable/changing CpG pair selection completed.")

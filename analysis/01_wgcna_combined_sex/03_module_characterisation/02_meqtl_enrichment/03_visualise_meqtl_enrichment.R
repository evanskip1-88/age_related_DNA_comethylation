#==============================================================================
# Visualise GoDMC meQTL enrichment
#
# Visualises enrichment of age-stable and age-variable co-methylation
# modules and probe pairs among cis- and trans-meQTL-associated CpGs.
#==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

rm(list = ls())

#------------------------------------------------------------------------------
# Path
#------------------------------------------------------------------------------

data_path <- "..."

#==============================================================================
# Load enrichment results
#==============================================================================

meqtls <- fread(
  paste0(
    data_path,
    "/01a_res_by_study_meqtl_enrichments.tsv"
  )
)

# Exclude unclassified WGCNA modules
meqtls <- meqtls[
  !grepl("Others", analysis)
]

# Rename QTL type column
setnames(
  meqtls,
  "eqtl_type",
  "qtl_type"
)

#==============================================================================
# Prepare enrichment statistics
#==============================================================================

# Plot odds ratios and confidence intervals on the natural-log scale
meqtls[
  ,
  `:=`(
    OR = log(OR),
    CI_lower = log(CI_lower),
    CI_upper = log(CI_upper)
  )
]

# OR and 95% CI labels
meqtls[
  ,
  OR_CI_label := ifelse(
    is.na(OR) |
      is.na(CI_lower) |
      is.na(CI_upper),
    "NA",
    paste0(
      round(OR, 2),
      "\n(",
      round(CI_lower, 2),
      "-",
      round(CI_upper, 2),
      ")"
    )
  )
]

# Format P-values for plotting
meqtls[
  ,
  pval_fmt := ifelse(
    is.na(pvalue),
    NA_character_,
    gsub(
      "e([+-]?[0-9]+)",
      " %*% 10^{\\1}",
      formatC(
        pvalue,
        format = "e",
        digits = 1
      )
    )
  )
]

#==============================================================================
# Construct heatmap matrices
#==============================================================================

row_names_order <- c(
  "WGCNA_Stable",
  "WGCNA_Changing",
  "Heteroskid_Stable",
  "Heteroskid_Changing"
)

col_names_order <- c(
  "Cis",
  "Trans"
)

row_names_labs <- c(
  "aSCMs",
  "aVCMs",
  "aSCPs",
  "aVCPs"
)

col_names_labs <- c(
  "cis-meQTLs",
  "trans-meQTLs"
)

row_split_labs <- factor(
  c(
    "modules",
    "modules",
    "pairs",
    "pairs"
  ),
  levels = c(
    "modules",
    "pairs"
  )
)

# Log odds ratios
or_matrix <- reshape2::acast(
  meqtls,
  analysis ~ qtl_type,
  value.var = "OR"
)

or_matrix <- or_matrix[
  row_names_order,
  col_names_order
]

# Log odds ratios and confidence intervals
ci_matrix <- reshape2::acast(
  meqtls,
  analysis ~ qtl_type,
  value.var = "OR_CI_label",
  fill = ""
)

ci_matrix <- ci_matrix[
  row_names_order,
  col_names_order
]

ci_matrix[
  ci_matrix == ""
] <- NA

# P-values
pval_matrix <- reshape2::acast(
  meqtls,
  analysis ~ qtl_type,
  value.var = "pval_fmt"
)

pval_matrix <- pval_matrix[
  row_names_order,
  col_names_order
]

#==============================================================================
# Colour scale
#==============================================================================

log_min <- floor(
  min(
    or_matrix,
    na.rm = TRUE
  )
)

log_max <- ceiling(
  max(
    or_matrix,
    na.rm = TRUE
  )
)

mat_col <- colorRamp2(
  c(
    log_min,
    0,
    log_max
  ),
  c(
    "skyblue",
    "white",
    "red"
  )
)

#==============================================================================
# Plot enrichment heatmap
#==============================================================================

Heatmap(
  or_matrix,
  name = "log(OR)",
  col = mat_col,
  na_col = "white",

  cluster_rows = FALSE,
  cluster_columns = FALSE,

  row_labels = row_names_labs,
  row_names_side = "left",

  column_labels = col_names_labs,
  column_names_rot = 45,

  row_split = row_split_labs,
  row_title = NULL,

  border_gp = gpar(
    col = "black",
    lty = 1,
    lwd = 1.5
  ),

  rect_gp = gpar(
    col = "black",
    lwd = 0.2
  ),

  row_names_gp = gpar(
    fontsize = 16,
    fontface = "bold"
  ),

  column_names_gp = gpar(
    fontsize = 16,
    fontface = "bold"
  ),

  column_gap = unit(
    2,
    "mm"
  ),

  row_gap = unit(
    1.5,
    "mm"
  ),

  cell_fun = function(
      j,
      i,
      x,
      y,
      width,
      height,
      fill) {

    ci_val <- ci_matrix[i, j]
    pval_val <- pval_matrix[i, j]

    if (!is.na(ci_val)) {
      grid.text(
        ci_val,
        x,
        y + unit(2.5, "mm"),
        gp = gpar(
          fontsize = 16
        )
      )
    }

    if (!is.na(pval_val)) {
      grid.text(
        parse(
          text = pval_val
        ),
        x,
        y - unit(6.5, "mm"),
        gp = gpar(
          fontsize = 16
        )
      )
    }
  }
)

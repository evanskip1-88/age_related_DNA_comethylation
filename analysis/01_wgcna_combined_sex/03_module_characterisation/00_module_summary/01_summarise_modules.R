#==============================================================================
# Summarise selected age-related WGCNA modules
#
# Generates descriptive statistics for age-stable co-methylation modules
# (aSCMs) and age-variable co-methylation modules (aVCMs), including:
#   - number of modules
#   - total number of CpGs
#   - distribution of module sizes
#   - number of GREAT-associated genes per module
#==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

rm(list = ls())

#------------------------------------------------------------------------------
# Paths
#------------------------------------------------------------------------------

input_path <- "results/01_wgcna_combined_sex"

#------------------------------------------------------------------------------
# Load selected module information
#------------------------------------------------------------------------------

sel_mods <- read_xlsx(
  paste0(
    input_path,
    "/01a_selected_mods_summary_stats.xlsx"
  )
)

setDT(sel_mods)

# Retain age-stable and age-variable modules
sel_mods <- sel_mods[
  Category %in% c(
    "aSCMs",
    "aVCMs"
  )
]

#==============================================================================
# Module and CpG summary
#==============================================================================

mods <- unique(
  sel_mods[
    ,
    .(
      moduleColor,
      moduleNumID,
      Category,
      moduleSize
    )
  ]
)

module_summary <- mods[
  ,
  .(
    n_modules = .N,
    total_CpGs = sum(moduleSize),
    median_CpGs = median(moduleSize),
    min_CpGs = min(moduleSize),
    max_CpGs = max(moduleSize)
  ),
  by = Category
]

print(module_summary)

#==============================================================================
# Gene summary
#==============================================================================

gene_summary <- unique(
  sel_mods[
    ,
    .(
      moduleColor,
      moduleNumID,
      Category,
      gene_rgreat
    )
  ]
)

# Count unique GREAT-associated genes within each module
gene_summary[
  ,
  genes_per_module := uniqueN(
    unlist(
      strsplit(
        paste(
          gene_rgreat,
          collapse = ";"
        ),
        ";"
      )
    )
  ),
  by = .(
    moduleColor,
    moduleNumID,
    Category
  )
]

gene_summary <- unique(
  gene_summary[
    ,
    .(
      moduleColor,
      moduleNumID,
      Category,
      genes_per_module
    )
  ]
)

gene_summary_stats <- gene_summary[
  ,
  .(
    median_genes = median(genes_per_module),
    min_genes = min(genes_per_module),
    max_genes = max(genes_per_module)
  ),
  by = Category
]

print(gene_summary_stats)

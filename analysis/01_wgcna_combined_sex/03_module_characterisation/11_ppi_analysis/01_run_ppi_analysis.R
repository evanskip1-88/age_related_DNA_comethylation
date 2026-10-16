
# =============================================================================
# Purpose: Assess STRING PPI enrichment and clustering in WGCNA modules.
# Inputs:  Selected module annotations and rGREAT CpG-to-gene annotations.
# Outputs: Module-level and pooled-group PPI enrichment statistics.
# =============================================================================

library(data.table)
library(readxl)
library(writexl)
library(STRINGdb)
library(igraph)

# -------------------------------------------------------------------------
# Paths and STRING database
# -------------------------------------------------------------------------

module_file < <- "/results/01_wgcna_combined_sex/02b_wgcna_top_select_all_info.xlsx"

annotation_file <- "/data/gsDNAm/GS_sub_qc1_gene_cpgs_assoc_rGREAT_h19_all_CpGs.opi"

output_file <- "results/01_wgcna_combined_sex/02a_ppi_enrichment.xlsx"

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

string_db <- STRINGdb$new(
  version = "12.0",
  species = 9606,
  score_threshold = 400
)

# -------------------------------------------------------------------------
# Load annotations and selected WGCNA modules
# -------------------------------------------------------------------------

opi <- fread(annotation_file)
setnames(opi, c(2, 4, 7), c("chr", "gene_rgreat", "orig_gene"))

sel_mods_wgcna <- as.data.table(read_xlsx(module_file))

# Count unique annotated genes per module
gene_counts <- sel_mods_wgcna[
  ,
  .(
    num_genes = length(unique(na.omit(
      unlist(strsplit(rgreat_genes, ";", fixed = TRUE))
    ))),
    unique_genes = paste(
      unique(na.omit(unlist(strsplit(rgreat_genes, ";", fixed = TRUE)))),
      collapse = ", "
    )
  ),
  by = moduleColor
]

sel_mods_wgcna[
  ,
  `:=`(
    analysis = paste0(presStats, "_", Category),
    analysis2 = paste0("WGCNA_", Category)
  )
]

# -------------------------------------------------------------------------
# Prepare pooled WGCNA gene sets and random reference
# -------------------------------------------------------------------------

background_genes <- unique(na.omit(
  unlist(strsplit(opi$gene_rgreat, ";", fixed = TRUE))
))
background_genes <- unique(trimws(background_genes))
background_genes <- background_genes[nzchar(background_genes)]

set.seed(123)
random_genes <- sample(background_genes, 200)

wgcna_gene_dt <- sel_mods_wgcna[
  analysis2 != "WGCNA_Others",
  .(
    genes = list(unique(trimws(
      unlist(strsplit(na.omit(rgreat_genes), ";", fixed = TRUE))
    )))
  ),
  by = analysis2
]

wgcna_gene_list <- setNames(
  wgcna_gene_dt$genes,
  wgcna_gene_dt$analysis2
)

wgcna_gene_list[["Random"]] <- random_genes

# -------------------------------------------------------------------------
# Prepare individual stable and changing modules
# -------------------------------------------------------------------------

dt_genes_per_mod_sub <- merge(
  sel_mods_wgcna,
  gene_counts[num_genes >= 2],
  by = "moduleColor"
)

dt_genes_per_mod_sub <- dt_genes_per_mod_sub[
  Category != "Others" & num_genes >= 2
]

# Stable modules (aSCMs)
stable_gene_dt <- dt_genes_per_mod_sub[
  Category == "Stable",
  .(
    genes = list(unique(trimws(
      unlist(strsplit(na.omit(unique_genes), ",", fixed = TRUE))
    )))
  ),
  by = moduleColor
]

wgcna_stable_sel_mods <- setNames(
  stable_gene_dt$genes,
  stable_gene_dt$moduleColor
)

# Changing modules (aVCMs)
changing_gene_dt <- dt_genes_per_mod_sub[
  Category == "Changing",
  .(
    genes = list(unique(trimws(
      unlist(strsplit(na.omit(unique_genes), ",", fixed = TRUE))
    )))
  ),
  by = moduleColor
]

wgcna_changing_sel_mods <- setNames(
  changing_gene_dt$genes,
  changing_gene_dt$moduleColor
)

# -------------------------------------------------------------------------
# Calculate STRING PPI enrichment and global clustering coefficients
# -------------------------------------------------------------------------

compute_ppi_clustering <- function(gene_list_per_module, study_set, string_db) {

  df_res_clust_coeff <- data.frame(
    moduleColor = character(),
    analysis = character(),
    n_input_genes = integer(),
    n_mapped_genes = integer(),
    n_string_ids = integer(),
    n_edges = integer(),
    glob_clst_coeff = numeric(),
    ppi_pval = numeric(),
    stringsAsFactors = FALSE
  )

  for (mod in names(gene_list_per_module)) {

    cat("PPI for module:", mod, "\n")

    sel_gene_lst <- unique(trimws(gene_list_per_module[[mod]]))
    sel_gene_lst <- sel_gene_lst[
      !is.na(sel_gene_lst) & nzchar(sel_gene_lst)
    ]

    df_unique_genes <- data.frame(
      gene_rgreat = sel_gene_lst,
      stringsAsFactors = FALSE
    )

    # Map genes to STRING proteins
    genes_mapped <- as.data.table(
      string_db$map(
        df_unique_genes,
        "gene_rgreat",
        removeUnmappedRows = TRUE
      )
    )

    if (nrow(genes_mapped) == 0) {

      df_res_clust_coeff <- rbind(
        df_res_clust_coeff,
        data.frame(
          moduleColor = mod,
          analysis = study_set,
          n_input_genes = length(sel_gene_lst),
          n_mapped_genes = 0,
          n_string_ids = 0,
          n_edges = 0,
          glob_clst_coeff = NA_real_,
          ppi_pval = NA_real_
        )
      )

      next
    }

    string_ids <- unique(genes_mapped$STRING_id)

    # STRING PPI enrichment
    ppi_pval <- tryCatch(
      string_db$get_ppi_enrichment(string_ids),
      error = function(e) NA_real_
    )

    # Retrieve interactions among mapped proteins
    ppi <- string_db$get_interactions(string_ids)

    if (is.null(ppi) || nrow(ppi) == 0) {

      df_res_clust_coeff <- rbind(
        df_res_clust_coeff,
        data.frame(
          moduleColor = mod,
          analysis = study_set,
          n_input_genes = length(sel_gene_lst),
          n_mapped_genes = nrow(genes_mapped),
          n_string_ids = length(string_ids),
          n_edges = 0,
          glob_clst_coeff = 0,
          ppi_pval = ppi_pval
        )
      )

      next
    }

    # Global clustering coefficient
    ppi_graph <- graph_from_data_frame(ppi, directed = FALSE)

    glob_clst_coeff <- transitivity(
      ppi_graph,
      type = "global",
      isolates = "zero"
    )

    df_res_clust_coeff <- rbind(
      df_res_clust_coeff,
      data.frame(
        moduleColor = mod,
        analysis = study_set,
        n_input_genes = length(sel_gene_lst),
        n_mapped_genes = nrow(genes_mapped),
        n_string_ids = length(string_ids),
        n_edges = nrow(ppi),
        glob_clst_coeff = glob_clst_coeff,
        ppi_pval = ppi_pval
      )
    )
  }

  return(df_res_clust_coeff)
}

# -------------------------------------------------------------------------
# Run PPI analyses
# -------------------------------------------------------------------------

# Stable modules
res_cluster_coeffs_wgcna_stable <- compute_ppi_clustering(
  gene_list_per_module = wgcna_stable_sel_mods,
  study_set = "aSCMs",
  string_db = string_db
)

# Changing modules
res_cluster_coeffs_wgcna_variable <- compute_ppi_clustering(
  gene_list_per_module = wgcna_changing_sel_mods,
  study_set = "aVCMs",
  string_db = string_db
)

# Pooled WGCNA groups and random reference
res_cluster_coeffs_combined_modules <- compute_ppi_clustering(
  gene_list_per_module = wgcna_gene_list,
  study_set = "combined_mods",
  string_db = string_db
)

# -------------------------------------------------------------------------
# Combine module-level and pooled-group results
# -------------------------------------------------------------------------

res_ppi <- rbind(
  res_cluster_coeffs_wgcna_stable,
  res_cluster_coeffs_wgcna_variable
)

res_ppi <- merge(
  unique(sel_mods_wgcna[, .(moduleColor, moduleNumID)]),
  res_ppi,
  by = "moduleColor"
)

# Preserve original output structure for pooled groups
setnames(
  res_cluster_coeffs_combined_modules,
  "moduleColor",
  "moduleNumID"
)

res_cluster_coeffs_combined_modules <- cbind(
  moduleColor = "NA",
  res_cluster_coeffs_combined_modules
)

res_ppi_c <- rbind(
  res_ppi,
  res_cluster_coeffs_combined_modules
)

# -------------------------------------------------------------------------
# Save results
# -------------------------------------------------------------------------

write_xlsx(res_ppi_c, output_file)


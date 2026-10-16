
# =============================================================================
# Purpose: Compare global PPI clustering coefficients across gene groups.
# Inputs:  WGCNA gene lists and background rGREAT gene annotations.
# Outputs: Supplementary Figure 11c and pairwise Wilcoxon test results.
# =============================================================================

library(data.table)
library(STRINGdb)
library(igraph)
library(ggplot2)
library(ggpubr)

# -------------------------------------------------------------------------
# Paths and settings
# -------------------------------------------------------------------------

gene_list_file <- "/results/01_wgcna_combined_sex/heteroskid_wgcna_gene_list.rds"

annotation_file <- "data/gsDNAm/GS_sub_qc1_gene_cpgs_assoc_rGREAT_h19_all_CpGs.opi"

outputFigs_dir <- "/figures"
output_tables_dir <- "/results/01_wgcna_combined_sex"
figure_file <- file.path(outputFigs_dir, "Figure_S11c_PPI_clustering.pdf")
stats_file <- file.path(output_tables_dir, "Figure_S11c_pairwise_tests.csv")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(123)

n_bootstrap <- 50
sample_size <- 100

# -------------------------------------------------------------------------
# Load gene lists
# -------------------------------------------------------------------------

opi <- fread(annotation_file)

heteroskid_wgcna_gene_list <- readRDS(gene_list_file)

whole_genes <- unique(na.omit(
  unlist(strsplit(opi$annotated_genes, ";", fixed = TRUE))
))

whole_genes <- trimws(whole_genes)
whole_genes <- whole_genes[nzchar(whole_genes)]

stable_gene_list <- heteroskid_wgcna_gene_list[["Heteroskid_Stable"]]
changing_gene_list <- heteroskid_wgcna_gene_list[["Heteroskid_Changing"]]

stopifnot(
  length(stable_gene_list) > 0,
  length(changing_gene_list) > 0,
  length(whole_genes) > 0
)

# Original background sampling procedure
random_sample_gene_list <- sample(
  whole_genes,
  length(changing_gene_list),
  replace = TRUE
)

# -------------------------------------------------------------------------
# Initialise STRING
# -------------------------------------------------------------------------

string_db <- STRINGdb$new(
  version = "12.0",
  species = 9606,
  score_threshold = 0.7
)

# -------------------------------------------------------------------------
# Calculate global PPI clustering coefficient
# -------------------------------------------------------------------------

compute_ppi_clustering <- function(gene_list, study_set) {

  results <- vector("list", length(gene_list))

  for (i in seq_along(gene_list)) {

    mod <- names(gene_list)[i]
    cat("PPI for gene list:", mod, "\n")

    sel_gene_lst <- unique(gene_list[[i]])
    sel_gene_lst <- sel_gene_lst[
      !is.na(sel_gene_lst) & nzchar(sel_gene_lst)
    ]

    df_unique_genes <- data.frame(
      gene_rgreat = sel_gene_lst
    )

    genes_mapped <- as.data.table(
      string_db$map(
        df_unique_genes,
        "gene_rgreat",
        removeUnmappedRows = TRUE
      )
    )

    glob_clst_coeff <- NA_real_

    if (nrow(genes_mapped) > 0) {

      string_ids <- genes_mapped$STRING_id
      ppi <- string_db$get_interactions(string_ids)

      if (!is.null(ppi) && nrow(ppi) > 0) {

        ppi_graph <- graph_from_data_frame(
          ppi,
          directed = FALSE
        )

        glob_clst_coeff <- transitivity(
          ppi_graph,
          type = "global",
          isolates = "zero"
        )
      }
    }

    results[[i]] <- data.table(
      moduleColor = mod,
      analysis = study_set,
      glob_clst_coeff = glob_clst_coeff
    )
  }

  rbindlist(results)
}

# -------------------------------------------------------------------------
# Resampling analysis
# -------------------------------------------------------------------------

res_stable <- vector("list", n_bootstrap)
res_changing <- vector("list", n_bootstrap)
res_random <- vector("list", n_bootstrap)

for (i in seq_len(n_bootstrap)) {

  cat("Bootstrap:", i, "\n")

  # Stable modules
  stable_sample <- sample(
    stable_gene_list,
    sample_size,
    replace = TRUE
  )

  res_stable[[i]] <- compute_ppi_clustering(
    gene_list = list(aSCPs = stable_sample),
    study_set = "Stable"
  )

  # Changing modules
  changing_sample <- sample(
    changing_gene_list,
    sample_size,
    replace = TRUE
  )

  res_changing[[i]] <- compute_ppi_clustering(
    gene_list = list(aVCPs = changing_sample),
    study_set = "Changing"
  )

  # Random background
  random_sample <- sample(
    random_sample_gene_list,
    sample_size,
    replace = TRUE
  )

  res_random[[i]] <- compute_ppi_clustering(
    gene_list = list(Random = random_sample),
    study_set = "Random"
  )
}

# -------------------------------------------------------------------------
# Combine results
# -------------------------------------------------------------------------

df_res <- rbindlist(
  c(res_stable, res_changing, res_random),
  use.names = TRUE
)

df_res[, study_set := factor(
  rep(c("aSCPs", "aVCPs", "Random"), each = n_bootstrap),
  levels = c("aSCPs", "aVCPs", "Random")
)]

# -------------------------------------------------------------------------
# Statistical comparisons
# -------------------------------------------------------------------------

pairwise_results <- pairwise.wilcox.test(
  df_res$glob_clst_coeff,
  df_res$study_set,
  p.adjust.method = "bonferroni"
)

print(pairwise_results)

comparisons <- list(
  c("aSCPs", "aVCPs"),
  c("aSCPs", "Random"),
  c("aVCPs", "Random")
)

# -------------------------------------------------------------------------
# Supplementary Figure 11c
# -------------------------------------------------------------------------

p <- ggplot(
  df_res,
  aes(
    x = study_set,
    y = glob_clst_coeff,
    fill = study_set
  )
) +
  geom_boxplot(
    width = 0.3,
    outlier.shape = NA,
    alpha = 0.6
  ) +
  geom_jitter(
    width = 0.1,
    alpha = 0.5
  ) +
  theme_bw() +
  labs(
    x = "",
    y = "Global clustering coefficient"
  ) +
  scale_fill_manual(
    values = c(
      "aSCPs" = "steelblue",
      "aVCPs" = "firebrick",
      "Random" = "gray"
    )
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold"
    ),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.text.y = element_text(
      size = 16,
      hjust = 0,
      face = "bold"
    )
  ) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    p.adjust.method = "bonferroni",
    label = "p.signif"
  )

print(p)

# -------------------------------------------------------------------------
# Save results
# -------------------------------------------------------------------------

ggsave(
  filename = figure_file,
  plot = p,
  width = 7,
  height = 6,
  units = "in"
)

fwrite(df_res, file.path(output_tables_dir, "Figure_S11c_clustering_data.csv"))

pairwise_table <- as.data.table(
  as.table(pairwise_results$p.value)
)

setnames(
  pairwise_table,
  c("Group1", "Group2", "adjusted_p_value")
)

pairwise_table <- pairwise_table[!is.na(adjusted_p_value)]

fwrite(pairwise_table, stats_file)



# 01_run_aging_gene_enrichment.R
# Gene-set enrichment used for GTEx (Fig. S17), curated aging genes,
# and longevity interventions (Fig. 6b / supplementary results).

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(dplyr)
})

whole_genes <- 17593L # Total number of genes used as background set
input_gene_sets <- "/results/01_wgcna_combined_sex/heteroskid_wgcna_gene_list.rds"
wangetal_dir <- "/data/aging_genes/wang_et_al_2024"
victoria_cr_path <- "/data/aging_genes/victoria_et_al"
sun_data_path <- "/data/aging_genes/sun_et_al"
yang_et_al_path <- "/data/aging_genes/Yang_et_al_Metformin"
hagr_path <- "/data/aging_genes"
atlas_genes_path <- "/data/aging_genes/AgingAtlas"
res_dir <- "/results/01_wgcna_combined_sex"
output_dir <- "/results/trait_associations/aging_gene_enrichment"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

heteroskid_wgcna_gene_list <- readRDS(input_gene_sets)
stopifnot(is.list(heteroskid_wgcna_gene_list),
          !is.null(names(heteroskid_wgcna_gene_list)))

# Function to perform enrichment analysis
perform_enrichment_fun <- function(study_hit_list, pub_hit_list, bg_set) {
  res_list <- list()
  for (studySet in names(study_hit_list)) {
    study_hits <- study_hit_list[[studySet]]
    pub_hits <- pub_hit_list
    overlaps <- intersect(study_hits, pub_hits)
    # Hit counts
    N <- bg_set  # Total number of genes in background
    K <- length(pub_hits)  # published hit list
    n <- length(study_hits)  # Study hit list
    k <- length(overlaps) # Overlapping hit list
    # Avoid division by zero in OR calculation
    if (k == 0 || K == k || n == k || (N - K - (n - k)) == 0) {
      OR <- NA
      SE_log_OR <- NA
      CI_lower <- NA
      CI_upper <- NA
    } else {
      OR <- (k / (K - k)) / ((n - k) / (N - K - (n - k)))
      SE_log_OR <- sqrt(1/k + 1/(n-k) + 1/(K-k) + 1/(N-K-(n-k)))
      Z_alpha <- 1.96  # 95% confidence interval Z-score
      log_OR <- log(OR)
      CI_lower <- exp(log_OR - Z_alpha * SE_log_OR)
      CI_upper <- exp(log_OR + Z_alpha * SE_log_OR)
    }
    res <- data.frame(
      analysis = studySet,
      publishedHits = K,
      study_hits = n,
      overlap = k,
      pvalue = phyper(k - 1, K, N - K, n, lower.tail = FALSE), 
      OR = OR,
      CI_lower = CI_lower,
      CI_upper = CI_upper
    )
    res_list[[studySet]] <- res
  }
  return(do.call(rbind, res_list))
}

run_enrich_fun <-function(pub_data, study_hits) {
  # Unique tissues
  tissue_list <- unique(pub_data$tissueName)
  sex_class <- unique(pub_data$Class)
  result_list <- list()
  # Iterate over all tissues
  for (tissue_id in tissue_list) {
    # Subset DEGs for the current tissue
    tissue <- pub_data[tissueName == tissue_id, ]
    # Iterate over sex classes
    sex_class <- unique(tissue$Class)
    for (sex in sex_class) {
      # Subset DEGs for the current sex class
      sex_tissue <- tissue[Class == sex, ]
      sex_tissue_gene_lst <- unique(sex_tissue$symbol)
      # Loop over each gene set and perform enrich
        res_enrich <- cbind(tissue = tissue_id,
                      sex_class = sex,
                      perform_enrichment_fun(
                      study_hit_list = study_hits,
                      pub_hit_list = sex_tissue_gene_lst,
                      bg_set = whole_genes ) )
        result_list[[paste0("Tissue_", tissue_id, "_", sex)]] <- res_enrich
      }
  }
  # Print summary of enrichment results
  res_enrichment <- setDT(do.call('rbind', result_list))
  return(res_enrichment)
}

# 1. GTEx age-related sex-stratified DEGs (Wang et al.; Fig. S17)
wang_file <- file.path(wangetal_dir, "TableS3_Sex_stratified_age_DGEs_DASs_list.xlsx")
wang_tissue_labs <- setDT(read_xlsx(wang_file, sheet = 1, skip = 1))
colnames(wang_tissue_labs) <- c("tissue", "tissueName")
wang_tissue_labs[, tissue := as.numeric(tissue)]
wang_DEGs <- merge(wang_tissue_labs,
                   setDT(read_xlsx(wang_file, sheet = 2, skip = 1)), by = "tissue")
res_DEGs_enrichment <- run_enrich_fun(wang_DEGs, heteroskid_wgcna_gene_list)
res_DEGs_enrichment[, Group := "DEGs"]
writexl::write_xlsx(res_DEGs_enrichment,
                    file.path(output_dir, "10a_GTEx_age_DEGs_enrichment.xlsx"))
# Only DEGs are exported: DAS was not used for submitted Fig. S17.

# 2. Curated aging genes
combine_files <- function(dataPath){
  lst.files <- list.files(path = file.path(atlas_genes_path, dataPath), full.names = TRUE)
  dt.lst.files <- list()
  for(f in lst.files){
    dt <- fread(f)
    dt$Gene_ID <- as.integer(dt$Gene_ID)
    dt.lst.files[[f]] <- dt
  }
  dt.lst.files.c <- bind_rows(dt.lst.files)
  return(dt.lst.files.c)
}
dt_humanGenes <- combine_files("humans")
humanGenes_list <- unique(dt_humanGenes$Symbol)

dt_modelGenes <- combine_files("musculus")
dt_modelGenes$Symbol <- toupper(dt_modelGenes$Symbol)
dt_modelGenes[Symbol %in% dt_humanGenes$Symbol, ]
modelGenes_list <- unique(dt_modelGenes$Symbol)

# Combine above
dt_all_atlas_genes <- bind_rows(dt_humanGenes, dt_modelGenes)
dt_all_atlas_genes_list <- unique(dt_all_atlas_genes$Symbol)
res_aging_atlas_gene_enrichment <- perform_enrichment_fun(study_hit_list = heteroskid_wgcna_gene_list, 
                                                          pub_hit_list = dt_all_atlas_genes_list, 
                                                          bg_set = whole_genes)
res_aging_atlas_gene_enrichment$Dataset <- "AgingAtlas"
res_aging_atlas_gene_enrichment$Group <- "DEGs"
genAgeHuman_path = file.path(hagr_path, "HAGR", "datasetGenAge", "human_genes")
dt_genAgeHuman  <- fread(paste0(genAgeHuman_path, "/genage_human.csv"))
digAge_Atlas_genes_list <- unique(dt_genAgeHuman$symbol)
head(digAge_Atlas_genes_list)

res_genAgeHuman_enrichment <- perform_enrichment_fun(study_hit_list = heteroskid_wgcna_gene_list, 
                                                     pub_hit_list = digAge_Atlas_genes_list, 
                                                     bg_set = whole_genes)
res_genAgeHuman_enrichment$Dataset <- "HAGR_genes"
res_genAgeHuman_enrichment$Group <- "DEGs"

# b) HAGR: Under- & Over- expressed genes
over_expr_genes <- setDT(read_xlsx(file.path(hagr_path, "HAGR", "datasetGeneExpr", "ageing_signatures.xlsx"), 
                                   sheet = "TableS3-Over.All", skip = 1))
over_expr_genes$Group <- "Over_expressed"

under_expr_genes <- setDT(read_xlsx(file.path(hagr_path, "HAGR", "datasetGeneExpr", "ageing_signatures.xlsx"), 
                                    sheet = "TableS7-Under.All", skip = 1))
under_expr_genes$Group <- "Under_expressed"
# combine
dt_all_expr <- bind_rows(over_expr_genes, under_expr_genes)
HAGR_all_expr_genes_list <- unique(dt_all_expr$Gene)

res_HAGR_twas_enrichment <- perform_enrichment_fun(study_hit_list = heteroskid_wgcna_gene_list, 
                                                   pub_hit_list = HAGR_all_expr_genes_list, 
                                                   bg_set = whole_genes)
res_HAGR_twas_enrichment$Dataset <- "HAGR_meta"
res_HAGR_twas_enrichment$Group <- "HAGR_DEGs"
twas_all_genes_list <-  unique( c(digAge_Atlas_genes_list, dt_all_atlas_genes_list, HAGR_all_expr_genes_list) )

res_all_twas_enrichment <- perform_enrichment_fun(study_hit_list = heteroskid_wgcna_gene_list, 
                                                  pub_hit_list = twas_all_genes_list, 
                                                  bg_set = whole_genes)
res_all_twas_enrichment$Dataset <- "curated-Aging-Genes"
res_all_twas_enrichment$Group <- "DEGs"

# Preserve separate sources for the exploratory curated-gene point plot.
res_curated_sources <- rbindlist(list(res_aging_atlas_gene_enrichment,
                                      res_genAgeHuman_enrichment,
                                      res_HAGR_twas_enrichment), fill = TRUE)
writexl::write_xlsx(res_curated_sources,
                    file.path(output_dir, "10c_curated_genes_enrich.xlsx"))
writexl::write_xlsx(as.data.frame(res_all_twas_enrichment),
                    file.path(output_dir, "10c_curated_aging_genes_combined.xlsx"))

# 3. Longevity interventions
# Caloric restriction only
victoria_et_al_cr <- setDT(read_xlsx( paste0(victoria_cr_path, "/NIHMS1820506-supplement-Supplemental_Data_table_2.xlsx"), skip = 5, sheet = 3))
victoria_et_al_cr_sub <- victoria_et_al_cr[, c(1:6)]
colnames(victoria_et_al_cr_sub) <- c("GeneID", "baseMean", "foldChange", "lfcSE", "padj", "diff_exp")
victoria_et_al_cr_sub[, Gene := toupper(sub(".*_", "", GeneID))]  # Create Gene column
setcolorder(victoria_et_al_cr_sub, c("Gene", setdiff(names(victoria_et_al_cr_sub), c("Gene", "GeneID"))))  # Reorder
victoria_et_al_cr_sub[, GeneID := NULL]  # Remove GeneID
cr_victoria_et_al_sub_genes <- victoria_et_al_cr_sub$Gene
res_victoria_cr_enrichment <- perform_enrichment_fun(study_hit_list = heteroskid_wgcna_gene_list, 
                                                     pub_hit_list = cr_victoria_et_al_sub_genes, 
                                                     bg_set = whole_genes)
res_victoria_cr_enrichment$Dataset <- "Victoria_et_al"
res_victoria_cr_enrichment$Group <- "DR"
# fasting
victoria_et_al_fast <- setDT(read_xlsx( paste0(victoria_cr_path, "/NIHMS1820506-supplement-Supplemental_Data_table_2.xlsx"), skip = 5, sheet = 4))
victoria_et_al_fast_sub <- victoria_et_al_fast[, c(1:6)]
colnames(victoria_et_al_fast_sub) <- c("GeneID", "baseMean", "foldChange", "lfcSE", "padj", "diff_exp")
victoria_et_al_fast_sub[, Gene := toupper(sub(".*_", "", GeneID))]  # Create Gene column
setcolorder(victoria_et_al_fast_sub, c("Gene", setdiff(names(victoria_et_al_fast_sub), c("Gene", "GeneID"))))  # Reorder
victoria_et_al_fast_sub[, GeneID := NULL]  # Remove GeneID
fast_victoria_et_al_sub_genes <- victoria_et_al_fast_sub$Gene
res_victoria_fast_enrichment <- perform_enrichment_fun(study_hit_list = heteroskid_wgcna_gene_list, 
                                                       pub_hit_list = fast_victoria_et_al_sub_genes, 
                                                       bg_set = whole_genes)
res_victoria_fast_enrichment$Dataset <- "Victoria_et_al"
res_victoria_fast_enrichment$Group <- "Fasting"
# Timing
victoria_et_al_time <- setDT(read_xlsx( paste0(victoria_cr_path, "/NIHMS1820506-supplement-Supplemental_Data_table_2.xlsx"), skip = 5, sheet = 5))
victoria_et_al_time$GeneID...28 <- toupper(victoria_et_al_time$GeneID...28)
victoria_et_al_time_sub_genes <- victoria_et_al_time$GeneID...28
res_victoria_time_enrichment <- perform_enrichment_fun(study_hit_list = heteroskid_wgcna_gene_list, 
                                                       pub_hit_list = victoria_et_al_time_sub_genes, 
                                                       bg_set = whole_genes)
res_victoria_time_enrichment$Dataset <- "Victoria_et_al"
res_victoria_time_enrichment$Group <- "Timing"
sun_et_al_2013 <- setDT(read_xlsx( paste0(sun_data_path, "/elife-01098-supp1-v1.xlsx"), skip = 3, sheet = 1))
sun_et_al_2013$`Gene Symbol` <- toupper(sun_et_al_2013$`Gene Symbol`)
sun_et_al_2013_sub_genes <- sun_et_al_2013$`Gene Symbol`
res_sun_ghr_enrichment <- perform_enrichment_fun(study_hit_list = heteroskid_wgcna_gene_list, 
                                                 pub_hit_list = sun_et_al_2013_sub_genes, 
                                                 bg_set = whole_genes)
res_sun_ghr_enrichment$Dataset <- "Sun_et_al"
res_sun_ghr_enrichment$Group <- "GHR_knock_out"
sheets <- "Rescue (Brain)"
res_rescue_list <- list()
for(sheet in sheets) {
  print(sheet)
  dt_yang_rescue <-  setDT(read_xlsx(  paste0(yang_et_al_path,"/Table_S4_mmc2.xlsx"), sheet = sheet, skip = 1))
  colnames(dt_yang_rescue)[1:3] <- c("Gene", "Cell_type", "Cluster")
  yang_rescue_gene_list <- dt_yang_rescue$Gene
  res_yang_rescue_enrichment <- perform_enrichment_fun(study_hit_list = heteroskid_wgcna_gene_list, 
                                                            pub_hit_list = yang_rescue_gene_list, 
                                                            bg_set = whole_genes)
  res_yang_rescue_enrichment$Tissue <- sheet
  res_rescue_list[[sheet]] <-  res_yang_rescue_enrichment
}
res_yang_rescue_enrichment_fnl <- setDT(do.call('rbind', res_rescue_list))
res_yang_rescue_enrichment_fnl$Dataset <- "Yang_et_al_2024"
res_yang_rescue_enrichment_fnl$Group <- "Metformin"
# Select only brain for repressing interven
res_yang_rescue_enrichment_fnl_brain <- res_yang_rescue_enrichment_fnl[Tissue == "Rescue (Brain)", ]
res_yang_rescue_enrichment_fnl_brain[, Tissue:=NULL]
# a) Combine all interventions 
res_all_intervent <- rbind(res_victoria_time_enrichment, res_victoria_cr_enrichment, res_victoria_fast_enrichment, 
                           res_sun_ghr_enrichment, res_yang_rescue_enrichment_fnl_brain)

writexl::write_xlsx(as.data.frame(res_all_intervent),
                    file.path(output_dir, "10b_intervention_enrich.xlsx"))


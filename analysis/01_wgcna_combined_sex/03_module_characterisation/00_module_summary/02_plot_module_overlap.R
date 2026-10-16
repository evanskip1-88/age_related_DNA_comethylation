
# =============================================================================
# Purpose: Visualise CpG and gene overlaps across WGCNA and pairwise groups.
# Inputs:  WGCNA/heteroskedasticity CpG and gene lists.
# Outputs: Supplementary Figures S8e–f.
# =============================================================================

library(VennDiagram)
library(grid)

# -------------------------------------------------------------------------
# Input and output paths
# -------------------------------------------------------------------------

input_dir <- "results/01_wgcna_combined_sex"

cpg_file <- file.path(
  input_dir,
  "heteroskid_wgcna_cpg_list.rds"
)

gene_file <- file.path(
  input_dir,
  "heteroskid_wgcna_gene_list.rds"
)

output_dir <- "figures"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Load CpG and gene lists
# -------------------------------------------------------------------------

heteroskid_cpg_list <- readRDS(cpg_file)
heteroskid_gene_list <- readRDS(gene_file)

# Explicitly select the four comparison groups
group_names <- c(
  "WGCNA_Stable",
  "WGCNA_Changing",
  "Heteroskid_Stable",
  "Heteroskid_Changing"
)

stopifnot(
  all(group_names %in% names(heteroskid_cpg_list)),
  all(group_names %in% names(heteroskid_gene_list))
)

cpg_list <- heteroskid_cpg_list[group_names]
gene_list <- heteroskid_gene_list[group_names]

# Figure labels
plot_labels <- c(
  "aSCMs",
  "aVCMs",
  "aSCPs",
  "aVCPs"
)

names(cpg_list) <- plot_labels
names(gene_list) <- plot_labels

# -------------------------------------------------------------------------
# Venn diagram function
# -------------------------------------------------------------------------

venn_function <- function(input_list) {

  venn.diagram(
    x = input_list,
    filename = NULL,
    output = TRUE,
    imagetype = "png",
    col = "black",
    fill = c("red", "blue", "green", "purple"),
    alpha = 0.5,
    cex = 1.5,
    fontface = "bold",
    cat.col = c("red", "blue", "green", "purple"),
    cat.cex = 1.2,
    cat.fontface = "bold",
    main = "",
    main.cex = 2,
    margin = 0.11
  )
}

# -------------------------------------------------------------------------
# Supplementary Figure S8e: CpG overlap
# -------------------------------------------------------------------------

pdf(
  file.path(output_dir, "Figure_S8e_CpG_overlap.pdf"),
  width = 8,
  height = 8
)

grid.newpage()
grid.draw(venn_function(cpg_list))

dev.off()

# -------------------------------------------------------------------------
# Supplementary Figure S8f: Gene overlap
# -------------------------------------------------------------------------

pdf(
  file.path(output_dir, "Figure_S8f_gene_overlap.pdf"),
  width = 8,
  height = 8
)

grid.newpage()
grid.draw(venn_function(gene_list))

dev.off()

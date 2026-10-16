# Associate Generation Scotland CpGs with NCBI build 37 genes.
#
# Original analysis:
# CpG position between gene start - 5 kb and gene end + 1 kb.
# This is a coordinate-based, non-strand-aware window.
#
# Note: Despite the original script's 50 kb comment,
# the implemented offsets were 5 kb and 1 kb.

library(data.table)
library(dplyr)

# Input files
gene_file <- "/data/gsDNAm/NCBI37.3.gene.loc"
opi_file <- "data/gsDNAm/GS_sub_qc1_h19_gene_cpgs_assoc_rGREAT.opi"

# Output
output_dir <- "data/gwas/magma"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load NCBI build 37 gene coordinates
build_37 <- fread(gene_file)

setnames(
  build_37,
  c("geneID", "chr", "bp_start", "bp_end", "strand", "gene")
)

build_37[, strand := NULL]

# Load Generation Scotland CpG annotations
opi <- fread(opi_file)

setnames(
  opi,
  old = names(opi)[1:3],
  new = c("cpg", "bp", "chr")
)

opi[, chr := sub("^chr", "", chr)]

# Associate CpGs with genes using the original window
dt_res <- opi %>%
  left_join(
    build_37,
    by = "chr",
    relationship = "many-to-many"
  ) %>%
  filter(
    bp >= bp_start - 5000,
    bp <= bp_end + 1000
  ) %>%
  select(-bp_start, -bp_end)

# Preserve original output column names
colnames(dt_res) <- c(
  "cpg", "bp", "chr", "start", "end", "width",
  "strand", "gene_opi", "dist_to_TSS",
  "geneID", "gene_build37"
)

# Export
fwrite(
  dt_res,
  file.path(
    output_dir,
    "opi_build37_cpg5kb_up_1kdown_window_merged_genes.txt"
  ),
  sep = "\t",
  quote = FALSE
)

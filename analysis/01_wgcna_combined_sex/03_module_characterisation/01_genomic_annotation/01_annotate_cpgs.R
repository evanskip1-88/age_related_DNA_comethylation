#==============================================================================
#
# Annotate CpG sites using annotatr
#
# CpG sites passing DNAm quality control are annotated with:
#   1. CpG context annotations
#   2. Gene-based annotations
#
# Genome build: hg19
#
#==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(annotatr)
  library(GenomicRanges)
})

rm(list = ls())

# Set paths
opi_path <- "data/gsDNAm"

#------------------------------------------------------------------------------
# Load CpG positions
#------------------------------------------------------------------------------

opi <- fread(
  paste0(opi_path, "/GS_sub_qc1_h19.bed")
)

colnames(opi) <- c(
  "chr",
  "position",
  "end",
  "probe_id",
  "gene",
  "strand"
)

# Convert CpG positions to GRanges
opi_gr <- GRanges(
  seqnames = opi$chr,
  ranges = IRanges(
    start = opi$position,
    end = opi$position
  ),
  probe_id = opi$probe_id
)

#------------------------------------------------------------------------------
# Build annotatr annotations
#------------------------------------------------------------------------------

# CpG-context annotations
annotations_cpgs <- build_annotations(
  genome = "hg19",
  annotations = "hg19_cpgs"
)

# Gene-based annotations
annotations_genic <- build_annotations(
  genome = "hg19",
  annotations = c(
    "hg19_basicgenes",
    "hg19_genes_intergenic"
  )
)

#------------------------------------------------------------------------------
# Annotate CpG context
#------------------------------------------------------------------------------

dnam_annotated_cpgs <- annotate_regions(
  regions = opi_gr,
  annotations = annotations_cpgs,
  ignore.strand = TRUE,
  quiet = FALSE
)

dt_dnam_annotated_cpgs <- as.data.table(
  dnam_annotated_cpgs
)

#------------------------------------------------------------------------------
# Annotate gene context
#------------------------------------------------------------------------------

dnam_annotated_genic <- annotate_regions(
  regions = opi_gr,
  annotations = annotations_genic,
  ignore.strand = TRUE,
  quiet = FALSE
)

dt_dnam_annotated_genic <- as.data.table(
  dnam_annotated_genic
)

# Retain relevant annotation fields
dt_annotated_genic_sub <- unique(
  dt_dnam_annotated_genic[
    ,
    .(
      probe_id,
      seqnames,
      start,
      annot.type,
      annot.symbol
    )
  ]
)

#------------------------------------------------------------------------------
# Save annotation tables
#------------------------------------------------------------------------------

fwrite(
  dt_annotated_genic_sub,
  paste0(opi_path, "/GS_sub_qc1_h19_genic.annotr"),
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t"
)

fwrite(
  dt_dnam_annotated_cpgs,
  paste0(opi_path, "/GS_sub_qc1_h19_cpgs.annotr"),
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t"
)

cat("CpG annotation complete.\n")

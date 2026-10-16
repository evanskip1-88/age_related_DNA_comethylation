#==============================================================================
# Prepare Roadmap Epigenomics histone-mark metadata
#
# Matches histone-mark peak filenames to tissue metadata and selects
# tissues used for downstream enrichment analysis.
#==============================================================================

library(data.table)

#------------------------------------------------------------------------------
# Paths
#------------------------------------------------------------------------------

input_path <- "/data/RoadMapHistoneMark"
output_path <- "/data/RoadMapHistoneMark"

#------------------------------------------------------------------------------
# Download histone-mark filename from RoadMap (NB: Not and Rscript; Run separately)
#------------------------------------------------------------------------------
#wget -qO- \
    #"http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/gappedPeak/" | \
    #grep -Eo 'E[0-9]+-(H3K4me1|H3K4me3|H3K9ac|H3K9me3|H3K27ac|H3K27me3|H3K36me3)\.gappedPeak\.gz' \
    #| sort | uniq > ${input_path}/histone_marks.txt

#------------------------------------------------------------------------------
# Load histone-mark filenames
#------------------------------------------------------------------------------

hm <- fread(
  file.path(input_path, "histone_marks.txt"),
  header = FALSE
)

setnames(hm, "V1", "file_name")

hm[
  ,
  histone_mark := sub(
    "^E[0-9]+-(.*)\\.gappedPeak\\.gz$",
    "\\1",
    file_name
  )
]

hm[
  ,
  tissue_id := sub(
    "^(E[0-9]+)-.*",
    "\\1",
    file_name
  )
]

#------------------------------------------------------------------------------
# Load Roadmap tissue metadata
#------------------------------------------------------------------------------

hm_meta <- fread(
  file.path(input_path, "histone_meta_data.txt")
)

hm_meta_sub <- hm_meta[
  ,
  .(
    histone_mark = MARK,
    tissue_id = EID,
    e_mnemonic = `E-Mnemonic`,
    standard_epigenome = `Standardised epigenome name`,
    edacc_epigenome_name =
      `Epigenome name (from EDACC Release 9 directory)`,
    file_name = FNAME
  )
]

hm_meta_sub[
  ,
  file_name := paste0(
    sub("\\*", "", file_name),
    "gappedPeak.gz"
  )
]

#------------------------------------------------------------------------------
# Match files to tissue metadata
#------------------------------------------------------------------------------

hm_meta <- merge(
  hm,
  hm_meta_sub[
    ,
    .(
      file_name,
      e_mnemonic,
      standard_epigenome,
      edacc_epigenome_name
    )
  ],
  by = "file_name"
)

#------------------------------------------------------------------------------
# Select tissues
#------------------------------------------------------------------------------

selected_tissues <- c(
  "E035", "E055", "E062", "E063",
  "E065", "E066", "E071", "E073",
  "E079", "E095", "E096", "E107",
  "E113", "E114", "E118", "E119"
)

hm_selected <- hm_meta[
  tissue_id %in% selected_tissues
]

#------------------------------------------------------------------------------
# Save selected metadata
#------------------------------------------------------------------------------

dir.create(
  output_path,
  recursive = TRUE,
  showWarnings = FALSE
)

fwrite(
  hm_selected,
  file.path(
    output_path,
    "histone_marks_final_meta_sel_tiss.txt"
  ),
  sep = "\t",
  quote = FALSE
)

cat("Selected histone-mark files:", nrow(hm_selected), "\n")
cat("Metadata preparation completed.\n")

# Convert OpenGWAS summary statistics to MAGMA SNP p-value input.

library(data.table)

gwas_id <- "ieu-b-5067"

input_file <- file.path(
  "data/gwas/sumstats",
  paste0(gwas_id, ".vcf.gz.sub")
)

output_dir <- "/data/gwas/magma"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load extracted VCF summary-statistics table
dt_gwas <- fread(input_file)

# Extract fields from the GWAS-specific column
dt_gwas[, c("ES", "SE", "LP", "AF", "ID") :=
          tstrsplit(get(gwas_id), ":", fixed = TRUE)]

# Convert -log10(p) to p-value
dt_gwas[, LP := as.numeric(LP)]
dt_gwas[, P := 10^(-LP)]

# Prepare MAGMA input
dt_pval <- dt_gwas[, .(SNP = ID, P)]

# Save
output_file <- file.path(
  output_dir,
  paste0(gwas_id, ".pval")
)

fwrite(
  dt_pval,
  output_file,
  sep = "\t",
  quote = FALSE
)

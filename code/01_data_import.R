
# ==========================================
# AMR Genotype Analysis in E. coli
# Script 01: Data Import
# ==========================================

# Load required packages
library(readr)

# Import cleaned E. coli dataset
isolates <- readRDS("data/ecoli_clean.rds")

# Import NCBI reference gene catalogue
reference_genes <- read.delim(
  "data/ReferenceGeneCatalog.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Display dataset dimensions
dim(isolates)
dim(reference_genes)

# Display dataset structure
str(isolates)
str(reference_genes)

# Display first six records
head(isolates)
head(reference_genes)

# Save imported data
saveRDS(isolates, "results/isolates_import_checkpoint.rds")

# Display completion message
cat("E. coli data import completed successfully!\n")

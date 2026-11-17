
# ==========================================
# AMR Genotype Analysis in E. coli
# Script 02: Data Cleaning
# ==========================================

# Import filtered E. coli dataset
isolates <- readRDS("data/ecoli_clean.rds")

# Check missing values
missing_summary <- data.frame(
  Column = names(isolates),
  Missing_Values = colSums(
    is.na(isolates) | isolates == ""
  )
)

print(missing_summary)

# Check duplicate records
duplicate_count <- sum(duplicated(isolates))

cat("Duplicate records:", duplicate_count, "\n")

# Check missing AMR genotype values
missing_genotypes <- sum(
  is.na(isolates$AMR.genotypes) |
    isolates$AMR.genotypes == ""
)

cat("Missing AMR genotype values:", missing_genotypes, "\n")

# Create cleaned working dataset
isolates_clean <- isolates

# Check remaining misplaced identifiers
remaining_errors <- sum(
  grepl("GCA_|SAMN", isolates_clean$AMR.genotypes),
  na.rm = TRUE
)

cat("Remaining identifier errors:", remaining_errors, "\n")

# Save cleaned dataset
saveRDS(
  isolates_clean,
  "results/isolates_clean.rds"
)

cat("E. coli data cleaning completed!\n")

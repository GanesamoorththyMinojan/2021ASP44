
# ==========================================
# AMR Genotype Analysis in E. coli
# Script 03: Exploratory Data Analysis
# ==========================================

# Load corrected E. coli dataset
isolates_clean <- readRDS("results/isolates_clean.rds")

# 1. Dataset dimensions
cat("Dataset Dimensions:\n")
print(dim(isolates_clean))

# 2. Dataset structure
cat("\nDataset Structure:\n")
str(isolates_clean)

# 3. First six records
cat("\nFirst Six Records:\n")
print(head(isolates_clean))

# 4. Missing values (NA and empty values)
missing_summary <- data.frame(
  Column = names(isolates_clean),
  Missing_Values = sapply(isolates_clean, function(x) {
    sum(is.na(x) | trimws(as.character(x)) == "")
  })
)

cat("\nMissing Values:\n")
print(missing_summary)

# Save missing value summary
write.csv(
  missing_summary,
  "results/missing_values_summary.csv",
  row.names = FALSE
)

# 5. Duplicate records
duplicate_count <- sum(duplicated(isolates_clean))

cat("\nDuplicate Records:\n")
print(duplicate_count)

# 6. Number of unique isolates
valid_isolates <- isolates_clean$Isolate[
  !is.na(isolates_clean$Isolate) &
    trimws(isolates_clean$Isolate) != ""
]

unique_isolates <- length(unique(valid_isolates))

cat("\nUnique Isolates:\n")
print(unique_isolates)

# 7. AMR genotype information
genotype_count <- sum(
  !is.na(isolates_clean$AMR.genotypes) &
    trimws(isolates_clean$AMR.genotypes) != ""
)

cat("\nRecords with AMR Genotype Information:\n")
print(genotype_count)

# 8. Save dataset summary
capture.output(
  summary(isolates_clean),
  file = "results/dataset_summary.txt"
)

cat("\nExploratory Data Analysis Completed!\n")

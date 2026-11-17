# ==========================================
# AMR Genotype Analysis in E. coli
# Script 04: AMR Gene Analysis
# ==========================================

# Create results folder if it does not exist
dir.create("results", showWarnings = FALSE)

# 1. Load cleaned E. coli dataset
isolates_clean <- readRDS("results/isolates_clean.rds")

# 2. Load NCBI reference gene catalogue
reference_genes <- read.delim(
  "data/ReferenceGeneCatalog.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 3. Extract valid AMR genotype records
amr_data <- isolates_clean[
  !is.na(isolates_clean$AMR.genotypes) &
    trimws(isolates_clean$AMR.genotypes) != "",
]

# 4. Split genotype information into individual genes
all_genes <- unlist(
  strsplit(
    amr_data$AMR.genotypes,
    split = ",",
    fixed = TRUE
  )
)

all_genes <- trimws(all_genes)

# Keep valid genotype entries
valid_tokens <- all_genes[
  grepl("=", all_genes, fixed = TRUE)
]

# Remove quotation marks
valid_tokens <- gsub('"', "", valid_tokens, fixed = TRUE)

# Extract gene names
gene_names <- sub("=.*", "", valid_tokens)
gene_names <- trimws(gene_names)

# Extract detection status
gene_status <- sub("^.*=", "", valid_tokens)
gene_status <- trimws(gene_status)

# 5. Calculate AMR gene frequencies
gene_counts <- sort(
  table(gene_names),
  decreasing = TRUE
)

gene_frequency_df <- data.frame(
  Gene = names(gene_counts),
  Frequency = as.numeric(gene_counts)
)

# Save gene frequency table
write.csv(
  gene_frequency_df,
  "results/amr_gene_frequencies.csv",
  row.names = FALSE
)

# 6. Display top 15 AMR genes
top_15_genes <- head(gene_frequency_df, 15)

cat("\nTop 15 AMR Genes:\n")
print(top_15_genes)

# 7. Plot top 15 AMR genes
png(
  "results/top_15_amr_genes.png",
  width = 1200,
  height = 800
)

barplot(
  rev(top_15_genes$Frequency),
  names.arg = rev(top_15_genes$Gene),
  horiz = TRUE,
  las = 1,
  main = "Top 15 Frequently Detected AMR Genes",
  xlab = "Number of Gene Detections",
  col = "steelblue"
)

dev.off()

# 8. Analyse AMR detection status
status_counts <- sort(
  table(gene_status),
  decreasing = TRUE
)

status_df <- data.frame(
  Status = names(status_counts),
  Frequency = as.numeric(status_counts)
)

cat("\nAMR Detection Status:\n")
print(status_df)

# Save detection status table
write.csv(
  status_df,
  "results/amr_detection_status.csv",
  row.names = FALSE
)

# Plot detection status
png(
  "results/amr_detection_status.png",
  width = 1200,
  height = 800
)

barplot(
  status_counts,
  las = 2,
  main = "AMR Gene Detection Status",
  ylab = "Number of Detections",
  col = "darkseagreen"
)

dev.off()

# 9. Match detected genes with NCBI reference catalogue
ncbi_reference <- unique(
  reference_genes[
    ,
    c("gene_family", "class", "subclass", "product_name")
  ]
)

matched_genes <- intersect(
  unique(gene_names),
  unique(ncbi_reference$gene_family)
)

gene_classification <- merge(
  data.frame(gene_family = matched_genes),
  ncbi_reference,
  by = "gene_family",
  all.x = TRUE
)

# Save gene classification
write.csv(
  gene_classification,
  "results/gene_classification.csv",
  row.names = FALSE
)

# 10. Calculate antimicrobial class frequencies
gene_frequency_df$gene_family <- gene_frequency_df$Gene

gene_class_frequency <- merge(
  gene_frequency_df,
  unique(
    gene_classification[
      ,
      c("gene_family", "class")
    ]
  ),
  by = "gene_family"
)

gene_class_frequency <- gene_class_frequency[
  !is.na(gene_class_frequency$class) &
    trimws(gene_class_frequency$class) != "",
]

class_counts <- sort(
  tapply(
    gene_class_frequency$Frequency,
    gene_class_frequency$class,
    sum
  ),
  decreasing = TRUE
)

class_df <- data.frame(
  Antimicrobial_Class = names(class_counts),
  Frequency = as.numeric(class_counts)
)

cat("\nAntimicrobial Class Frequencies:\n")
print(class_df)

# Save antimicrobial class frequencies
write.csv(
  class_df,
  "results/antimicrobial_class_frequencies.csv",
  row.names = FALSE
)

# 11. Plot top 10 antimicrobial classes
top_10_classes <- head(class_df, 10)

png(
  "results/top_10_antimicrobial_classes.png",
  width = 1200,
  height = 800
)

barplot(
  rev(top_10_classes$Frequency),
  names.arg = rev(top_10_classes$Antimicrobial_Class),
  horiz = TRUE,
  las = 1,
  main = "Top 10 Antimicrobial Resistance Classes",
  xlab = "Number of Gene Detections",
  col = "coral"
)

dev.off()

cat("\nAMR Gene Analysis Completed Successfully!\n")
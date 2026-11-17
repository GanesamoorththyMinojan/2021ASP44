# ==========================================
# AMR Genotype Analysis in E. coli
# Script 05: Hierarchical Clustering Analysis
# ==========================================

# 1. Load required package
library(vegan)

# 2. Load cleaned dataset
isolates_clean <- readRDS("results/isolates_clean.rds")

# 3. Identify the top 10 AMR genes
amr_data <- isolates_clean[
  !is.na(isolates_clean$AMR.genotypes) &
    trimws(isolates_clean$AMR.genotypes) != "",
]

all_genes <- unlist(
  strsplit(
    amr_data$AMR.genotypes,
    split = ",",
    fixed = TRUE
  )
)

all_genes <- trimws(all_genes)
all_genes <- gsub('"', "", all_genes, fixed = TRUE)

valid_tokens <- all_genes[
  grepl("=", all_genes, fixed = TRUE)
]

gene_names <- sub("=.*", "", valid_tokens)
gene_names <- trimws(gene_names)

gene_counts <- sort(
  table(gene_names),
  decreasing = TRUE
)

top_genes <- names(head(gene_counts, 10))

cat("Top 10 AMR Genes:\n")
print(top_genes)

# 4. Select 500 random isolates
valid_rows <- which(
  !is.na(isolates_clean$AMR.genotypes) &
    trimws(isolates_clean$AMR.genotypes) != ""
)

set.seed(123)

sample_rows <- sample(valid_rows, 500)

sample_genotypes <- isolates_clean$AMR.genotypes[
  sample_rows
]

# 5. Create gene presence/absence matrix
gene_matrix <- matrix(
  0,
  nrow = 500,
  ncol = length(top_genes)
)

colnames(gene_matrix) <- top_genes

for (i in seq_along(sample_genotypes)) {
  
  tokens <- unlist(
    strsplit(
      sample_genotypes[i],
      split = ",",
      fixed = TRUE
    )
  )
  
  tokens <- trimws(tokens)
  tokens <- gsub('"', "", tokens, fixed = TRUE)
  
  detected_genes <- tokens[
    grepl("=", tokens, fixed = TRUE)
  ]
  
  detected_genes <- sub("=.*", "", detected_genes)
  detected_genes <- trimws(detected_genes)
  
  gene_matrix[i, colnames(gene_matrix) %in% detected_genes] <- 1
}

cat("\nGene Presence Counts:\n")
print(colSums(gene_matrix))

# Save presence/absence matrix
write.csv(
  gene_matrix,
  "results/gene_presence_absence_matrix.csv"
)

# 6. Calculate Jaccard distance
gene_distance <- vegdist(
  gene_matrix,
  method = "jaccard",
  binary = TRUE
)

# 7. Perform hierarchical clustering
gene_cluster <- hclust(
  gene_distance,
  method = "average"
)

# 8. Plot dendrogram
png(
  "results/hierarchical_clustering.png",
  width = 1200,
  height = 800
)

plot(
  gene_cluster,
  labels = FALSE,
  main = "Hierarchical Clustering of E. coli Isolates",
  xlab = "500 Isolates",
  ylab = "Jaccard Distance",
  sub = ""
)

dev.off()

# 9. Divide isolates into three clusters
cluster_groups <- cutree(
  gene_cluster,
  k = 3
)

cluster_counts <- table(cluster_groups)

cat("\nCluster Sizes:\n")
print(cluster_counts)

# Save cluster assignments
cluster_results <- data.frame(
  Sample_Row = sample_rows,
  Cluster = cluster_groups
)

write.csv(
  cluster_results,
  "results/isolate_clusters.csv",
  row.names = FALSE
)

# 10. Calculate cluster gene profiles
cluster_data <- as.data.frame(gene_matrix)

cluster_data$Cluster <- factor(cluster_groups)

cluster_profile <- aggregate(
  . ~ Cluster,
  data = cluster_data,
  FUN = mean
)

cat("\nCluster Gene Profiles:\n")
print(cluster_profile)

write.csv(
  cluster_profile,
  "results/cluster_gene_profile.csv",
  row.names = FALSE
)

# 11. Create improved cluster heatmap
png(
  "results/amr_cluster_heatmap.png",
  width = 1800,
  height = 1100,
  res = 200
)

heatmap(
  t(gene_matrix),
  Rowv = NA,
  Colv = as.dendrogram(gene_cluster),
  scale = "none",
  col = c("white", "steelblue"),
  main = "AMR Gene Presence/Absence Heatmap",
  xlab = "500 Selected E. coli Isolates",
  ylab = "",
  labCol = rep("", nrow(gene_matrix)),
  margins = c(6, 10),
  cexRow = 1.1
)

dev.off()

cat("\nHierarchical Clustering Analysis Completed!\n")
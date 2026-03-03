#!/usr/bin/env Rscript

# Purpose: Cluster top 200 variable yeast stress-response genes and export heatmap + cluster assignments.

suppressPackageStartupMessages({
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required. Install it with: install.packages('pheatmap')")
  }
  library(pheatmap)
})

# Input and output paths
input_file <- "results/yeast_stress_cv_top200.tsv"
output_dir <- "results"
heatmap_file <- file.path(output_dir, "cv_top200_cluster_heatmap_02.png")
cluster_file <- file.path(output_dir, "cluster_assignment.tsv")

# Create results/ directory if needed
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load data
df <- read.delim(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# Validate required column
if (!("gene_id" %in% colnames(df))) {
  stop("Input file must contain a 'gene_id' column.")
}

# Identify expression columns: all numeric columns after gene_id, then take first 30 in original order
numeric_cols <- sapply(df, is.numeric)
# Ensure gene_id isn't accidentally numeric
numeric_cols[which(colnames(df) == "gene_id")] <- FALSE

expr_colnames <- colnames(df)[numeric_cols]

# If expression columns are missing (e.g., summary-only file), load from Gasch 2000 matrix
if (length(expr_colnames) < 30) {
  gasch_file <- "data/gasch2000.txt"
  if (!file.exists(gasch_file)) {
    stop("Fewer than 30 numeric expression columns found and data/gasch2000.txt is missing.")
  }

  gasch <- read.delim(gasch_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  if (!("UID" %in% colnames(gasch))) {
    stop("data/gasch2000.txt must contain a 'UID' column for gene IDs.")
  }

  # Stress-condition columns: all columns except metadata, preserving original order
  meta_cols <- c("UID", "NAME", "GWEIGHT")
  cond_cols <- setdiff(colnames(gasch), meta_cols)
  if (length(cond_cols) < 30) {
    stop("data/gasch2000.txt has fewer than 30 condition columns.")
  }
  cond_cols <- cond_cols[1:30]

  # Subset to the requested genes and coerce to numeric
  idx <- match(df$gene_id, gasch$UID)
  if (any(is.na(idx))) {
    missing <- sum(is.na(idx))
    stop(paste0("Missing ", missing, " gene_id values in data/gasch2000.txt."))
  }

  expr_mat <- as.matrix(gasch[idx, cond_cols, drop = FALSE])
  expr_mat <- apply(expr_mat, 2, function(x) as.numeric(x))
  rownames(expr_mat) <- df$gene_id
} else {
  # Use expression columns from input file (first 30 in original order)
  expr_colnames <- expr_colnames[1:30]
  expr_mat <- as.matrix(df[, expr_colnames, drop = FALSE])
  # Ensure numeric conversion in case columns were read as character
  expr_mat <- apply(expr_mat, 2, function(x) as.numeric(x))
  rownames(expr_mat) <- df$gene_id
}

# Row-wise Z-score normalization using exact formula (x - row_mean) / row_sd
row_mean <- rowMeans(expr_mat, na.rm = TRUE)
row_sd <- apply(expr_mat, 1, sd, na.rm = TRUE)

z_mat <- expr_mat
for (i in seq_len(nrow(expr_mat))) {
  if (is.na(row_sd[i]) || row_sd[i] == 0) {
    # If sd is zero, set entire row to 0 for consistent normalization
    z_mat[i, ] <- 0
  } else {
    z_mat[i, ] <- (expr_mat[i, ] - row_mean[i]) / row_sd[i]
  }
}

# Hierarchical clustering on rows only using Euclidean distance and ward.D2
row_dist <- dist(z_mat, method = "euclidean")
hc <- hclust(row_dist, method = "ward.D2")
clusters <- cutree(hc, k = 4)

# Prepare cluster assignment table
cluster_df <- data.frame(
  gene_id = rownames(z_mat),
  cluster = as.integer(clusters),
  stringsAsFactors = FALSE
)
cluster_df <- cluster_df[order(cluster_df$cluster, cluster_df$gene_id), ]

# Write cluster assignment file
write.table(cluster_df, file = cluster_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Annotation for heatmap
annotation_row <- data.frame(cluster = factor(clusters, levels = 1:4))
rownames(annotation_row) <- rownames(z_mat)

annotation_colors <- list(
  cluster = setNames(c("#1b9e77", "#d95f02", "#7570b3", "#e7298a"), levels(annotation_row$cluster))
)

# Generate heatmap PDF
png(file = heatmap_file, width = 8, height = 12, units = "in", res = 300, bg = "white")
pheatmap(
  z_mat,
  cluster_rows = hc,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100),
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  border_color = NA,
  fontsize = 10,
  fontsize_row = 6,
  fontsize_col = 9,
  main = "Top 200 Variable Yeast Stress Genes (Z-score, First 30 Conditions)",
  silent = TRUE
)
dev.off()

message("Done. Outputs written to: ", heatmap_file, " and ", cluster_file)

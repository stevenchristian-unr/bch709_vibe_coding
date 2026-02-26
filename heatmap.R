#!/usr/bin/env Rscript

# Gene heatmap from Gasch 2000 data

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(RColorBrewer)
})

url <- "https://www.shackett.org/files/gasch2000.txt"
local_path <- file.path("data", "gasch2000.txt")

if (!dir.exists("data")) dir.create("data", recursive = TRUE)
if (!file.exists(local_path)) {
  message("Downloading data...")
  download.file(url, destfile = local_path, mode = "wb", quiet = TRUE)
}

# Read data
# Read data with robust parsing (quoted header fields present)
raw <- read.delim(
  local_path,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  fill = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Identify gene and data columns
lower_names <- tolower(colnames(raw))
gene_col <- if (any(lower_names == "uid")) which(lower_names == "uid")[1] else 1

gweight_col <- if (any(lower_names == "gweight")) which(lower_names == "gweight")[1] else gene_col

genes <- raw[[gene_col]]
mat <- raw[, (gweight_col + 1):ncol(raw), drop = FALSE]
mat <- as.matrix(mat)

# Coerce to numeric
mat <- apply(mat, 2, function(x) as.numeric(as.character(x)))
rownames(mat) <- genes

# Select the first 30 genes and first 10 conditions
n_rows <- min(10, nrow(mat))
n_cols <- min(30, ncol(mat))
mat_top <- mat[seq_len(n_rows), seq_len(n_cols), drop = FALSE]

# Data are typically log2 ratios already; keep as-is and label accordingly
log_label <- "Gene Expression (log scale)"

# Prepare data for plotting
plot_df <- as.data.frame(mat_top)
plot_df$gene <- rownames(mat_top)
plot_long <- pivot_longer(plot_df, -gene, names_to = "condition", values_to = "value")

# Ensure factors for ordering
plot_long$gene <- factor(plot_long$gene, levels = rownames(mat_top))
plot_long$condition <- factor(plot_long$condition, levels = colnames(mat_top))

# Color palette (fallback to manual PuGn if unavailable)
if ("PuGn" %in% rownames(brewer.pal.info)) {
  pal <- brewer.pal(10, "PuGn")
} else {
  pal <- c("#f6eff7", "#d0d1e6", "#a6bddb", "#74a9cf", "#3690c0", "#0570b0", "#045a8d", "#023858", "#016c59")
}

p <- ggplot(plot_long, aes(x = gene, y = condition, fill = value)) +
  geom_tile(color = "black", linewidth = 0.2) +
  scale_fill_gradientn(colors = pal, name = log_label) +
  labs(title = "Gene Heatmap", x = "gene", y = "conditions") +
  theme_minimal(base_family = "Arial", base_size = 20) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = c(-0.75, 0.05),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  ) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

if (!dir.exists("results")) dir.create("results", recursive = TRUE)

out_file <- file.path("results", "gene_heatmap.png")

ggsave(out_file, plot = p, width = 10, height = 10, units = "in", dpi = 300, bg = "white")

message("Saved heatmap to: ", out_file)

# Short description
message("Description: Heatmap of the first 30 genes across the first 10 conditions, using a PuGn palette on a log scale.")

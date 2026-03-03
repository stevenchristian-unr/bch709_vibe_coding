data <- read.table("data/gasch2000.txt", header=TRUE, sep="\t", row.names=1)
vars <- apply(data, 1, var, na.rm=TRUE)
top200 <- head(sort(vars, decreasing=TRUE), 200)
heatmap(as.matrix(data[names(top200), ]))

library(data.table)
library(pheatmap)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

f <- "data/gasch2000.txt"
dt <- fread(f, header = TRUE)

# Extract gene_id and numeric condition columns
gene_ids <- dt[[1]]
meta_cols <- c("UID", "NAME", "GWEIGHT")
desc_col <- names(dt)[3]
skip_cols <- c(meta_cols, desc_col)
cond_cols <- setdiff(names(dt), skip_cols)
mat <- as.matrix(dt[, ..cond_cols])
mode(mat) <- "numeric"

mean_expr <- rowMeans(mat, na.rm = TRUE)
sd_expr   <- apply(mat, 1, sd, na.rm = TRUE)

keep <- !is.na(sd_expr) & sd_expr > 0
res <- data.table(
  gene_id   = gene_ids[keep],
  mean_expr = mean_expr[keep],
  sd_expr   = sd_expr[keep]
)
res[, cv := sd_expr / abs(mean_expr)]
res <- res[is.finite(cv)]
setorder(res, -cv)
top200 <- res[1:min(200, .N)]

# Summary TSV
top200_out <- copy(top200)
top200_out[, `:=`(
  mean_expr = round(mean_expr, 4),
  sd_expr   = round(sd_expr, 4),
  cv        = round(cv, 4)
)]
fwrite(top200_out, "results/yeast_stress_cv_top200.tsv", sep = "\t")

# Heatmap matrix
idx <- match(top200$gene_id, gene_ids)
submat <- mat[idx, , drop = FALSE]
rownames(submat) <- top200$gene_id

# Save (1800×1200 px, dpi 200)
png("results/yeast_stress_cv_top200_heatmap.png",
    width = 1800, height = 1200, res = 200)
pheatmap(
  submat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_col = 4,
  fontsize_row = 3,
  angle_col = 90,
  main = "Yeast stress response, CV top200 (Gasch et al. 2000)"
)
dev.off()

print(top200_out[1:min(10, .N)])
cat("Saved: results/yeast_stress_cv_top200.tsv\n")
cat("Saved: results/yeast_stress_cv_top200_heatmap.png\n")
library(data.table)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

f <- "data/gasch2000.txt"
dt <- fread(f, header = TRUE)

# Extract gene_id and numeric condition columns (skip NAME, description, GWEIGHT)
gene_ids <- dt[[1]]
meta_cols <- c("UID", "NAME", "GWEIGHT")
desc_col <- names(dt)[3]  # description column
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
top200[, `:=`(
  mean_expr = round(mean_expr, 4),
  sd_expr   = round(sd_expr, 4),
  cv        = round(cv, 4)
)]

fwrite(top200, "results/yeast_stress_cv_top200.tsv", sep = "\t")
print(top200[1:min(10, .N)])
cat("Saved: results/yeast_stress_cv_top200.tsv\n")
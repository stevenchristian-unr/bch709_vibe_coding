# scripts/variable_stress_genes.R
options(stringsAsFactors = FALSE)

path <- "data/gasch2000.txt"

# Load expression matrix (Gasch 2000 stress data)
x <- read.delim(path, sep = "\t", header = TRUE, quote = "", na.strings = c("", "NA", "NaN"))

# Numeric expression columns start at column 4
num <- x[, 4:ncol(x)]
num <- as.data.frame(lapply(num, function(col) as.numeric(col)))

# Per-gene variability (SD across conditions)
sdv <- apply(num, 1, sd, na.rm = TRUE)

# Extract standard gene name from NAME field (2nd token)
name_tokens <- strsplit(x$NAME, "\\s+")
name_tokens <- lapply(name_tokens, function(t) t[t != ""])
std_name <- vapply(name_tokens, function(t) if (length(t) >= 2) t[2] else NA_character_, character(1))

res <- data.frame(UID = x$UID, STD = std_name, SD = sdv, stringsAsFactors = FALSE)
res <- res[order(res$SD, decreasing = TRUE), ]

cat("Top 20 most variable genes:\\n")
print(head(res, 20), row.names = FALSE)

# Coordination test: pairwise correlation among top 50 variable genes
k <- 50
idx <- order(sdv, decreasing = TRUE)[1:k]
mat <- as.matrix(num[idx, ])
cor_mat <- cor(t(mat), use = "pairwise.complete.obs")
upper <- cor_mat[upper.tri(cor_mat)]

cat("\\nTop50 pairwise correlation summary (Pearson):\\n")
print(summary(upper))
cat("Top50 fraction corr > 0.5:", mean(upper > 0.5, na.rm = TRUE), "\\n")

# Random baseline for comparison
set.seed(1)
rand_idx <- sample(1:nrow(num), k)
mat2 <- as.matrix(num[rand_idx, ])
cor_mat2 <- cor(t(mat2), use = "pairwise.complete.obs")
upper2 <- cor_mat2[upper.tri(cor_mat2)]

cat("\\nRandom50 pairwise correlation summary (Pearson):\\n")
print(summary(upper2))
cat("Random50 fraction corr > 0.5:", mean(upper2 > 0.5, na.rm = TRUE), "\\n")

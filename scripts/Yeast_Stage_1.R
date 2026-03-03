data <- read.table("data/gasch2000.txt", header=TRUE, sep="\t", row.names=1)
vars <- apply(data, 1, var, na.rm=TRUE)
top200 <- head(sort(vars, decreasing=TRUE), 200)
heatmap(as.matrix(data[names(top200), ]))
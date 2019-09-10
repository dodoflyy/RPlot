library(tidyverse)
library(pheatmap)
library(BuenColors)

print("[Expr], [DEGs], [OutputName]")
args <- commandArgs(TRUE)
stopifnot(length(args) >= 3)
expr_path <- args[1]
deg_path <- args[2]
output <- args[3]

pdfOut <- paste(output, "pdf", sep = ".")
pngOut <- paste(output, "png", sep = ".")

expr <- read_tsv(expr_path)
deg <- read_tsv(deg_path)
# 要求至少2样品
stopifnot(length(colnames(expr)) >= 3)
stopifnot("ENSEMBL" %in% colnames(expr), "ENSEMBL" %in% colnames(deg))

expr2 <- filter(expr, ENSEMBL %in% deg$ENSEMBL) %>% as.data.frame()
dim(expr2) %>% print()
rownames(expr2) <- expr2$ENSEMBL
expr2$ENSEMBL <- NULL

expr3 <- t(expr2)
expr4 <- scale(expr3)
expr2 <- t(expr4)
head(expr2, n = 3)

exprV <- as.vector(expr2)
maxV <- max(exprV)
minV <- min(exprV)

# color_color <- c(colorRampPalette(colors = c("#28B463", "#D0D3D4"))(100), colorRampPalette(colors=c("#D0D3D4", "#A93226"))(100))
color_color <- jdb_palette("Zissou", type = "continuous")
pheatmap(expr2, cluster_rows = TRUE, cluster_cols = TRUE, legend = TRUE, show_colnames = TRUE, show_rownames = FALSE, color = color_color, scale = "none", main = "DEGs", legend_breaks = c(minV, maxV), legend_labels = c("Low", "High"), filename = pdfOut)
pheatmap(expr2, cluster_rows = TRUE, cluster_cols = TRUE, legend = TRUE, show_colnames = TRUE, show_rownames = FALSE, color = color_color, scale = "none", main = "DEGs", legend_breaks = c(minV, maxV), legend_labels = c("Low", "High"), filename = pngOut)
print("完成")

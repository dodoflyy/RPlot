library(tidyverse)
library(pheatmap)
library(BuenColors)

print("[表达文件], [差异基因文件], [输出目录] [输出文件名] [图片标题(可选)]")
args <- commandArgs(TRUE)
stopifnot(length(args) >= 4)
expr_path <- args[1]
deg_path <- args[2]
outDir <- args[3]
name <- args[4]

if(length(args) >= 5){
  titleText <- args[5]
  }else{
  titleText <- "Gene expression heatmap"
  }

pdfOut <- str_glue("{outDir}/{name}.pdf")
pngOut <- str_glue("{outDir}/{name}.png")

expr <- read_tsv(expr_path)
deg <- read_tsv(deg_path)
dim(expr) %>% print()
# 要求至少2样品
stopifnot(length(colnames(expr)) >= 3)

if(all("ENSEMBL" %in% colnames(expr), "ENSEMBL" %in% colnames(deg))){
  expr2 <- filter(expr, ENSEMBL %in% deg$ENSEMBL) %>% as.data.frame()
  }else if(all("entrezgene_id" %in% colnames(expr), "entrezgene_id" %in% colnames(deg))){
  expr2 <- filter(expr, entrezgene_id %in% deg$entrezgene_id) %>% as.data.frame()
  }else{
  stop("输入2表格要求同时有 ENSEMBL 列或者 entrezgene_id 列")
  }

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
pheatmap(expr2, cluster_rows = TRUE, cluster_cols = TRUE, legend = TRUE, show_colnames = TRUE, show_rownames = FALSE, color = color_color, scale = "none", main = titleText, legend_breaks = c(minV, maxV), legend_labels = c("Low", "High"), filename = pdfOut)
pheatmap(expr2, cluster_rows = TRUE, cluster_cols = TRUE, legend = TRUE, show_colnames = TRUE, show_rownames = FALSE, color = color_color, scale = "none", main = titleText, legend_breaks = c(minV, maxV), legend_labels = c("Low", "High"), filename = pngOut)
print("完成")

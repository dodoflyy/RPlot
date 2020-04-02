# 用 DESeq2 产生的结果文件生成差异基因表达热图
# 默认结果文件包含 ensembl_gene_id, hgnc_symbol, entrezgene_id 3列基因ID，可以选择任意一个去做热图
# 脚本在 R 3.6 环境测试通过
# 需要以下 R 包支持
# tidyverse, pheatmap, BuenColors

writeLines("Rscript pHeatmap.R Expression.csv DEG.csv OutputDir Filename GeneIdType [\"Image Title\"]\n")
args <- commandArgs(TRUE)
stopifnot(length(args) >= 5)

library(tidyverse)
library(pheatmap)
library(BuenColors)

expr_path <- args[1]
deg_path <- args[2]
outDir <- args[3]
name <- args[4]
geneIdType <- args[5]

if(length(args) >= 6){
  titleText <- args[6]
  }else{
  titleText <- "Gene expression heatmap"
  }

pdfOut <- str_glue("{outDir}/{name}.pdf")
pngOut <- str_glue("{outDir}/{name}.png")

expr <- read_csv(expr_path)
deg <- read_csv(deg_path)
dim(expr) %>% print()
# 要求至少2样品
stopifnot(ncol(expr) >= 5)

if (geneIdType == "ensembl_gene_id") {
  expr2 <- dplyr::select(expr, -hgnc_symbol, -entrezgene_id) %>% dplyr::filter(ensembl_gene_id %in% deg$ensembl_gene_id) %>% dplyr::distinct(ensembl_gene_id, .keep_all=TRUE) %>% dplyr::rename(gene_id=ensembl_gene_id) %>% as.data.frame()
} else if (geneIdType == "hgnc_symbol") {
  expr2 <- dplyr::select(expr, -ensembl_gene_id, -entrezgene_id) %>% dplyr::filter(!is.na(hgnc_symbol) & hgnc_symbol != "") %>% dplyr::filter(hgnc_symbol %in% deg$hgnc_symbol) %>% dplyr::distinct(hgnc_symbol, .keep_all=TRUE) %>% dplyr::rename(gene_id=hgnc_symbol) %>% as.data.frame()
} else {
  expr2 <- dplyr::select(expr, -ensembl_gene_id, -hgnc_symbol) %>% dplyr::filter(!is.na(entrezgene_id)) %>% dplyr::filter(entrezgene_id %in% deg$entrezgene_id) %>% dplyr::distinct(entrezgene_id, .keep_all=TRUE) %>% dplyr::rename(gene_id=entrezgene_id) %>% as.data.frame()
}

head(expr2, n=3)
dim(expr2)
rownames(expr2) <- expr2$gene_id
expr2$gene_id <- NULL

# 针对最终的行(基因)进行 scale
expr3 <- t(expr2)
expr4 <- scale(expr3)
expr2 <- t(expr4)

exprV <- as.vector(expr2)
maxV <- max(exprV)
minV <- min(exprV)

# color_color <- c(colorRampPalette(colors = c("#28B463", "#D0D3D4"))(100), colorRampPalette(colors=c("#D0D3D4", "#A93226"))(100))
color_color <- jdb_palette("Zissou", type = "continuous")
pheatmap(expr2, cluster_rows = TRUE, cluster_cols = TRUE, legend = TRUE, show_colnames = TRUE, show_rownames = FALSE, color = color_color, scale = "none", main = titleText, legend_breaks = c(minV, maxV), legend_labels = c("Low", "High"), filename = pdfOut)
pheatmap(expr2, cluster_rows = TRUE, cluster_cols = TRUE, legend = TRUE, show_colnames = TRUE, show_rownames = FALSE, color = color_color, scale = "none", main = titleText, legend_breaks = c(minV, maxV), legend_labels = c("Low", "High"), filename = pngOut)
writeLines("完成")

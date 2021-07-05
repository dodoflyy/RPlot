# 用 DESeq2 产生的结果文件生成差异基因表达热图
# 默认结果文件包含 ensembl_gene_id, hgnc_symbol, entrezgene_id 3列基因ID，可以选择任意一个去做热图
# 脚本在 R 3.6 环境测试通过
# 需要以下 R 包支持
# argparse, tidyverse, pheatmap, BuenColors

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(BuenColors))

scriptDescription="热图展示 DESeq2 差异基因分析结果"
parser <- ArgumentParser(description=scriptDescription, add_help=TRUE)
parser$add_argument("--Expression", dest="EXPR", help="csv 格式基因表达数据", required=TRUE)
parser$add_argument("--DEGs", dest="DEGs", help="DESeq2 差异基因分析结果", required=TRUE)
parser$add_argument("--OutputDir", dest="OUT", help="输出目录", required=TRUE)
parser$add_argument("--basename", dest="BASE", help="输出文件名", required=TRUE)
parser$add_argument("--Title", dest="TITLE", help="热图标题", default="Gene expression heatmap")
parser$add_argument("--GeneType", dest="GT", help="基因名类型，可选值为：hgnc_symbol, ensembl_gene_id, entrezgene_id", 
                    choices=c("hgnc_symbol", "ensembl_gene_id", "entrezgene_id"), required=TRUE)

argvs <- parser$parse_args()
exprPath <- file.path(argvs$EXPR)
degPath <- file.path(argvs$DEGs)
outDir <- file.path(argvs$OUT)
baseName <- argvs$BASE
titleText <- argvs$TITLE
geneIdType <- argvs$GT

pdfName <- paste(baseName, "Heatmap", "pdf", sep=".")
pngName <- paste(baseName, "Heatmap", "png", sep=".")
pdfOut <- file.path(outDir, pdfName)
pngOut <- file.path(outDir, pngName)

expr <- read_csv(exprPath)
deg <- read_csv(degPath)
# 要求至少 4 样品
if (ncol(expr) < 4) {
  writeLines("\n样品数目小于 4 请检查数据和参数！")
  q(save="no")
}

if (geneIdType == "ensembl_gene_id") {
  expr2 <- dplyr::select(expr, -hgnc_symbol, -entrezgene_id) %>% dplyr::filter(ensembl_gene_id %in% deg$ensembl_gene_id) %>% 
    dplyr::distinct(ensembl_gene_id, .keep_all=TRUE) %>% dplyr::rename(gene_id=ensembl_gene_id) %>% as.data.frame()
} else if (geneIdType == "hgnc_symbol") {
  expr2 <- dplyr::select(expr, -ensembl_gene_id, -entrezgene_id) %>% dplyr::filter(!is.na(hgnc_symbol) & hgnc_symbol != "") %>% 
    dplyr::filter(hgnc_symbol %in% deg$hgnc_symbol) %>% dplyr::distinct(hgnc_symbol, .keep_all=TRUE) %>% dplyr::rename(gene_id=hgnc_symbol) %>% as.data.frame()
} else {
  expr2 <- dplyr::select(expr, -ensembl_gene_id, -hgnc_symbol) %>% dplyr::filter(!is.na(entrezgene_id)) %>% 
    dplyr::filter(entrezgene_id %in% deg$entrezgene_id) %>% dplyr::distinct(entrezgene_id, .keep_all=TRUE) %>% 
    dplyr::rename(gene_id=entrezgene_id) %>% as.data.frame()
}

rownames(expr2) <- expr2$gene_id
expr2$gene_id <- NULL

# 针对最终的行(基因)进行 scale
expr3 <- t(expr2)
expr4 <- scale(expr3)
expr2 <- t(expr4)

exprV <- as.vector(expr2)
maxV <- max(exprV)
minV <- min(exprV)

color_color <- jdb_palette("Zissou", type = "continuous")
pheatmap(expr2, cluster_rows = TRUE, cluster_cols = TRUE, legend = TRUE, show_colnames = TRUE, 
         show_rownames = FALSE, color = color_color, scale = "none", main = titleText, 
         legend_breaks = c(minV, maxV), legend_labels = c("Low", "High"), filename = pdfOut)
pheatmap(expr2, cluster_rows = TRUE, cluster_cols = TRUE, legend = TRUE, show_colnames = TRUE, 
         show_rownames = FALSE, color = color_color, scale = "none", main = titleText, 
         legend_breaks = c(minV, maxV), legend_labels = c("Low", "High"), filename = pngOut)

writeLines("\n完成！")

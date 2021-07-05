# 从表达矩阵绘制热图，默认行为基因，列为样品
# 第一列为基因名，命名为 gene_id
# 输出 pdf 
# 在 R 3.6 环境测试通过
# 需要以下 R 包支持
# argparse, tidyverse, ComplexHeatmap, BuenColors

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(BuenColors))

scriptDescription <- "绘制表达热图"
parser <- ArgumentParser(description=scriptDescription, add_help=TRUE)
parser$add_argument("--Expression", dest="EXPR", help="csv 格式表达数据，第一列为基因名 gene_id", required=TRUE)
parser$add_argument("--OutputDir", dest="OUT", help="输出目录", default=".")
parser$add_argument("--Basename", dest="BASE", help="输出文件名，默认 \"Unknow\"", default="Unknow")
parser$add_argument("--Title", dest="TITLE", help="图像标题，默认 \"Gene Expression\"", default="Gene Expression")
parser$add_argument("--GeneID", dest="GENEID", help="基因名的列名，默认 \"gene_id\"", default="gene_id")
parser$add_argument("--Zscore", dest="ZSCORE", action="store_true", help="是否对数据进行 Z score 处理？")

argvs <- parser$parse_args()
expressionPath <- file.path(argvs$EXPR)
outputDir <- file.path(argvs$OUT)
titleText <- argvs$TITLE
baseName <- argvs$BASE
geneColumnName <- argvs$GENEID
zScore <- argvs$ZSCORE
fileName <- paste(baseName, "Heatmap", "pdf", sep=".")
heatmapPath <- file.path(outputDir, fileName)

exprData <- read_csv(expressionPath) %>% dplyr::distinct(gene_id, .keep_all=TRUE) %>% as.data.frame()
rownames(exprData) <- exprData$gene_id
exprData$gene_id <- NULL

# 是否 Z Score 转换
if (zScore) {
  plotData <- as.matrix(exprData) %>% t() %>% scale() %>% t()
} else {
  plotData <- as.matrix(exprData)
}

colorFun <- jdb_palette("Zissou", type="continuous")
hm <- Heatmap(plotData, name="Normalized Expression", col=colorFun, cluster_rows = TRUE, cluster_columns = FALSE, 
              show_row_names = FALSE, column_names_side = "top", column_title = titleText)

pdf(heatmapPath)
draw(hm)
dev.off()

writeLines("\n完成！")


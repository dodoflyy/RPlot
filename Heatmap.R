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
parser$add_argument("--Expression", "-E", dest="EXPR", help="csv 格式表达数据，第一列为基因名 gene_id", required=TRUE)
parser$add_argument("--OutputDir", "-O", dest="OUT", help="输出目录", required=TRUE)
parser$add_argument("--Basename", "-B", dest="BASE", help="输出文件名", required=TRUE)
parser$add_argument("--Title", "-T", dest="TITLE", help="图像标题", default="Gene Expression")

argvs <- parser$parse_args()
expressionPath <- file.path(argvs$EXPR)
outputDir <- file.path(argvs$OUT)
fileName <- paste(argvs$BASE, "Heatmap", "pdf", sep=".")
heatmapPath <- file.path(outputDir, fileName)
titleText <- argvs$TITLE

writeLines("\n====== 读取表达数据 ======")
Data <- read_csv(expressionPath) %>% dplyr::distinct(gene_id, .keep_all=TRUE) %>% as.data.frame()
rownames(Data) <- Data$gene_id
Data$gene_id <- NULL
Data2 <- t(Data) %>% scale() %>% t()
head(Data2, n=3)

color <- jdb_palette("Zissou", type="continuous")
hm <- Heatmap(Data2, name="Expression(scaled)", col=color, cluster_rows = TRUE, cluster_columns = FALSE, 
              show_row_names = FALSE, column_names_side = "top", column_title = titleText)

pdf(heatmapPath)
draw(hm)
dev.off()

writeLines("\n完成！")


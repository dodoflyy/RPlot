# 从表达矩阵生成样本聚类树
# 默认输入矩阵行为基因，列为样本，第一列为 gene_id
# CutHeight 参数如果设置为 0 表示不进行 Cutree 只出一个聚类树
# 脚本在 R 3.6 测试通过
# 需要下列依赖包
# argparse, tidyverse

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))

scriptDescription="表达矩阵生成样本聚类树"
parser <- ArgumentParser(description=scriptDescription, add_help=TRUE)
parser$add_argument("--Expression", "-E", dest="EXPR", help="csv 格式表达数据", required=TRUE)
parser$add_argument("--CutHeight", "-C", dest="CUT", help="Cut Height", required=TRUE)
parser$add_argument("--OutputDir", "-O", dest="OUT", help="输出目录", required=TRUE)
parser$add_argument("--basename", "-B", dest="BASE", help="输出文件名", required=TRUE)
parser$add_argument("--Title", "-T", dest="TITLE", help="图像标题", default="Sample Cluster Dendrogram")

argvs <- parser$parse_args()
expressionPath <- file.path(argvs$EXPR)
cutHeight <- as.double(argvs$CUT)
outDir <- file.path(argvs$OUT)
baseName <- argvs$BASE
plotTitle <- argvs$TITLE

Data1 <- read_csv(expressionPath)
glimpse(Data1)
Data2 <- as.data.frame(Data1)
rownames(Data2) <- Data2$gene_id
Data2$gene_id <- NULL
Data3 <- t(Data2)

dd <- dist(Data3, method = "euclidean")
hc <- hclust(dd, method = "average")
pdfName <- paste(baseName, "CutTree", "pdf", sep=".")
pdfOut <- file.path(outDir, pdfName)
w <- (nrow(Data3) * 0.2) %>% as.integer()

pdf(pdfOut, width = w, height = 9)
plot(hc, xlab = "Sample", main = plotTitle)
if (cutHeight != 0) {
  abline(h = cutHeight, col = "red")
}
dev.off()

exprFileName <- paste(baseName, "KeepTree", "csv", sep=".")
exprOut <- file.path(outDir, exprFileName)
rmFileName <- paste(baseName, "RemoveSample", "txt", sep=".")
rmOut <- file.path(outDir, rmFileName)

if (cutHeight != 0) {
  # 只考虑分2组，那么最多的那一组就是需要保留的组，其余移除
  cutTree <- cutree(hc, h = cutHeight)
  cutTreeStat <- table(cutTree) %>% as.data.frame() %>% as_tibble() %>% dplyr::arrange(desc(Freq))
  keepGroup <- cutTreeStat[[1, 1]]
  keepSample <- cutTree[cutTree==keepGroup] %>% names()
  keepData <- dplyr::select(Data1, gene_id, keepSample)
  removeList <- setdiff(rownames(Data3), keepSample)
  write_csv(keepData, exprOut)
  write.table(removeList, file = rmOut, sep = "\n", row.names = FALSE, col.names = FALSE)
  print(removeList)
}


writeLines("\nヽ(✿ﾟ▽ﾟ)ノ")

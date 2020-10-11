# 从表达矩阵生成样本聚类树
# 默认输入矩阵行为基因，列为样本，第一列为 gene_id
# CutHeight 参数如果设置为 0 表示不进行 Cutree 只出一个聚类树
# 脚本在 R 3.6 测试通过
# 需要 tidyverse 包

writeLines("\n表达矩阵生成样本聚类树\n")
writeLines("Usage:")
writeLines("\nRscript SampleCutTree.R Expression.csv CutHeight OutputDir Filename [\"Plot Title\"]")
argvs <- commandArgs(trailingOnly = TRUE)
stopifnot(length(argvs) >= 4)

expressionPath <- file.path(argvs[1])
cutHeight <- argvs[2] %>% as.double()
outDir <- argvs[3]
fileName <- argvs[4]

writeLines(stringr::str_glue("表达数据：{expressionPath}"))
writeLines(stringr::str_glue("Cut Height: {cutHeight}"))
writeLines(stringr::str_glue("输出目录：{outDir}"))
writeLines(stringr::str_glue("输出文件名：{fileName}"))

if (length(argvs) >= 5) {
  plotTitle <- argvs[5]
  writeLines(stringr::str_glue("图像标题：{plotTitle}"))
} else {
  plotTitle <- "Sample Cluster Dendrogram"
}

library(tidyverse, quietly = TRUE, warn.conflicts = FALSE)

Data1 <- read_csv(expressionPath)
glimpse(Data1)
Data2 <- as.data.frame(Data1)
rownames(Data2) <- Data2$gene_id
Data2$gene_id <- NULL
Data3 <- t(Data2)

dd <- dist(Data3, method = "euclidean")
hc <- hclust(dd, method = "average")
pdfOut <- stringr::str_glue("{outDir}/{fileName}_CutTree.pdf")
w <- (nrow(Data3) * 0.2) %>% as.integer()
pdf(pdfOut, width = w, height = 9)
plot(hc, xlab = "Sample", main = plotTitle)
if (cutHeight != 0) {
  abline(h = cutHeight, col = "red")
}
dev.off()

if (cutHeight != 0) {
  # 只考虑分2组，那么最多的那一组就是需要保留的组，其余移除
  cutTree <- cutree(hc, h = cutHeight)
  cutTreeStat <- table(cutTree) %>% as.data.frame() %>% as_tibble() %>% dplyr::arrange(desc(Freq))
  keepGroup <- cutTreeStat[[1, 1]]
  keepSample <- cutTree[cutTree==keepGroup] %>% names()
  keepData <- dplyr::select(Data1, gene_id, keepSample)
  exprOut <- stringr::str_glue("{outDir}/{fileName}_TreeKeep.csv")
  write_csv(keepData, exprOut)
  removeList <- setdiff(rownames(Data3), keepSample)
  rmOut <- exprOut <- stringr::str_glue("{outDir}/{fileName}_RemoveSample.txt")
  write.table(removeList, file = rmOut, sep = "\n", row.names = FALSE, col.names = FALSE)
  writeLines("移除以下样本数据")
  print(removeList)
}


writeLines("\n完成")

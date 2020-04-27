## 脚本效果不好。因为 RNA-seq 数据，尤其是样本多且差异大时难以用 PCA 2个维度进行展示。
# 从 RNA-seq 数据用 Isolation Forests 检测离群样本
# 默认输入矩阵行为基因列为样本，第一列是 gene_id 
# 脚本在 R 3.6 测试通过
# 需要 tidyverse, solitude, ggrepel 包

writeLines("\nRscript IsolationForestsOutlier.R Expression.csv AnomalyScoreCutoff OutputDir Filename [\"Plot Title\"]\n")
argvs <- commandArgs(trailingOnly = TRUE)
stopifnot(length(argvs) >= 4)

library(tidyverse, quietly = TRUE)
library(solitude, quietly = TRUE)
library(ggrepel, quietly = TRUE)

if (length(argvs) >= 5) {
  plotTitle <- argvs[5]
} else {
  plotTitle <- "Outlier Sample"
}
outDir <- argvs[3]
fileName <- argvs[4]
asCutoff <- argvs[2] %>% as.double() %>% round(3)

Data <- read_csv(argvs[1])
Data2 <- as.data.frame(Data)
rownames(Data2) <- Data2$gene_id
Data2$gene_id <- NULL
Data3 <- t(Data2) %>% as.data.frame()

outlierStatus <- function(x) {
  if (x >= asCutoff) {
    return("Yes")
  } else {
    return("No")
  }
}

dataPca <- prcomp(Data3)
pcaSummary <- summary(dataPca)
# 取得 PC1,PC2 解释占比
summaryTable <- pcaSummary$importance
summaryTable [, 1:6]
pc1Por <- summaryTable[[2, 1]] %>% round(3)
pc2Por <- summaryTable[[2, 2]] %>% round(3)
pcaData <- dataPca$x %>% as_tibble(rownames="Sample") %>% dplyr::select(Sample, PC1, PC2)

# 因为基因一般非常多，把 sampleSize 设置大些
iforest <- isolationForest$new(mtry = 5, nproc=12, replace = TRUE, num_trees = 200)
iforest$fit(Data3)
ifsResult <- as_tibble(iforest$scores) %>% dplyr::mutate(Sample=rownames(Data3))
glimpse(ifsResult)
plotData <- dplyr::left_join(pcaData, ifsResult, by="Sample") %>% dplyr::mutate(Outlier=map_chr(anomaly_score, outlierStatus))
glimpse(plotData)

repelData <- data.frame(x=(max(pcaData$PC1) * 0.75), y=(max(pcaData$PC2) * 0.75), z=stringr::str_glue("Anomaly Score Cutoff: {asCutoff}"))
colFun <- c("Yes"="#F22300", "No"="#3C9AB2")
xLab <- stringr::str_glue("PC1({pc1Por})")
yLab <- stringr::str_glue("PC2({pc2Por})")
pcaPlot <- ggplot2::ggplot(plotData, aes(x=PC1, y=PC2)) +
           ggplot2::geom_point(aes(color=Outlier)) +
           ggplot2::scale_color_manual(values = colFun) +
           ggplot2::labs(title = plotTitle, x = xLab, y = yLab) +
           ggrepel::geom_text_repel(data = repelData, mapping = aes(x, y, label=z)) +
           ggplot2::theme_bw()

pdfOut <- stringr::str_glue("{fileName}_IsolationForest.pdf")
ggplot2::ggsave(filename = pdfOut, plot = pcaPlot, device = "pdf", path = outDir)

normalSample <- dplyr::filter(plotData, Outlier=="No")$Sample
outlierSample <- tibble(Sample=(dplyr::filter(plotData, Outlier=="Yes")$Sample))
normalData <- dplyr::select(Data, gene_id, normalSample)
dim(normalData)
outlierOut <- stringr::str_glue("{outDir}/{fileName}_OutlierList.txt")
normalOut <- stringr::str_glue("{outDir}/{fileName}_NormalSample.csv")
write_csv(outlierSample, outlierOut)
write_csv(normalData, normalOut)

writeLines("完成")


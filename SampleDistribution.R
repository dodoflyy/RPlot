# 提供 RNA-seq 表达矩阵，行为基因，列为样本。出图展示样品关系。
# 应输入 Log 后矩阵，比如 DESeq2 的 rlog 结果
# 假定输入第一列为基因ID，且列名为 gene_id
# 假定SampleGroup 有 Sample 和 Group 2列
# 在 R3.6 测试通过。
# 需要 tidyverse, BuenColors, ComplexHeatmap, ggrepel 包

writeLines("\n出图展示 RNA-seq 样品分布\n")
writeLines("Usage:")
writeLines("\nRscript SampleDistribution.R Expression.csv SampleGroup.csv Output.pdf\n")
argvs <- commandArgs(trailingOnly = TRUE)
stopifnot(length(argvs) >= 3)

expressionPath <- file.path(argvs[1])
groupPath <- file.path(argvs[2])
outPath <- file.path(argvs[3])
writeLines(stringr::str_glue("输入表达文件：{expressionPath}"))
writeLines(stringr::str_glue("样品分组文件：{groupPath}"))
writeLines(stringr::str_glue("输出图片文件：{outPath}"))

library(tidyverse, quietly = TRUE, warn.conflicts = FALSE)
library(BuenColors, quietly = TRUE, warn.conflicts = FALSE)
library(ComplexHeatmap, quietly = TRUE, warn.conflicts = FALSE)
library(ggrepel, quietly = TRUE, warn.conflicts = FALSE)

sampleGroup <- read_csv(groupPath)
Data0 <- read_csv(expressionPath)
if (all(c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol") %in% colnames(Data0))) {
  Data0 <- dplyr::select(Data0, -entrezgene_id, -hgnc_symbol) %>% dplyr::rename(gene_id=ensembl_gene_id)
}
Data <-  as.data.frame(Data0)
rownames(Data) <- Data$gene_id
Data$gene_id <- NULL
head(Data, n = 3)

dataPca <- prcomp(t(Data))
pcaSummary <- summary(dataPca)
print(pcaSummary)
# 取得 PC1,PC2 解释占比
summaryTable <- pcaSummary$importance
pc1Por <- summaryTable[[2, 1]] %>% round(3)
pc2Por <- summaryTable[[2, 2]] %>% round(3)
pcaData <- dataPca$x %>% as_tibble(rownames="Sample") %>% dplyr::select(Sample, PC1, PC2) %>% dplyr::left_join(sampleGroup, by="Sample")
glimpse(pcaData)
xLab <- stringr::str_glue("PC1({pc1Por})")
yLab <- stringr::str_glue("PC2({pc2Por})")
pcaPlot <- ggplot2::ggplot(pcaData, aes(x=PC1, y=PC2)) +
           geom_point(aes(shape=Group)) +
           ggrepel::geom_text_repel(aes(label=Sample)) +
           ggplot2::labs(title = "Sample PCA", x = xLab, y = yLab) +
           ggplot2::theme_bw()

dataDist <- dist(t(Data), method = "euclidean")
dataHc <- hclust(dataDist, method = "average")

col_fun <- jdb_palette("brewer_fire", type = "continuous") %>% rev()
hmData <- as.matrix(dataDist)
hm <- Heatmap(hmData, name = "Euclidean", col = col_fun, cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = TRUE, row_names_side = "right", show_row_dend = FALSE, show_column_dend = FALSE, column_title = "Sample Distance", column_title_side = "top")

pdf(outPath)
print(pcaPlot)
plot(dataHc, xlab = "Sample")
draw(hm)
dev.off()

writeLines("\n完成")
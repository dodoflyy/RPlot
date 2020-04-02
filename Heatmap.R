# 从表达矩阵绘制热图，默认行为基因，列为样品
# 第一列为基因名，命名为 gene_id
# 输出 pdf 
# 在 R 3.6 环境测试通过
# 需要以下 R 包支持
# tidyverse, ComplexHeatmap, BuenColors

writeLines("Rscript Heatmap.R ExpressionData.csv Heatmap.pdf [Title]\n")
argvs <- commandArgs(trailingOnly = TRUE)
stopifnot(length(argvs) >= 2)

library(tidyverse, quietly = TRUE)
library(ComplexHeatmap, quietly = TRUE)
library(BuenColors, quietly = TRUE)

Data <- read_csv(argvs[1]) %>% dplyr::distinct(gene_id, .keep_all=TRUE) %>% as.data.frame()
rownames(Data) <- Data$gene_id
Data$gene_id <- NULL
Data2 <- t(Data) %>% scale() %>% t()
head(Data2, n=3)

color <- jdb_palette("brewer_celsius", type="continuous")
if (length(argvs) >= 3) {
  titleText <- argvs[3]
  hm <- Heatmap(Data2, name="Expression(scaled)", col=color, cluster_rows = TRUE, cluster_columns = FALSE, show_row_names = FALSE, column_names_side = "top", column_title = titleText)
} else {
  hm <- Heatmap(Data2, name="Expression(scaled)", col=color, cluster_rows = TRUE, cluster_columns = FALSE, show_row_names = FALSE, column_names_side = "top")
}

pdf(argvs[2])
draw(hm)
dev.off()

writeLines("完成")


# 展示 DESeq2 分析的差异基因总体情况
# 包含将差异基因（padj < 0.05）排序后点图
# 以及每个基因差异倍数与总体差异倍数中位值的差值
# 可以看总体差异基因情况，可以用来挑选 Outlier 阈值（或许）
# 脚本在 R 3.6 测试通过
# 需要下列包依赖
# argparse, tidyverse

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))

scriptDescription <- "展示 DESeq2 分析的差异基因总体情况，可用来挑选合适的 Outlier 阈值"
parser <- ArgumentParser(description=scriptDescription, add_help=TRUE)
parser$add_argument("--DEGs", dest="DEGs", help="差异基因数据", required=TRUE)
parser$add_argument("--PdfPath", dest="PDF", help="输出 PDF 路径", required=TRUE)
parser$add_argument("--Pvalue", dest="PVAL", help="P 值阈值，默认 0.05", default=0.05)

argvs <- parser$parse_args()
degPath <- file.path(argvs$DEGs)
plotPath <- file.path(argvs$PDF)
pVal <- as.double(argvs$PVAL)

degData <- readr::read_csv(degPath) %>% dplyr::filter(!is.na(log2FoldChange) & padj < pVal)
if (nrow(degData) < 25) {
  writeLines("╮（╯＿╰）╭")
  writeLines("差异基因数目过少，请检查数据和参数！")
  q(save="no")
}

# Y 轴坐标固定 0.5 间隔
nmads_plot <- function(plot_data, plot_title) {
  yMin <- as.integer(min(plot_data$Nmads)) - 0.5
  yMax <- as.integer(max(plot_data$Nmads)) + 0.5
  plotObj <- ggplot2::ggplot(data = plot_data, mapping = aes(x = Gene, y = Nmads)) +
    ggplot2::geom_point(color = "#000000") +
    ggplot2::scale_y_continuous(breaks = seq(from = yMin, to = yMax, by = 0.5)) +
    ggplot2::labs(title = plot_title, x = "Gene", y = "Nmads") +
    ggplot2::theme_bw()
  return(plotObj)
}

# 不用与中位值的差的绝对值
nmads_data <- function(degData) {
  degMad <- median(degData$log2FoldChange)
  nmadsData <- dplyr::mutate(degData, Nmads=log2FoldChange - degMad) %>% 
    dplyr::arrange(Nmads) %>% dplyr::mutate(Gene = 1:nrow(degData))
  return(nmadsData)
}

plotData1 <- dplyr::arrange(degData, abs(log2FoldChange)) %>% dplyr::mutate(Gene = 1:nrow(degData))
plot1 <- ggplot2::ggplot(data = plotData1, mapping = aes(x = Gene, y = log2FoldChange)) +
  ggplot2::geom_point(color = "#000000") +
  ggplot2::labs(title = "Ranked DEGs", x = "Gene", y = "Log2FoldChange") +
  ggplot2::theme_bw() +
  ggplot2::theme(panel.grid = element_blank())

plotData2 <- nmads_data(degData = degData)
plot2 <- nmads_plot(plot_data = plotData2, plot_title = "All DEGs Nmads")
degData2 <- dplyr::filter(degData, log2FoldChange >= 0)
plotData3 <- nmads_data(degData = degData2)
plot3 <- nmads_plot(plot_data = plotData3, plot_title = "Up DEGs Nmads")
degData3 <- dplyr::filter(degData, log2FoldChange < 0)
plotData4 <- nmads_data(degData = degData3)
plot4 <- nmads_plot(plot_data = plotData4, plot_title = "Down DEGs Nmads")

# 将多个图像保存在一个 PDF 文件里
plotList <- list(plot1, plot2,  plot3, plot4)
pdf(plotPath)
plotList
dev.off()

writeLines("o(^▽^)o")
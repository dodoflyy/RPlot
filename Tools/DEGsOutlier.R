# 根据提供的 Nmads 阈值提取属于 Outlier 的差异基因
# 可以用来选择特征基因
# 默认为 DESeq2 结果
# 脚本在 R 3.6 测试通过
# 依赖以下 R 包
# tidyverse

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))

scriptDescription <- "从 DESeq2 差异基因结果里根据提供的 Nmads 阈值提取 Outlier 基因"
parser <- ArgumentParser(description=scriptDescription, add_help=TRUE)
parser$add_argument("--DEGs", "-D", dest="DEGs", help="csv 格式 DESeq 差异基因结果", required=TRUE)
parser$add_argument("--Output", "-O", dest="OUT", help="输出路径", required=TRUE)
parser$add_argument("--UpNmads", "-UN", dest="UN", help="上调 DEGs 离群阈值", required=TRUE)
parser$add_argument("--DownNmads", "-DN", dest="DN", help="下调 DEGs 离群阈值", required=TRUE)

argvs <- parser$parse_args()
degPath <- file.path(argvs$DEGs)
outlierPath <- file.path(argvs$OUT)
upperCutoff <- as.numeric(argvs$UN)
lowerCutoff <- as.numeric(argvs$DN)

# 不用与中位值的差的绝对值
nmads_data <- function(deg_data) {
  degMad <- median(deg_data$log2FoldChange)
  nmadsData <- dplyr::mutate(deg_data, Nmads=log2FoldChange - degMad)
  return(nmadsData)
}

# 包含阈值本身
degData <- read_csv(degPath) %>% dplyr::filter(!is.na(log2FoldChange) & padj < 0.05)
upDeg <- dplyr::filter(degData, log2FoldChange >= 0)
downDeg <- dplyr::filter(degData, log2FoldChange < 0)
upNmads <- nmads_data(upDeg) %>% dplyr::mutate(Outlier = (Nmads >= upperCutoff))
downNmads <- nmads_data(downDeg) %>% dplyr::mutate(Outlier = (Nmads <= lowerCutoff))
nmadsData <- dplyr::bind_rows(upNmads, downNmads) %>% dplyr::filter(Outlier)
write_csv(nmadsData, outlierPath)

message("\n完成！")
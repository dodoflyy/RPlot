# 为病毒画深度覆盖图，用 geom_area 而不是 geom_line
# 输入数据有 Position, Depth 2列分别是位置和深度
# 脚本在 R3.6 测试通过
# 需要下列 R 包依赖
# argparse, tidyverse

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))

scriptDescription <- "测序深度图"
parser <- ArgumentParser(description=scriptDescription, add_help=TRUE)
parser$add_argument("--Input", "-I", dest="INPUT", help="csv 格式每个位置深度数据", required=TRUE)
parser$add_argument("--OutputDir", "-O", dest="OUT", help="输出目录", required=TRUE)
parser$add_argument("--Basename", "-B", dest="BASE", help="输出文件名", required=TRUE)
parser$add_argument("--Title", "-T", dest="TITLE", help="图像标题", default="Read Depth")

argvs <- parser$parse_args()
depthPath <- file.path(argvs$INPUT)
outDir <- file.path(argvs$OUT)
baseName <- argvs$BASE
titleText <- argvs$TITLE

writeLines("\n====== 读取深度数据 ======")
depthData <- read_csv(depthPath)
p <- ggplot2::ggplot(data = depthData, aes(Position, Depth)) +
  ggplot2::geom_area(fill="black") +
  ggplot2::labs(title = titleText, x = "Position", y = "Depth") +
  ggplot2::scale_x_continuous(labels = unit_format(scale = 0.001, unit = "K")) +
  ggplot2::theme_bw()

p10 <- ggplot2::ggplot(data = depthData, aes(Position, log10(Depth + 1))) +
  ggplot2::geom_area(fill="black") +
  ggplot2::labs(title = titleText, x = "Position", y = "log10(Depth + 1)") +
  ggplot2::scale_x_continuous(labels = unit_format(scale = 0.001, unit = "K")) +
  ggplot2::theme_bw()

pdfName1 <- paste(baseName, "Depth", "pdf", sep=".")
pdfName2 <- paste(baseName, "Depth", "Log10", "pdf", sep=".")
epsName1 <- paste(baseName, "Depth", "eps", sep=".")
epsName2 <- paste(baseName, "Depth", "Log10", "eps", sep=".")

writeLines("\n====== 保存图像 ======")
ggplot2::ggsave(filename = pdfName1, plot = p, device = "pdf", path = outDir, width = 200, height = 150, units = "mm")
ggplot2::ggsave(filename = epsName1, plot = p, device = "eps", path = outDir, width = 200, height = 150, units = "mm")
ggplot2::ggsave(filename = pdfName2, plot = p10, device = "pdf", path = outDir, width = 200, height = 150, units = "mm")
ggplot2::ggsave(filename = epsName2, plot = p10, device = "eps", path = outDir, width = 200, height = 150, units = "mm")

writeLines("\n完成！")
# WGBS M-bias plot
# 输入为 csv 格式，输出为 pdf
# 脚本在 R3.6 环境测试通过
# 需要的 R 包依赖：
# tidyverse, gridExtra, argparse

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gridExtra))

scriptDescription <- "WGBS 测序的 M-bias plot"
parser <- ArgumentParser(description=scriptDescription, add_help=TRUE)
parser$add_argument("--input", dest="INPUT", help="csv 格式输入数据", required=TRUE)
parser$add_argument("--output", dest="OUTPUT", help="输出 pdf 路径", required=TRUE)

argvs <- parser$parse_args()
inPath <- file.path(argvs$INPUT)
outPath <- file.path(argvs$OUTPUT)

plotData <- read_csv(inPath)
colFun <- c("#FB0007", "#139177", "#ED9E08", "#000000", "#169F0F")
p1 <- ggplot(plotData) +
  geom_line(aes(Position, Methylated_Percent, color=Context), size = 2, show.legend = FALSE) +
  scale_color_manual(values=colFun) +
  labs(title = "M-bias", x="Position", y="Methylated Percent") +
  theme_bw(base_size = 13) +
  facet_wrap(vars(Read), ncol = 1)

p2 <- ggplot(plotData) +
  geom_line(aes(Position, Total_Count, color=Context), size = 2) +
  scale_color_manual(values=colFun) +
  labs(title = "M-bias", x="Position", y="Total Count") +
  theme_bw(base_size = 13) +
  facet_wrap(vars(Read), ncol = 1)

p <- grid.arrange(p1, p2, nrow=1)
ggsave(filename = outPath, plot = p, device = "pdf", width = 800, height = 400, units = "mm")
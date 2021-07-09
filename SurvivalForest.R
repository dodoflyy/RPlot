# 从 cox 分析结果画森林图
# 脚本在 R 3.6 环境测试通过
# 需要以下包支持
# tidyverse, argparse

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))

scriptDescription <- "Cox 生存分析的森林图"
parser <- ArgumentParser(description=scriptDescription, add_help=TRUE)
parser$add_argument("--input", dest="COX", help="csv 格式 COX 分析结果", required=TRUE)
parser$add_argument("--output", dest="OUT", help="输出 pdf 图像路径", required=TRUE)
parser$add_argument("--y-title", dest="YLAB", help="Y 轴标题", required=TRUE)
parser$add_argument("--color-low", dest="COL1", help="P 值 0 对应的颜色，默认 \"#f03b20\"", default="#f03b20")
parser$add_argument("--color-mid", dest="COL2", help="P 值 0.5 对应颜色，默认 \"#f03b20\"", default="#ffeda0")
parser$add_argument("--color-high", dest="COL3", help="P 值 1 对应颜色，默认 \"#2b8cbe\"", default="#2b8cbe")

argvs <- parser$parse_args()
inPath <- file.path(argvs$COX)
outPath <- file.path(argvs$OUT)z
yTitle <- argvs$YLAB
col1 <- argvs$COL1
col2 <- argvs$COL2
col3 <- argvs$COL3

coxData <- read_csv(inPath)
p <- ggplot(coxData) +
  geom_point(aes(x = HR, y = term, colour = P, fill = P), shape = 23, stroke = 3) +
  geom_segment(aes(x = lower_95, xend = upper_95, y = term, yend = term, colour = P), size = 2) +
  geom_vline(xintercept = 1, linetype = "longdash", size = 1.5) +
  scale_colour_gradient2(low = col1, mid = col2, high = col3, breaks = c(0, 0.5, 1)) +
  scale_fill_gradient2(low = col1, mid = col2, high = col3, breaks = c(0, 0.5, 1)) +
  xlab("Hazard ratios") + ylab(yTitle) + 
  labs(colour = "P value", fill = "P value") +
  theme(axis.text = element_text(size = 12), panel.border = element_rect(fill = NA, linetype = "solid"), 
        panel.grid = element_blank(), panel.background = element_rect(fill = "white"), 
        title = element_text(size = 12, face = "bold"))

ggsave(filename = outPath, plot = p, device = "pdf")

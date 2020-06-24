# WGBS M-bias plot
# 输入为 csv 格式，输出为 pdf
# 脚本在 R3.6 环境测试通过
# 需要的 R 包依赖：tidyverse, gridExtra

writeLines("\nRscript MbiasPlot.R Input.csv Output.pdf \n")
argvs <- commandArgs(trailingOnly = TRUE)
stopifnot(length(argvs) >= 2)

library(tidyverse, quietly = TRUE)
library(gridExtra, quietly = TRUE)
inPath <- argvs[1]
outPath <- argvs[2]
plotData <- read_csv(inPath)
colFun <- c("#FB0007", "#139177", "#ED9E08", "#000000", "#169F0F")
p1 <- ggplot2::ggplot(plotData) +
  ggplot2::geom_line(aes(Position, Methylated_Percent, color=Context), show.legend = FALSE) +
  ggplot2::scale_color_manual(values=colFun) +
  ggplot2::labs(title = "M-bias", x="Position", y="Methylated Percent") +
  ggplot2::theme_bw() +
  ggplot2::facet_wrap(vars(Read), ncol = 1)

p2 <- ggplot2::ggplot(plotData) +
  ggplot2::geom_line(aes(Position, Total_Count, color=Context)) +
  ggplot2::scale_color_manual(values=colFun) +
  ggplot2::labs(title = "M-bias", x="Position", y="Total Count") +
  ggplot2::theme_bw() +
  ggplot2::facet_wrap(vars(Read), ncol = 1)

p <- grid.arrange(p1, p2, nrow=1)
ggplot2::ggsave(filename = outPath, plot = p, device = "pdf", width = 800, height = 400, units = "mm")
writeLines("完成")
# 为病毒画深度覆盖图，用 geom_area 而不是 geom_line
# 输入数据有 Position, Depth 2列分别是位置和深度
# 脚本在 R3.6 测试通过
# 需要 tidyverse

writeLines("Rscript DepthPlot.R Depth.csv OutputDir Filename [\"PlotTitle\"]\n")
argvs <- commandArgs(trailingOnly = TRUE)
stopifnot(length(argvs) >= 3)

library(tidyverse, quietly = TRUE)
library(scales, quietly = TRUE)

if (length(argvs) >= 4) {
  titleText <- argvs[4]
} else {
  titleText <- "Read Depth"
}

Data <- read_csv(argvs[1])
p <- ggplot2::ggplot(data = Data, aes(Position, Depth)) +
  ggplot2::geom_area(fill="black") +
  ggplot2::labs(title = titleText, x = "Position", y = "Depth") +
  ggplot2::scale_x_continuous(labels = unit_format(scale = 0.001, unit = "K")) +
  ggplot2::theme_bw()

p10 <- ggplot2::ggplot(data = Data, aes(Position, log10(Depth + 1))) +
  ggplot2::geom_area(fill="black") +
  ggplot2::labs(title = titleText, x = "Position", y = "log10(Depth + 1)") +
  ggplot2::scale_x_continuous(labels = unit_format(scale = 0.001, unit = "K")) +
  ggplot2::theme_bw()

fileName <- argvs[3]
pdfOut <- stringr::str_glue("{fileName}.pdf")
epsOut <- stringr::str_glue("{fileName}.eps")
pdfOut10 <- stringr::str_glue("{fileName}_Log10.pdf")
epsOut10 <- stringr::str_glue("{fileName}_Log10.eps")

outDir <- argvs[2]
ggplot2::ggsave(filename = pdfOut, plot = p, device = "pdf", path = outDir, width = 200, height = 150, units = "mm")
ggplot2::ggsave(filename = epsOut, plot = p, device = "eps", path = outDir, width = 200, height = 150, units = "mm")
ggplot2::ggsave(filename = pdfOut10, plot = p10, device = "pdf", path = outDir, width = 200, height = 150, units = "mm")
ggplot2::ggsave(filename = epsOut10, plot = p10, device = "eps", path = outDir, width = 200, height = 150, units = "mm")

writeLines("完成")
# 为病毒画深度覆盖图，用 geom_area 而不是 geom_line
# 输入数据有 Position, Depth 2列分别是位置和深度
# 脚本在 R3.6 测试通过
# 需要 tidyverse

writeLines("\n测序深度图\n")
writeLines("Usage:")
writeLines("Rscript DepthPlot.R Depth.csv OutputDir Filename [\"PlotTitle\"]\n")
argvs <- commandArgs(trailingOnly = TRUE)
stopifnot(length(argvs) >= 3)

depthPath <- file.path(argvs[1])
outDir <- argvs[2]
fileName <- argvs[3]
writeLines(stringr::str_glue("深度数据：{depthPath}"))
writeLines(stringr::str_glue("输出目录：{outDir}"))
writeLines(stringr::str_glue("输出文件名：{fileName}"))
if (length(argvs) >= 4) {
  titleText <- argvs[4]
  writeLines(stringr::str_glue("图像标题：{titleText}"))
} else {
  titleText <- "Read Depth"
}

library(tidyverse, quietly = TRUE, warn.conflicts = FALSE)
library(scales, quietly = TRUE)

Data <- read_csv(depthPath)
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

pdfOut <- stringr::str_glue("{fileName}.pdf")
epsOut <- stringr::str_glue("{fileName}.eps")
pdfOut10 <- stringr::str_glue("{fileName}_Log10.pdf")
epsOut10 <- stringr::str_glue("{fileName}_Log10.eps")

ggplot2::ggsave(filename = pdfOut, plot = p, device = "pdf", path = outDir, width = 200, height = 150, units = "mm")
ggplot2::ggsave(filename = epsOut, plot = p, device = "eps", path = outDir, width = 200, height = 150, units = "mm")
ggplot2::ggsave(filename = pdfOut10, plot = p10, device = "pdf", path = outDir, width = 200, height = 150, units = "mm")
ggplot2::ggsave(filename = epsOut10, plot = p10, device = "eps", path = outDir, width = 200, height = 150, units = "mm")

writeLines("\n完成")
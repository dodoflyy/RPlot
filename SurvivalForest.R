# 从 cox 分析结果画森林图
# 脚本在 R 3.6 环境测试通过
# 需要以下包支持
# tidyverse, BuenColors

writeLines("\nCox 生存分析的森林图\n")
writeLines("Usage:")
writeLines("Rscript SurvivalForest.R Input.csv OutputDir Filename AxisYLable\n")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 4)

inPath <- file.path(args[1])
outDir <- args[2]
fileName <- args[3]
yTitle <- args[4]

writeLines(stringr::str_glue("输入文件：{inPath}"))
writeLines(stringr::str_glue("输出目录：{outDir}"))
writeLines(stringr::str_glue("输出文件名：{fileName}"))
writeLines(stringr::str_glue("Y 轴标题：{yTitle}"))

library(tidyverse, quietly = TRUE, warn.conflicts = FALSE)
library(BuenColors, quietly = TRUE, warn.conflicts = FALSE)

coxData <- read_csv(inPath)
head(coxData) %>% print()

p <- ggplot(coxData) +
  geom_point(aes(x = HR, y = term, colour = P, fill = P), shape = 23, stroke = 3) +
  geom_segment(aes(x = lower_95, xend = upper_95, y = term, yend = term, colour = P), size = 2) +
  geom_vline(xintercept = 1, linetype = "longdash", size = 1.5) +
  scale_colour_gradientn(colours = rev(jdb_palette("Zissou"))) +
  scale_fill_gradientn(colours = rev(jdb_palette("Zissou"))) +
  xlab("Hazard ratios") + ylab(yTitle) + 
  labs(colour = "P value", fill = "P value") +
  theme(axis.text = element_text(size = 12), panel.border = element_rect(fill = NA, linetype = "solid"), panel.grid = element_blank(), panel.background = element_rect(fill = "white"), title = element_text(size = 12, face = "bold"))

pngPlot <- str_glue("{outDir}/{fileName}.png")
pdfPlot <- str_glue("{outDir}/{fileName}.pdf")
ggsave(filename = pdfPlot, plot = p, device = "pdf")
ggsave(filename = pngPlot, plot = p, dpi = 600, device = "png")

writeLines("\N完成")

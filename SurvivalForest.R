library(tidyverse, quietly = TRUE)
library(BuenColors, quietly = TRUE)

print("Rscript SurvivalForest.R 输入文件路径 输出目录 文件名 Y轴标题")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 4)
inPath <- args[1]
outDir <- args[2]
fileName <- args[3]
yTitle <- args[4]

coxData <- read_tsv(inPath)
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
ggsave(filename = pdfPlot, plot = p)
ggsave(filename = pngPlot, plot = p, dpi = 600)
print("完成")

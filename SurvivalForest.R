# 从 cox 分析结果画森林图
# 脚本在 R 3.6 环境测试通过
# 需要以下包支持
# tidyverse, BuenColors

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(BuenColors))

scriptDescription <- "Cox 生存分析的森林图"
parser <- ArgumentParser(description=scriptDescription, add_help=TRUE)
parser$add_argument("--Cox", "-C", dest="COX", help="csv 格式 COX 分析结果", required=TRUE)
parser$add_argument("--OutputDir", "-O", dest="OUT", help="输出目录", required=TRUE)
parser$add_argument("--Basename", "-B", dest="BASE", help="输出文件名", required=TRUE)
parser$add_argument("--Ylab", "-Y", dest="YLAB", help="Y 轴标题", required=TRUE)

argvs <- parser$parse_args()
inPath <- file.path(argvs$COX)
outDir <- file.path(argvs$OUT)
baseName <- argvs$BASE
yTitle <- argvs$YLAB

writeLines("\n====== 读取数据 ======")
coxData <- read_csv(inPath)
p <- ggplot(coxData) +
  geom_point(aes(x = HR, y = term, colour = P, fill = P), shape = 23, stroke = 3) +
  geom_segment(aes(x = lower_95, xend = upper_95, y = term, yend = term, colour = P), size = 2) +
  geom_vline(xintercept = 1, linetype = "longdash", size = 1.5) +
  scale_colour_gradientn(colours = rev(jdb_palette("Zissou"))) +
  scale_fill_gradientn(colours = rev(jdb_palette("Zissou"))) +
  xlab("Hazard ratios") + ylab(yTitle) + 
  labs(colour = "P value", fill = "P value") +
  theme(axis.text = element_text(size = 12), panel.border = element_rect(fill = NA, linetype = "solid"), 
        panel.grid = element_blank(), panel.background = element_rect(fill = "white"), 
        title = element_text(size = 12, face = "bold"))

writeLines("\n====== 保存图像 ======")
pngName <- paste(baseName, "Forest", "png", sep=".")
pdfName <- paste(baseName, "Forest", "pdf", sep=".")
ggsave(filename = pdfName, plot = p, device = "pdf", path=outDir)
ggsave(filename = pngName, plot = p, dpi = 600, device = "png", path=outDir)

writeLines("\N完成！")

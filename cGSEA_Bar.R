# 取 clusterProfiler GSEA 结果画条形图展示，默认展示 P < 0.05 的通路
# 默认输出2图片，分别使用通路 ID 和 Description
# 脚本在 R 3.6 环境测试通过
# 需要的依赖包
# argparse, tidyverse

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))

scriptDescription <- "对 clusterProfiler GSEA 分析结果画条形图"
parser <- ArgumentParser(description=scriptDescription, add_help=TRUE)
parser$add_argument("--GSEA", dest="GSEA", help="csv 格式 GSEA 结果", required=TRUE)
parser$add_argument("--OutputDir", dest="OUT", help="输出目录，默认当前目录", default=".")
parser$add_argument("--Basename", dest="BASE", help="输出文件名，默认 \"GSEA\"", default="GSEA")
parser$add_argument("--Title", dest="TITLE", help="图像标题，默认 \"GSEA\"", default="GSEA")
parser$add_argument("--Pvalue", dest="PVAL", help="P 值阈值，默认 0.25", default=0.25)
parser$add_argument("--Color1", dest="COL1", help="最大 P 值颜色，默认 \"#F9886B\"", default="#F9886B")
parser$add_argument("--Color2", dest="COL2", help="最小 P 值颜色，默认 \"#9E0013\"", default="#9E0013")

argvs <- parser$parse_args()
gseaPath <- file.path(argvs$GSEA)
outDir <- file.path(argvs$OUT)
baseName <- argvs$BASE
plotTitle <- argvs$TITLE
pVal <- as.double(argvs$PVAL)
color1 <- argvs$COL1
color2 <- argvs$COL2

plotData <- read_csv(gseaPath) %>% dplyr::filter(`p.adjust` < pVal) %>% dplyr::arrange(desc(enrichmentScore)) %>% 
  mutate(Pathway=factor(Description, levels=Description))
if (nrow(plotData) < 5) {
  writeLines("X﹏X")
  writeLines("通路数目太少，请检查数据和参数！")
  q(save="no")
}

labelWrap <- function(x) {
  return(str_wrap(x, width = 35))
}

barPlot <- ggplot(plotData, aes(Pathway, enrichmentScore)) +
  geom_bar(aes(fill = `p.adjust`), stat = "identity") +
  scale_fill_gradient(high=color1, low=color2) +
  labs(y="Enrichment Score", title=plotTitle, x="Pathway Name", fill="P value") +
  scale_x_discrete(labels = labelWrap) +
  theme(panel.background=element_rect(fill="#FFFFFF", color="#000000", linetype="solid"), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  coord_flip()

pdfName <- paste(baseName, "Barplot", "pdf", sep=".")
plotHeight <- nrow(plotData) * 6
if (plotHeight < 100) {
  plotHeight = 100
}
ggsave(filename = pdfName, plot = barPlot, device = "pdf", path = outDir, width = 220, 
       height = plotHeight, units = "mm", limitsize = FALSE)

writeLines("\n╰(*°▽°*)╯")
# 从 DESeq2 差异基因结果生成火山图
# 在 R 3.6 环境测试通过
# 需要 tidyverse 包支持

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))

scriptDescription <- "基因差异表达火山图"
parser <- ArgumentParser(description=scriptDescription, add_help=TRUE)
parser$add_argument("--DEGs", "-D", dest="DEGs", help="DESeq2 差异基因分析结果", required=TRUE)
parser$add_argument("--OutputDir", "-O", dest="OUT", help="输出目录", required=TRUE)
parser$add_argument("--Basename", "-B", dest="BASE", help="输出文件名", required=TRUE)
parser$add_argument("--Title", "-T", dest="TITLE", help="图像标题", default="DEGs")
parser$add_argument("--Pvalue", "-P", dest="PVAL", help="P 值阈值", default=0.05)
parser$add_argument("--Ratio", "-R", dest="RATIO", help="差异倍数绝对值阈值", default=1)

argvs <- parser$parse_args()
inPath <- file.path(argvs$DEGs)
outDir <- file.path(argvs$OUT)
baseName <- argvs$BASE
titleText <- argvs$TITLE
pVal <- as.double(argvs$PVAL)
ratioCut <- as.double(argvs$RATIO)

writeLines("\n====== 读取数据 ======")
degsData <- read_csv(inPath) %>% filter(degs, !is.na(padj))
vp <- ggplot(degsData, aes(log2FoldChange, -log10(padj))) + 
      geom_point(aes(colour = (abs(log2FoldChange) >= ratioCut & padj < pVal)), alpha = 0.5, show.legend = FALSE) + 
	  scale_colour_manual(values = c("TRUE" = "#CD6155", "FALSE" = "#566573")) + 
	  geom_vline(xintercept = c(-ratioCut, ratioCut), alpha = 0.8, linetype = "dashed") + 
	  geom_hline(yintercept = -log10(pVal), alpha = 0.8, linetype = "dashed") + 
	  labs(x = "log2FoldChange", y = "-log10(P.adj)", title = titleText) + 
	  theme_bw() + 
	  theme(panel.grid = element_blank())

writeLines("\n====== 保存图像 ======")
pdfName <- paste(baseName, "Volcano", "pdf", sep=".")
pngName <- paste(baseName, "Volcano", "png", sep=".")
ggsave(filename = pdfName, plot = vp, device = "pdf", path=outDir)
ggsave(filename = pngName, plot = vp, dpi = 600, device = "png", path=outDir)

writeLines("\no(^▽^)o")



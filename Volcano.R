# 从 DESeq2 差异基因结果生成火山图
# 在 R 3.6 环境测试通过
# 需要 tidyverse 包支持


writeLines("Rscript Volcano.R DEG.csv OutputDir Filename \"PlotTitle\"\n")
args <- commandArgs(TRUE)
stopifnot(length(args) >= 4)

library(tidyverse)
in_path <- args[1]
out_path <- args[2]
name <- args[3]
titleText <- args[4]

degs <- read_csv(in_path)
stopifnot("log2FoldChange" %in% colnames(degs), "padj" %in% colnames(degs))
degs2 <- filter(degs, !is.na(padj))
vp <- ggplot(degs2, aes(log2FoldChange, -log10(padj))) + 
      geom_point(aes(colour = (abs(log2FoldChange) >= 1 & padj < 0.01)), alpha = 0.5, show.legend = FALSE) + 
	  scale_colour_manual(values = c("TRUE" = "#CD6155", "FALSE" = "#566573")) + 
	  geom_vline(xintercept = c(-log2(2), log2(2)), alpha = 0.8, linetype = "dashed") + 
	  geom_hline(yintercept = -log10(0.01), alpha = 0.8, linetype = "dashed") + 
	  labs(x = "log2FoldChange", y = "-log10(P.adj)", title = titleText) + 
	  theme_bw() + 
	  theme(panel.grid = element_blank())

pdfOut <- str_glue("{out_path}/{name}.pdf")
pngOut <- str_glue("{out_path}/{name}.png")
ggsave(filename = pdfOut, plot = vp, device = "pdf")
ggsave(filename = pngOut, plot = vp, dpi = 600, device = "png")
writeLines("完成")



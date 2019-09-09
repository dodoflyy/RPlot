library(tidyverse)

args <- commandArgs(TRUE)
stopifnot(length(args) >= 4)
in_path <- args[1]
out_path <- args[2]
group2 <- args[4]
group1 <- args[3]
titleText <- paste(group2, group1, sep = " VS ")

degs <- read_tsv(in_path)
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

pdfOut <- paste(out_path, "pdf", sep = ".")
epsOut <- paste(out_path, "eps", sep = ".")
pngOut <- paste(out_path, "png", sep = ".")
ggsave(filename = pdfOut, plot = vp)
ggsave(filename = epsOut, plot = vp)
ggsave(filename = pngOut, plot = vp)
print("完成")



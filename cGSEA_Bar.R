# 取 clusterProfiler GSEA 结果画条形图展示，默认展示 P < 0.05 的通路
# 脚本在 R 3.6 环境测试通过
# 需要 tidyverse 包支持

writeLines("Rscript cGSEA_Bar.R GSEA.csv OutputDir Filename [\"Plot Title\"]\n")
argvs <- commandArgs(trailingOnly = TRUE)
stopifnot(length(argvs) >= 3)

library(tidyverse, quietly = TRUE)
plotData <- read_csv(argvs[1]) %>% dplyr::filter(`p.adjust` < 0.05) %>% dplyr::arrange(desc(enrichmentScore))
plotData$Description <- factor(plotData$Description, levels = plotData$Description)
stopifnot(length(ncol(plotData)) >= 1)

if (length(argvs) >= 4) {
  plotTitle <- argvs[4]
} else {
  plotTitle <- "Pathway GSEA"
}

p <- ggplot2::ggplot(plotData, aes(Description, enrichmentScore)) +
  ggplot2::geom_bar(aes(fill = `p.adjust`), stat = "identity") +
  scale_fill_gradient(high="#1873CC", low="#ED2C2C") +
  labs(y="Enrichment Score", title=plotTitle, x="Pathway Name", fill="P value") +
  theme(panel.background=element_rect(fill="#FFFFFF", color="#000000", linetype="solid"), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  coord_flip()

filename <- argvs[3]
pdfName <- stringr::str_glue("{filename}.pdf")
pngName <- stringr::str_glue("{filename}.png")
ggplot2::ggsave(filename = pdfName, plot = p, device = "pdf", path = argvs[2])
ggplot2::ggsave(filename = pngName, plot = p, device = "png", path = argvs[2], dpi = 600)

writeLines("完成")
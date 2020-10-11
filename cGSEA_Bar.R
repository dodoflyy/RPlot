# 取 clusterProfiler GSEA 结果画条形图展示，默认展示 P < 0.05 的通路
# 默认输出2图片，分别使用通路 ID 和 Description
# 脚本在 R 3.6 环境测试通过
# 需要 tidyverse 包支持

writeLines("\n对 clusterProfiler GSEA 分析结果画条形图\n")
writeLines("Usage:")
writeLines("Rscript cGSEA_Bar.R GSEA.csv OutputDir Filename [\"Plot Title\"]\n")
argvs <- commandArgs(trailingOnly = TRUE)
stopifnot(length(argvs) >= 3)
gseaPath <- file.path(argvs[1])
outDir <- file.path(argvs[2])
fileName <- argvs[3]
writeLines(stringr::str_glue("GSEA 数据：{gseaPath}"))
writeLines(stringr::str_glue("输出目录：{outDir}"))
writeLines(stringr::str_glue("输出文件名：{fileName}"))
if (length(argvs) >= 4) {
  plotTitle <- argvs[4]
  writeLines(stringr::str_glue("图像标题：{plotTitle}"))
} else {
  plotTitle <- "Pathway GSEA"
}

library(tidyverse, quietly = TRUE, warn.conflicts = FALSE)

plotData <- read_csv(gseaPath) %>% dplyr::filter(`p.adjust` < 0.05) %>% dplyr::arrange(desc(enrichmentScore))
plotData$Description <- factor(plotData$Description, levels = plotData$Description)
stopifnot(length(ncol(plotData)) >= 1)

labelWrap <- function(x) {
  return(str_wrap(x, width = 45))
}

p1 <- ggplot2::ggplot(plotData, aes(Description, enrichmentScore)) +
  ggplot2::geom_bar(aes(fill = `p.adjust`), stat = "identity") +
  scale_fill_gradient(high="#1873CC", low="#ED2C2C") +
  labs(y="Enrichment Score", title=plotTitle, x="Pathway Name", fill="P value") +
  ggplot2::scale_x_discrete(labels = labelWrap) +
  theme(panel.background=element_rect(fill="#FFFFFF", color="#000000", linetype="solid"), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  coord_flip()

plotData$ID<- factor(plotData$ID, levels = plotData$ID)
p2 <- ggplot2::ggplot(plotData, aes(ID, enrichmentScore)) +
  ggplot2::geom_bar(aes(fill = `p.adjust`), stat = "identity") +
  scale_fill_gradient(high="#1873CC", low="#ED2C2C") +
  labs(y="Enrichment Score", title=plotTitle, x="Pathway Name", fill="P value") +
  theme(panel.background=element_rect(fill="#FFFFFF", color="#000000", linetype="solid"), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  coord_flip()

pdfName1 <- stringr::str_glue("{fileName}_Bar1.pdf")
plotHeight1 <- nrow(plotData) * 6
if (plotHeight1 < 100) {
  plotHeight1 = 100
}
ggplot2::ggsave(filename = pdfName1, plot = p1, device = "pdf", path = outDir, width = 260, height = plotHeight1, units = "mm", limitsize = FALSE)

pdfName2 <- stringr::str_glue("{fileName}_Bar2.pdf")
plotHeight2 <- nrow(plotData) * 2.7
if (plotHeight2 < 100) {
  plotHeight2 = 100
}
ggplot2::ggsave(filename = pdfName2, plot = p2, device = "pdf", path = outDir, width = 200, height = plotHeight2, units = "mm", limitsize = FALSE)

writeLines("\n完成！")
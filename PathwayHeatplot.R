# 行为通路，列为基因的通路富集热图
# 固定展示 P < 0.05 的通路，基因展示前 50 个
# 默认差异基因是 DESeq2 的结果，如果不是要修改表头
# 默认 GSEA 结果是 clusterProfiler 产生的
# 在 R 3.6 环境测试通过
# 需要以下 R 包支持
# tidyverse, ComplexHeatmap, BuenColors


writeLines("Rscript PathwayHeatplot.R Pathway.csv DEG.csv Heatplot.pdf \"TitleText\" \n")
argvs <- commandArgs(trailingOnly = TRUE)
stopifnot(length(argvs) >= 4)

options(stringsAsFactors = FALSE)
library(tidyverse, quietly = TRUE)
library(ComplexHeatmap, quietly = TRUE)
library(BuenColors, quietly=TRUE)

# 默认用 entrezgene_id 进行 GSEA
degData <- read_csv(argvs[2]) %>% dplyr::filter(!is.na(entrezgene_id)) %>% dplyr::distinct(entrezgene_id, .keep_all=TRUE)
deg <- degData$log2FoldChange
names(deg) <- degData$entrezgene_id

pathwayData <- read_csv(argvs[1]) %>% dplyr::filter(`p.adjust` < 0.05)
stopifnot(nrow(pathwayData) > 1)
geneList <- stringr::str_c(pathwayData$core_enrichment, collapse = "/") %>% strsplit(split = "/", fixed = TRUE) %>% unlist() %>% table() %>% unlist()
rankGeneList <- geneList[order(geneList, decreasing = TRUE)] %>% names()
if (length(rankGeneList) > 50) {
  rankGeneList <- rankGeneList[1:50]
}

DataList <- list()
pathwayGene <- pathwayData$core_enrichment
Description <- pathwayData$Description
names(pathwayGene) <- pathwayData$Description
for (i in 1:length(pathwayGene)) {
  # 按照 rankGeneList 的基因顺序找到相应基因差异倍数，如果相应基因不在 core enrichment 里，那么就赋值NA
  # 然后每一条通路的结果存到 list
  pathwayName <- Description[i]
  coreGene <- strsplit(pathwayGene[i], split = "/", fixed = TRUE) %>% unlist()
  geneMatch <- match(rankGeneList, coreGene)
  matchCore <- coreGene[geneMatch]
  geneDeg <- deg[match(matchCore, names(deg))]
  names(geneDeg) <- rankGeneList
  DataList[[pathwayName]] <- geneDeg
}
# head(DataList, n=3)
Data <- as.data.frame(DataList)
# head(Data, n=3)
geneMap <- tibble(entrezgene_id=as.double(rankGeneList)) %>% dplyr::left_join(degData, by="entrezgene_id") %>% dplyr::select(entrezgene_id, hgnc_symbol)
head(geneMap, n=3)
head(rankGeneList)
dim(geneMap)
length(rankGeneList)
rownames(Data) <- geneMap$hgnc_symbol
colnames(Data) <- Description
Data <- Data[, colSums(!is.na(Data)) > 0]

# 转换到行为通路，列为基因
PlotData <- t(Data)
head(PlotData, n=3)
color <- jdb_palette("brewer_celsius", type="continuous")
if (nrow(PlotData) <= 10) {
  h <- 4
  w <- 22
} else if (nrow(PlotData) > 10 && nrow(PlotData) <= 20) {
  h <- 7
  w <- 22
} else if (nrow(PlotData) > 20 && nrow(PlotData) <= 30) {
  h <- 9
  w <- 22
}  else if (nrow(PlotData) > 30 && nrow(PlotData) <= 50) {
  h <- 14
  w <- 22
} else {
  h <- 18
  w <- 22
}
 
hm <- Heatmap(PlotData, name="FoldChange(Log2)", col=color, column_title = argvs[4], column_title_side = "top", cluster_rows = FALSE, cluster_columns = FALSE, show_column_dend = FALSE, show_row_dend = FALSE, show_row_names = TRUE, show_column_names = TRUE, row_names_side = "left", column_names_side = "bottom", column_names_rot = 45, na_col = "white", column_order = order(colSums(is.na(PlotData))), row_order = order(rowSums(!is.na(PlotData))), row_names_max_width = max_text_width(rownames(PlotData)), border = TRUE, rect_gp = gpar(col="white"), heatmap_legend_param = list(grid_height = unit(8, "mm"), grid_width = unit(6, "mm")))

pdf(argvs[3], width = w, height = h)
draw(hm)
dev.off()

writeLines("完成")
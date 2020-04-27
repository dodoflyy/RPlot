# 行为通路，列为基因的通路富集热图
# 基因展示前 50 个
# 默认差异基因是 DESeq2 的结果，如果不是要修改表头
# 默认 GSEA 结果是 clusterProfiler 产生的
# 在 R 3.6 环境测试通过
# 需要以下 R 包支持
# tidyverse, ComplexHeatmap, BuenColors

writeLines("\nRscript PathwayHeatplot.R Pathway.csv CutoffP DEG.csv Heatplot.pdf \"TitleText\" \n")
argvs <- commandArgs(trailingOnly = TRUE)
stopifnot(length(argvs) >= 4)

options(stringsAsFactors = FALSE)
library(tidyverse, quietly = TRUE)
library(ComplexHeatmap, quietly = TRUE)
library(BuenColors, quietly=TRUE)
library(circlize, quietly = TRUE)

# 默认用 entrezgene_id 进行 GSEA
# 用 SYMBOL 进行排序，这样重复 entrezgene_id 时能优先选择有 symbol 的
degData <- read_csv(argvs[3]) %>% dplyr::arrange(hgnc_symbol) %>% dplyr::filter(!is.na(entrezgene_id)) %>% dplyr::distinct(entrezgene_id, .keep_all=TRUE)
deg <- degData$log2FoldChange
names(deg) <- degData$entrezgene_id

pCutoff <- as.double(argvs[2])
pathwayData <- read_csv(argvs[1]) %>% dplyr::filter(`p.adjust` < pCutoff)
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
head(rankGeneList)
glimpse(geneMap)
length(rankGeneList)
rownames(Data) <- geneMap$hgnc_symbol
colnames(Data) <- Description
Data <- Data[, colSums(!is.na(Data)) > 0]

# 转换到行为通路，列为基因
PlotData <- t(Data)
head(PlotData, n=3)
w <- 22
h <- nrow(PlotData) * 0.25
if (h < 7) {
  h = 7
}

minL <- min(PlotData, na.rm = TRUE)
maxL <- max(PlotData, na.rm = TRUE)
color <- jdb_palette("ocean_brick", type="continuous")
if (all(abs(minL) <= 1 , abs(maxL) <= 1)) {
  fromL <- -1
  toL <- 1
  breaksL <- c(-1, 0, 1)
  labelsL <- c("-1", "0", "1")
  heightL <- 8
} else if (all(abs(minL) <= 2 , abs(maxL) <= 2)) {
  fromL <- -2
  toL <- 2
  breaksL <- c(-2, 0, 2)
  labelsL <- c("-2", "0", "2")
  heightL <- 8
} else {
  fromL <- -4
  toL <- 4
  breaksL <- c(-4, -2, 0, 2, 4)
  labelsL <- c("<= -4", "-2", "0", "2", ">= 4")
  heightL <- 12
}

color_fun <- colorRamp2(seq(from = fromL, to = toL, length.out = length(color)), color)
hm <- Heatmap(PlotData, name="FoldChange(Log2)", col=color_fun, column_title = argvs[4], column_title_side = "top", cluster_rows = FALSE, cluster_columns = FALSE, show_column_dend = FALSE, show_row_dend = FALSE, show_row_names = TRUE, show_column_names = TRUE, row_names_side = "left", column_names_side = "bottom", column_names_rot = 45, na_col = "white", column_order = order(colSums(is.na(PlotData))), row_order = order(rowSums(!is.na(PlotData))), row_names_max_width = max_text_width(rownames(PlotData)), border = TRUE, rect_gp = gpar(col="white"), heatmap_legend_param = list(grid_height = unit(heightL, "mm"), grid_width = unit(6, "mm"), at = breaksL, labels = labelsL))

pdf(argvs[4], width = w, height = h)
draw(hm)
dev.off()

writeLines("完成")
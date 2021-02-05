# 从 DESeq2 的差异基因和 rlog 结果，生成图片展示样品间距离
# 差异基因结果指过滤后的数据，只包含差异基因
# 假定输入文件都包含 ensembl_gene_id, entrezgene_id, hgnc_symbol
# 将选定 ensembl_gene_id 作为基因 ID 进行分析
# 假定SampleGroup 有 Sample 和 Group 2列
# 在 R3.6 测试通过。
# 需要的依赖包
# argparse, tidyverse, BuenColors, ComplexHeatmap, ggrepel

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(BuenColors))
suppressPackageStartupMessages(library(ggrepel))

scriptDescription="出图展示 RNA-seq 样品分布"
parser <- ArgumentParser(description=scriptDescription, add_help=TRUE)
parser$add_argument("--Expression", "-E", dest="EXPR", help="csv 格式的表达数据", required=TRUE)
parser$add_argument("--Group", "-G", dest="GROUP", help="样品分组文件，包含 Sample, Group 2 列", required=TRUE)
parser$add_argument("--DEGs", "-D", dest="DEGs", help="csv 格式 DESeq2 差异基因分析结果", required=TRUE)
parser$add_argument("--OutputDir", "-O", dest="OUT", help="输出目录", required=TRUE)
parser$add_argument("--basename", "-B", dest="BASE", help="输出文件名", required=TRUE)

argvs <- parser$parse_args()
expressionPath <- file.path(argvs$EXPR)
groupPath <- file.path(argvs$GROUP)
outDir <- file.path(argvs$OUT)
baseName <- argvs$BASE
degsPath <- file.path(argvs$DEGs)

writeLines("\n====== 读取数据 ======")
sampleGroup <- read_csv(groupPath)
degsData <- read_csv(degsPath)
exprData <- read_csv(expressionPath) %>% dplyr::select(-entrezgene_id, -hgnc_symbol) %>% 
  dplyr::filter(ensembl_gene_id %in% degsData$ensembl_gene_id) %>% 
  dplyr::distinct(ensembl_gene_id, .keep_all = TRUE) %>% as.data.frame()
rownames(exprData) <- exprData$ensembl_gene_id
exprData$ensembl_gene_id <- NULL

writeLines("\n====== PCA 分析 ======")
dataPca <- prcomp(t(exprData))
pcaSummary <- summary(dataPca)
print(pcaSummary)
# 取得 PC1,PC2 解释占比
summaryTable <- pcaSummary$importance
pc1Por <- summaryTable[[2, 1]] %>% round(3)
pc2Por <- summaryTable[[2, 2]] %>% round(3)
pcaData <- dataPca$x %>% as_tibble(rownames="Sample") %>% dplyr::select(Sample, PC1, PC2) %>% 
  dplyr::left_join(sampleGroup, by="Sample")
glimpse(pcaData)
xLab <- stringr::str_glue("PC1({pc1Por})")
yLab <- stringr::str_glue("PC2({pc2Por})")
pcaPlot <- ggplot2::ggplot(pcaData, aes(x=PC1, y=PC2)) +
           geom_point(aes(shape=Group)) +
           ggrepel::geom_text_repel(aes(label=Sample)) +
           ggplot2::labs(title = "Sample PCA", x = xLab, y = yLab) +
           ggplot2::theme_bw()

writeLines("\n====== 聚类分析 ======")
dataDist <- dist(t(exprData), method = "euclidean")
dataHc <- hclust(dataDist, method = "average")

col_fun <- jdb_palette("brewer_fire", type = "continuous") %>% rev()
hmData <- as.matrix(dataDist)
hm <- Heatmap(hmData, name = "Euclidean", col = col_fun, cluster_rows = TRUE, 
              cluster_columns = TRUE, show_row_names = TRUE, row_names_side = "right", 
              show_row_dend = FALSE, show_column_dend = FALSE, column_title = "Sample Distance", 
              column_title_side = "top")

writeLines("\n====== 保存图片 ======")
pcaName <- paste(baseName, "PCA", "pdf", sep=".")
treeName <- paste(baseName, "clust", "pdf", sep=".")
hmName <- paste(baseName, "HeatPlot", "pdf", sep=".")
treePath <- file.path(outDir, treeName)
pcaPath <- file.path(outDir, pcaName)
ggplot2::ggsave(filename =  pcaName, plot = pcaPlot, device = "pdf", path = outDir)

pdf(treePath)
plot(dataHc, xlab = "Sample")
dev.off()

pdf(pcaPath)
draw(hm)
dev.off()


writeLines("\n完成！")
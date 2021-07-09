# 用 clusterProfiler 进行 GSEA 分析后一般结果文件 core_enrichment 列为 ENTREZ GENE ID
# 用这个脚本将 ENTREZ 修改为 SYMBOL 方便阅读
# 不同的基因 ID 采用 Ensembl biomart 数据库进行转换
# 脚本在 R 3.6 测试通过
# 需要下列包支持
# argparse, tidyverse, biomaRt

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(biomaRt))

scriptDescription <- "将 GSEA 结果的 core_enrichment 列从 ENTREZID 改为 SYMBOL"
parser <- ArgumentParser(description=scriptDescription, add_help=TRUE)
parser$add_argument("--input", dest="GSEA", help="csv 格式 GSEA 结果", required=TRUE)
parser$add_argument("--output", dest="OUT", help="输出路径", required=TRUE)

argvs <- parser$parse_args()
inPath <- file.path(argvs$GSEA)
outPath <- file.path(argvs$OUT)
ens <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

# 不考虑一个 entrezgene 对应多个 SYMBOL 情况，全部返回
mapGene <- function(entrez_vector) {
  geneMap <- getBM(filters = "entrezgene_id", attributes = c("entrezgene_id", "hgnc_symbol"), 
                  values = entrez_vector, mart = ens)
  symbolVector <- geneMap$hgnc_symbol %>% unique()
  return(symbolVector)
}

mapGeneString <- function(entrez_string) {
  geneList1 <- strsplit(entrez_string, split="/", fixed=TRUE) %>% unlist() %>% as.numeric()
  geneList2 <- mapGene(geneList1)
  symbolString <- paste(geneList2, collapse = "/")
  return(symbolString)
}

gseaData1 <- read_csv(inPath)
gseaData2 <- dplyr::mutate(gseaData1, core_enrichment_SYMBOL = purrr::map_chr(core_enrichment, mapGeneString))
readr::write_csv(gseaData2, outPath)

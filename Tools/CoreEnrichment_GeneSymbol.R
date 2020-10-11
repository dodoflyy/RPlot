# 用 clusterProfiler 进行 GSEA 分析后一般结果文件 core_enrichment 列为 ENTREZ GENE ID
# 用这个脚本将 ENTREZ 修改为 SYMBOL 方便阅读
# 不同的基因 ID 采用 Ensembl biomart 数据库进行转换
# 为了速度先将 biomart 数据库下载到本地，下载地址是 http://asia.ensembl.org/biomart/martview/a90bd86f40a8f019a194eee419189845
# 这个仓库里也提供了下载好了的 biomart 文件，格式是 csv
# 脚本在 R 3.6 测试通过
# 需要下列包支持
# tidyverse

writeLines("\n将 GSEA 结果的 core_enrichment 列从 ENTREZID 改为 SYMBOL\n")
writeLines("Usage:")
writeLines("Rscript CoreEnrichment_GeneSymbol.R GSEA.csv Biomart.csv Output.csv\n")

argvs <- commandArgs(trailingOnly = TRUE)
stopifnot(length(argvs) >= 3)

library(tidyverse, warn.conflicts = FALSE, quietly = TRUE)

inPath <- file.path(argvs[1])
biomartPath <- file.path(argvs[2])
outPath <- file.path(argvs[3])
writeLines(stringr::str_glue("输入 GSEA 结果：{inPath}"))
writeLines(stringr::str_glue("Biomart 数据库：{biomartPath}"))
writeLines(stringr::str_glue("结果输出：{outPath}"))

geneIDs <- readr::read_csv(biomartPath) %>% dplyr::select(`HGNC symbol`, `NCBI gene (formerly Entrezgene) ID`) %>% dplyr::rename(hgnc_symbol=`HGNC symbol`, entrezgene_id=`NCBI gene (formerly Entrezgene) ID`) %>% dplyr::filter(!is.na(entrezgene_id)) %>% dplyr::arrange(hgnc_symbol) %>% dplyr::distinct(entrezgene_id, .keep_all=TRUE)
dplyr::glimpse(geneIDs)

# 不考虑一个 entrezgene 对应多个 SYMBOL 情况，全部返回
# 如果出现 SYMBOL 为 NA 那么移除
mapGene <- function(entrezVector) {
  mapTable <- dplyr::filter(geneIDs, entrezgene_id %in% entrezVector) %>% dplyr::filter(!is.na(hgnc_symbol))
  symbolVector <- mapTable$hgnc_symbol
  return(symbolVector)
}

mapGeneString <- function(entrezString) {
  geneList1 <- stringr::str_split(entrezString, "/") %>% unlist() %>% as.numeric()
  geneList2 <- mapGene(geneList1)
  symbolString <- stringr::str_c(geneList2, collapse = "/")
  return(symbolString)
}

gseaData1 <- readr::read_csv(inPath)
gseaData2 <- dplyr::mutate(gseaData1, core_enrichment_SYMBOL = purrr::map_chr(core_enrichment, mapGeneString))
readr::write_csv(gseaData2, outPath)

writeLines("\n完成！")
# 从 clusterProfiler 通路富集结果输出泡泡图和柱状图
# 可以选择展示的通路数目，默认按P值排序后筛选相应数目
# 脚本在 R 3.6 测试通过
# 需要以下包依赖
# argparse, tidyverse, BuenColors


suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(BuenColors))
suppressPackageStartupMessages(library(tidyverse))

scriptDescription <- "生成 clusterProfiler 通路富集结果的泡泡图和柱状图"
parser <- ArgumentParser(description=scriptDescription, add_help=TRUE)
parser$add_argument("--Enrich", "-E", dest="ENRICH", help="csv 格式通路富集结果", required=TRUE)
parser$add_argument("--OutputDir", "-O", dest="OUT", help="输出目录，默认为当前目录", default=".")
parser$add_argument("--Basename", "-B", dest="BASE", help="输出文件名，默认 \"Enrichment\"", default="Enrichment")
parser$add_argument("--Title", "-T", dest="TITLE", help="图像标题，默认 \"Pathway Enrichment\"", default="Pathway Enrichment")
parser$add_argument("--Number", "-N", dest="NUM", help="展示的通路数目，默认 20", default=20)

argvs <- parser$parse_args()
inPath <- file.path(argvs$ENRICH)
outDir <- file.path(argvs$OUT)
baseName <- argvs$BASE
showNum <- as.integer(argvs$NUM) 
plotTitle <- argvs$TITLE

pathway <- read_csv(inPath) %>% arrange(`p.adjust`) %>% separate(GeneRatio, into=c("k", "n"), sep="/") %>% 
  separate(BgRatio, into=c("M", "N"), sep="/") %>% mutate(RichFactor=as.numeric(k)/as.numeric(M))
if (nrow(pathway) > showNum) {
  pathway <- dplyr::slice(pathway, 1:showNum)
}
glimpse(pathway)

dotPathway <- arrange(pathway, RichFactor)
dotPathway <- mutate(dotPathway, Description=factor(Description, levels=dotPathway$Description))
barPathway <- arrange(pathway, desc(Count))
barPathway <- mutate(barPathway, ID=factor(ID, levels=barPathway$ID))

brewerRed <- rev(jdb_palette("brewer_red"))
lowerCut <- as.integer(length(brewerRed) * 0.2)
upperCut <- as.integer(length(brewerRed) * 0.8)
brewerRed <- brewerRed[lowerCut:upperCut]

dotPlot <- ggplot(dotPathway, aes(x=Description, y=RichFactor)) +
  geom_point(aes(size=Count, colour=`p.adjust`)) +
  scale_colour_gradientn(colors=brewerRed) +
  scale_x_discrete(labels=scales::wrap_format(30)) +
  labs(y="RichFactor", title=plotTitle, size="Gene count", colour="P value") +
  theme(axis.title.y=element_blank(), panel.background=element_rect(fill="#FFFFFF", color="#000000", linetype="solid"), 
        panel.grid.major=element_line(color="#DCDCDC", linetype="solid")) +
  coord_flip()


barPlot <- ggplot(barPathway, aes(x=ID, y=Count)) +
  geom_bar(aes(fill=`p.adjust`), stat="identity") +
  scale_fill_gradientn(colors=brewerRed) +
  labs(y="Gene count", title=plotTitle, x="Pathway ID", fill="P value") +
  theme(axis.text.x=element_text(angle=45, hjust=1), panel.background=element_rect(fill="#FFFFFF", color="#000000", linetype="solid"), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank())

dotPng <- paste(baseName, "Dot", "png", sep = ".")
dotPdf <- paste(baseName, "Dot", "pdf", sep = ".")
barPng <- paste(baseName, "Bar", "png", sep = ".")
barPdf <- paste(baseName, "Bar", "pdf", sep = ".")

dotHeight <- showNum * 10 + 30
if (dotHeight < 150) {
  dotHeight<- 150
}
barWidth <- showNum * 10 + 20
if (barWidth < 120) {
  barWidth <- 120
}

ggsave(filename=dotPng, dpi=600, plot=dotPlot, device = "png", width = 150, height = dotHeight, unit = "mm", path=outDir)
ggsave(filename=dotPdf, plot=dotPlot, device = "pdf", width = 150, height = dotHeight, unit = "mm", path=outDir)
ggsave(filename=barPng, dpi=600, plot=barPlot, device = "png", width = barWidth, height = 150, unit = "mm", path=outDir)
ggsave(filename=barPdf, plot=barPlot, device = "pdf", width = barWidth, height = 150, unit = "mm", path=outDir)

writeLines("\nヽ(✿ﾟ▽ﾟ)ノ")


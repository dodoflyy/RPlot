# 从clusterProfiler 通路富集结果输出泡泡图和柱状图
# 可以选择展示的通路数目，默认按P值排序后筛选相应数目
# 脚本在 R 3.6 测试通过
# 需要以下包支持
# tidyverse, BuenColors

writeLines("\nRscript EnrichPathewy.R Input.csv OutputDir Filename PlotPathwayNum [PlotTitle]\n")
args <- commandArgs(TRUE)
stopifnot(length(args) >= 3)
inPath <- args[1]
outDir <- args[2]
fileName <- args[3]
showNum <- as.integer(args[4]) 

library(tidyverse)
library(BuenColors)

if(length(args) >= 5){
  plotTitle <- args[5]
  }else{
  plotTitle <- "Pathway enrichment"
  }

pathway <- read_csv(inPath) %>% arrange(`p.adjust`) %>% separate(GeneRatio, into=c("k", "n"), sep="/") %>% separate(BgRatio, into=c("M", "N"), sep="/") %>% mutate(RichFactor=as.numeric(k)/as.numeric(M))
if (nrow(pathway) > showNum) {
  pathway <- dplyr::slice(pathway, showNum)
}
glimpse(pathway)

dotPathway <- arrange(pathway, RichFactor)
dotPathway <- mutate(dotPathway, Description=factor(Description, levels=dotPathway$Description))
barPathway <- arrange(pathway, desc(Count))
barPathway <- mutate(barPathway, ID=factor(ID, levels=barPathway$ID))

brewerRed <- rev(jdb_palette("brewer_red"))
lowerCut <- as.integer(length(brewerRed) * 0.25)
upperCut <- as.integer(length(brewerRed) * 0.7)
brewerRed <- brewerRed[lowerCut:upperCut]

dotPlot <- ggplot(dotPathway, aes(x=Description, y=RichFactor)) +
  geom_point(aes(size=Count, colour=`p.adjust`)) +
  scale_colour_gradientn(colors=brewerRed) +
  scale_x_discrete(labels=scales::wrap_format(30)) +
  labs(y="RichFactor", title=plotTitle, size="Gene count", colour="P value") +
  theme(axis.title.y=element_blank(), panel.background=element_rect(fill="#FFFFFF", color="#000000", linetype="solid"), panel.grid.major=element_line(color="#DCDCDC", linetype="solid")) +
  coord_flip()

barPlot <- ggplot(barPathway, aes(x=ID, y=Count)) +
  geom_bar(aes(fill=`p.adjust`), stat="identity") +
  scale_fill_gradient(high="#539AC9", low="#0A3773") +
  labs(y="Gene count", title=plotTitle, x="Pathway ID", fill="P value") +
  theme(axis.text.x=element_text(angle=45, hjust=1), panel.background=element_rect(fill="#FFFFFF", color="#000000", linetype="solid"), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

dotPng <- str_glue("{outDir}/{fileName}_Dot.png")
dotPdf <- str_glue("{outDir}/{fileName}_Dot.pdf")
barPng <- str_glue("{outDir}/{fileName}_Bar.png")
barPdf <- str_glue("{outDir}/{fileName}_Bar.pdf")

dot_h <- showNum * 10 + 30
if (dot_h < 150) {
  dot_h <- 150
}
bar_w <- showNum * 10 + 20
if (bar_w < 120) {
  bar_w <- 120
}
ggsave(filename=dotPng, dpi=600, plot=dotPlot, device = "png", width = 150, height = dot_h, unit = "mm")
ggsave(filename=dotPdf, plot=dotPlot, device = "pdf", width = 150, height = dot_h, unit = "mm")
ggsave(filename=barPng, dpi=600, plot=barPlot, device = "png", width = bar_w, height = 150, unit = "mm")
ggsave(filename=barPdf, plot=barPlot, device = "pdf", width = bar_w, height = 150, unit = "mm")

writeLines("完成")


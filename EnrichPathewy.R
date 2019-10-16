library(tidyverse)
library(BuenColors)

print("Rscript EnrichPathewy.R 输入文件路径 输出目录 输出文件名 展示数目 [标题]")
args <- commandArgs(TRUE)
stopifnot(length(args) >= 3)
inPath <- args[1]
outDir <- args[2]
fileName <- args[3]
showNum <- args[4]

if(length(args) >= 5){
  plotTitle <- args[5]
  }else{
  plotTitle <- "Pathway enrichment"
  }

pathway <- read_tsv(inPath) %>% arrange(`p.adjust`) %>% separate(GeneRatio, into=c("k", "n"), sep="/") %>% separate(BgRatio, into=c("M", "N"), sep="/") %>% mutate(RichFactor=as.numeric(k)/as.numeric(M)) %>% slice(1:showNum)
head(pathway) %>% print()

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

dotPng <- str_glue("{outDir}/{fileName}-Dot.png")
dotPdf <- str_glue("{outDir}/{fileName}-Dot.pdf")
barPng <- str_glue("{outDir}/{fileName}-Bar.png")
barPdf <- str_glue("{outDir}/{fileName}-Bar.pdf")

ggsave(filename=dotPng, dpi=600, plot=dotPlot)
ggsave(filename=dotPdf, dpi=600, plot=dotPlot)
ggsave(filename=barPng, dpi=600, plot=barPlot)
ggsave(filename=barPdf, dpi=600, plot=barPlot)

print("完成")


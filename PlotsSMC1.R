install.packages('DESeq2')
install.packages('reshape2')
install.packages("vqv/ggbiplot")
ininstall.packages('plotly')
install.packages("devtools")
install.packages('ggfortify')
install.packages("Rcpp", dependencies = TRUE)
install.packages("dplyr")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

BiocManager::install("DESeq")
browseVignettes("DESeq")

BiocManager::install('EnhancedVolcano')


library(EnhancedVolcano)
library(ggbiplot)
library(ggfortify)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(plotly)
library(plyr)
library(scales)
library(grid)
library(devtools)
library(vctrs)

wd <- '/Users/gabriela/Desktop/Name'
setwd(wd)

readIn <- function (fileName, sampleName)
{
  dF <- read.csv(fileName, header = F, sep = '\t')
  dF <- as_tibble(dF)
  dF <- dF %>%
    select(-2)
  colnames(dF) <- c("Gene", sampleName)
  
  return (dF)
}

descGrabber <- function (fileName, sampleName)
{
  dF <- read.csv(fileName, header = F, sep = '\t')
  dF <- as_tibble(dF)
  dF <- dF %>%
    select(1:2)
  colnames(dF) <- c("Gene", "Description")
  
  return (dF)
  
}

descWhole <- descGrabber('01_name.txt', "Overexpression_A")

OEA <- readIn('01_name.txt', "Overexpression_A")
OEB <- readIn('02_name.txt', "Overexpression_B")
OEC <- readIn('03_name.txt', "Overexpression_C")
WTA <- readIn('04_name.txt', "WildType_A")
WTB <- readIn('05_name.txt', "WildType_B")
WTC <- readIn('06_name.txt', "WildType_C")

sapply(OEA, class)
OEA$Overexpression_A <- as.numeric(as.numeric(OEA$Overexpression_A))

OEASum <- aggregate(Overexpression_A~Gene,data=OEA,FUN=sum)
OEBSum <- aggregate(Overexpression_B~Gene,data=OEB,FUN=sum)
OECSum <- aggregate(Overexpression_C~Gene,data=OEC,FUN=sum)
WTASum <- aggregate(WildType_A~Gene,data=WTA,FUN=sum)
WTBSum <- aggregate(WildType_B~Gene,data=WTB,FUN=sum)
WTCSum <- aggregate(WildType_C~Gene,data=WTC,FUN=sum)

mergeByGene <- function(a,b,c,d,e,f)
{
  dF <- merge(a, b, by = "Gene")
  dF <- merge(dF, c, by = "Gene")
  dF <- merge(dF, d, by = "Gene")
  dF <- merge(dF, e, by = "Gene")
  dF <- merge(dF, f, by = "Gene")
  row.names(dF) <- dF$Gene
  
  return (dF)
}

AllData <- mergeByGene (OEASum, OEBSum, OECSum, WTASum, WTBSum, WTCSum)

make_ddsDataFrame <- function (names, data) 
{
  data <- data[6:nrow(data),]
  data <- as.matrix(data[1:nrow(data), 2:ncol(data)])
  colnames(data) <- colnames(data[2:ncol(data)])
  
  types <- c("OE", "OE", "OE", "WT", "WT", "WT") #edit this for your samples
  factor(types, levels = c("WT", "OE")) #force them to be listed in this order
  expData <- data.frame(SampleName = names,
                        type = types)
  ddsMat <- DESeqDataSetFromMatrix(countData = data,
                                   colData =expData,
                                   design = ~ type)
  ddsMat$type <- relevel(ddsMat$type, "WT")
  ddsMat <- ddsMat[ rowSums(counts(ddsMat)) > 1,]
  return (ddsMat)
}

Names <- c( "OverexpressionA", "OverexpressionB", "OverexpressionC", "WildTypeA", "WildTypeB", "WildTypeC")

dds <- DESeq(ddsMat_wholeGenome)
res <- results(dds)
summary(res)
res <- res[order(res$padj),]
table(res$padj <.01)


resSig <- subset(res, padj < 0.05)
resSig <- resSig[order(resSig$log2FoldChange),]

EnhancedVolcano (res, 
                 lab = NA, 
                 x = 'log2FoldChange', 
                 y = 'pvalue',
                 FCcutoff = 0.25,
                 pCutoff = 10e-5,
                 col = c("gray", "gray", "gray", "blue"))


vsd <- vst (ddsMat_wholeGenome, blind = FALSE)
rdl <- rlog (ddsMat_wholeGenome, blind = FALSE)
head (assay(vsd), 3)
plotPCA (vsd, intgroup=c ("type"))



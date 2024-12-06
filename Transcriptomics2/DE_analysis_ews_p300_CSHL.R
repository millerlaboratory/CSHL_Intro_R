library(tidyverse) 
library(DESeq2) #our DE analysis tool
library(pheatmap) # easy heatmap maker 
#library(ComplexHeatmap) better package to make more complicated heatmaps
library(EnhancedVolcano) #volcano plots
library(UpSetR) #better than a Venn Diagram
library(fgsea) #gsea via R

setwd("~/Google Drive/My Drive/Other/CSHL Comp Genomics/CSHL 2024/RNA-seq tutorial/")
geneNames <- read_tsv("ensGene.txt")
pcGenes <- filter(geneNames, `Gene type` == "protein_coding")

counts <- read_tsv("subread_counts.txt", comment = "#")
#counts <- mutate(counts, Geneid = str_match(Geneid,"ENSG\\w+")[,1])
pcCounts <- filter(counts, Geneid %in% pcGenes$`Gene stable ID version`)

meta <- readxl::read_xlsx("p300w_sampleSheet.xlsx")

##take out CIC-dux experiments
meta <- filter( meta, Cancer != "CIC-DUX4")

##make the tables you need for DESeq2
#combing all EWS and other cancer samples
pcCounts <- pcCounts[,c("Geneid",meta$SampleID)]
pcCounts <- column_to_rownames(pcCounts, "Geneid")
meta <- column_to_rownames(meta, "SampleID")
ddsAll <- DESeqDataSetFromMatrix(countData = pcCounts, 
                                 colData = meta,
                                 design = ~treatment + CellLine) ##why did I include cell line? 
ddsAll <- DESeq(ddsAll)

#create normalized gene level expression values
vstAll <- as.data.frame(assay(vst(ddsAll)))
vstAll$Geneid <- rownames(vstAll)
vstAll$var <- apply(vstAll, 1, function(x) var(x[1:34]))
vstAll <- arrange(vstAll, desc(var))

pheatmap(vstAll[1:5000,1:34],
         scale = "row",
         show_rownames = F,
         show_colnames = T,
         annotation_col = meta)


##Remove cell line effects with limma
vstAllCell <- limma::removeBatchEffect(vstAll[,1:34], batch = meta$CellLine)
vstAllCelldf <- as.data.frame(vstAllCell)
vstAllCelldf$geneid <- rownames(vstAllCell)
write_tsv(vstAllCelldf, file = "VSTnorm_cellLineRemoved_allSamples.txt" )

annotColors <- list(Cancer = c(EWS = "#01665e", Other = "#c7eae5"),
                    treatment = c(p300 = "#b2182b", Control = "#999999"),
                    CellLine = c(A4573 = "#a6cee3", A673 = "#1f78b4", 
                                 SKES = "#b2df8a", TC71 = "#33a02c"))
pheatmap(vstAllCell[1:5000,1:34],
         scale = "row",
         show_rownames = F,
         show_colnames = T,
         annotation_col = meta[,c("Cancer","treatment")],
         annotation_colors = annotColors)

##pca with corrected gene expression 
pca <- prcomp(t(vstAllCell))
pcaDF <- as.data.frame(pca$x)
pcaDF <- cbind(pcaDF, meta)

ggplot(pcaDF, aes(x = PC1, y = PC2, color = treatment, shape = Cancer)) +
  geom_point(size = 5) +
  scale_color_manual(values = annotColors$treatment) + 
  theme_bw() 


ggsave(filename = "PCA_allSamples.pdf", height = 5, width = 7)  

##DE genes in all EWS lines
metaEWS <- filter(meta, Cancer == "EWS")
countEWS <- pcCounts[,rownames(metaEWS)]
ddsEWS <- DESeqDataSetFromMatrix(countData = countEWS,
                                 colData = metaEWS,
                                 design = ~treatment + CellLine)
ddsEWS <- DESeq(ddsEWS)
resultsEWS <- results(ddsEWS, contrast = c("treatment", "p300", "Control"), tidy = T) %>%
  arrange(padj)
#resultsEWS <- mutate(resultsEWS, Geneid = str_extract(row, "ENSG\\d+") )
resultsEWS <- left_join(resultsEWS, geneNames, by = c("row" = "Gene stable ID version"))

write_tsv(resultsEWS, file = "DEresults_EWScellLines_p300-Ctl.txt")
EnhancedVolcano(resultsEWS, 
                x = "log2FoldChange",
                 y = "padj",
                lab = resultsEWS$`Gene name`,
                FCcutoff = 2,
                title = "Ewings Cancer Cell Lines",
                subtitle = "p300w - Control, log2FC >2, padj < 0.05"
                )
ggsave(filename = "volcano_EWScellLines.pdf", height = 9, width = 7)


##DE genes all Other Cell lines
metaOther <- filter(meta, Cancer == "Other")
countOther <- pcCounts[,rownames(metaOther)]
ddsOther <- DESeqDataSetFromMatrix(countData = countOther,
                                   colData = metaOther,
                                   design = ~treatment + CellLine)
ddsOther <- DESeq(ddsOther)
resultsOther <- results(ddsOther, 
                        contrast = c("treatment", "p300", "Control"), 
                        tidy = T) %>% arrange(padj)


#resultsOther <- mutate(resultsOther, Geneid = str_extract(row, "ENSG\\d+") )
resultsOther <- left_join(resultsOther, geneNames, by = c("row" = "Gene stable ID version"))
write_tsv(resultsOther, file = "DEresults_otherCellLines_p300-Ctl.txt")
EnhancedVolcano(resultsOther, 
                x = "log2FoldChange",
                y = "padj",
                FCcutoff = 2,
                lab = resultsEWS$`Gene name`,
                title = "Other Cancer Cell Lines",
                subtitle = "p300w - Control,log2FC >2, padj < 0.05"
)
ggsave("volcano_OthercellLines.pdf", height = 9, width = 7)

#####upset plot for DE results in EWS and Other Cancers 
#log2FC > 2
ewsDE <- filter(resultsEWS, padj < 0.05 & abs(log2FoldChange) >= 2)$"Gene name" %>% unique()
otherDE <- filter(resultsOther, padj < 0.05 & abs(log2FoldChange) >= 2)$"Gene name" %>% unique()
ewsUp <- filter(resultsEWS, padj < 0.05 & log2FoldChange >= 2)$row
ewsDown <- filter(resultsEWS, padj < 0.05 & log2FoldChange <= -2)$row
otherUp <- filter(resultsOther, padj < 0.05 & log2FoldChange >= 2)$row
otherDown <- filter(resultsOther, padj < 0.05 & log2FoldChange <= -2)$row

upDownlist <- list(EWS_Up = ewsUp,
                   EWS_Down = ewsDown,
                   Other_Up = otherUp,
                   Other_Down = otherDown)

pdf(file = "UpSet_EWSvsOther_logFC2cutoff.pdf", width = 7, height = 5 )
upset(fromList(upDownlist), order.by = "freq", point.size = 3.5, line.size = 2, 
      mainbar.y.label = "Gene Overlap", sets.x.label = "Number of DE Genes",
      text.scale = c(2, 1, 1.3, 1, 2, 2), ) 
dev.off()

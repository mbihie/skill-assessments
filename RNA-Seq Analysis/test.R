
#BiocManager::install("DESeq2")
library(DESeq2)
#BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
#BiocManager::install("pheatmap")
library(pheatmap)
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
#BiocManager::install("KEGGREST")
library(KEGGREST)
#BiocManager::install("pathview")
library(pathview)
#BiocManager::install("gage")
library(gage)
#BiocManager::install("gageData")
library(gageData)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("recount")
library(recount)
#install.packages("RColorBrewer")
library(RColorBrewer)
#install.packages("DT")
library(DT)
#install.packages("enrichR")
library(enrichR)
#install.packages("kableExtra")
library(kableExtra)




#read the data into R
#   class(rse): RangedSummarizedExperiment
rse <- readRDS("EwS.rds")

#remove the version number from the geneID
rownames(rse) <- gsub(rownames(rse),
                      pattern = "\\..+", replacement = "")

#make the dds
dds <- DESeqDataSet(rse, design = ~condition)

#set the factor level
dds$condition <- relevel(dds$condition, ref = "shCTR")

#creat DESeq 2 object
dds <- DESeq(dds)

#decrease the size of the DESeq object to make DESeq2 functions faster
dds <- estimateSizeFactors(dds)
idx <- rowSums( counts(dds, normalized=TRUE) >= 10 ) >= 3
dds <- dds[idx,]

#remove rows with low gene counts
keep <- as.data.frame(counts(dds)) %>%
  rowSums(counts(dds, normalized=TRUE)) >= 10
dds <- dds[keep,]


#create rlog for PCA
rld <- rlog(dds)

#plot the PCA
plotadj <- plotPCA(rld)

plotadj 

plotadj + coord_fixed(ratio = 1,
                      xlim=c(-50, 50),
                      ylim=c(-10, 10))


#extract the results
res <- results(dds,
               contrast = c("condition", 
                            "shCTR",      
                            "shEF1"),  
               alpha = 0.05)

#shrink data
resNorm <- lfcShrink(dds = dds,
                     res = res,
                     type = "normal",
                     coef = 2)

#MA Plot showing the relationship between mean count and log2 fold change.
plotMA(resNorm,
       main = "MA-plot of Normalized Ewing Sarcoma RNA-seq Data")

#make a dataframe of the results to view them
resdf <- as.data.frame(resNorm)

#extract gene symbol from EnsDb.Hsapiens.v86
ens2sym <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                 keys = keys(EnsDb.Hsapiens.v86), 
                                 columns = c("SYMBOL"))

#join ens2sym to resdf by shared column (GENEID)
resdfsym <- resdf %>%
  rownames_to_column() %>%
  mutate(GENEID = gsub(rowname, pattern = "\\..+", replacement = "")) %>% 
  inner_join(y = ens2sym, by = "GENEID") %>% 
  dplyr::select(-rowname) %>%
  mutate(padj = case_when(padj == 0 ~ .Machine$double.xmin,
                          TRUE ~ padj)) #replacing 0s with minimum R value

#DF to specify level of signficance
DEgenes <- resdfsym %>%
  arrange(padj) %>%
  dplyr::filter(padj < 0.05 & log2FoldChange > 2) %>%  
  #slice_head(n=10) %>%
  select("SYMBOL", "GENEID", "padj", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue") 

#data table of DE genes
datatable(DEgenes,
          class = 'cell-border stripe',
          caption = htmltools::tags$caption(
            style = 'caption-side: bottom; text-align: center;',
            'Table 1: ',
            width = 3,
            htmltools::em('Differentially Expressed Genes.')),
          filter = 'top',
          extensions = 'FixedColumns',
          options = list(
            pageLength = 5,
            autoWidth = TRUE
            ,
            scrollX = TRUE
          )
)

EnhancedVolcano(resdfsym, lab = resdfsym$SYMBOL, pCutoff = 0.05,
                FCcutoff = 2,
                x = "log2FoldChange",
                y = "padj",
                #xlim = c(-30, 30),
                title = "EWSR1-FLI1 Suppressed Ewing Sarcoma DEGs")

#filter df to significant DEGs that are both over expressed and underexpressed
orderedSig <-resdfsym %>%
  arrange(padj) %>%
  dplyr::filter(padj < 0.05 & log2FoldChange > 2 | log2FoldChange < -2) 

#create matrix for heat map
id1 <- orderedSig$GENEID 
id2 <- orderedSig$SYMBOL 
mat <- assay(rld) #extract mtx of nrmlzd counts
DE <- mat[id1,]
rownames(DE) <- id2 #rmv ensemble iDs from row name

#transform matrix to df to filter to top 10 DEGs
#overexpressed
overDE <- as.data.frame(DE) %>% #filter to top 10 overexpressed genes
  slice_head(n=10)
overDE <- as.matrix(overDE) #make it a matrix for pheatmap
#underexpressed
underDE <- as.data.frame(DE) %>% #filter to top 10 underrexpressed genes
  slice_tail(n=10)
underDE <- as.matrix(underDE) #make it a matrix for pheatmap


DE <- rbind(overDE, underDE) #bind both matrices into 1
annotation <- as.data.frame(colData(rld)[, c("run","condition")]) #add annotation onto map

#heatmap
pheatmap(DE,
         scale = "row", 
         clustering_distance_rows = "correlation", 
         annotation_col = annotation,
         show_colnames = FALSE, # removes coloumn names from run at the btm
         annotation_names_col = FALSE, #removes names beside the colours at the top
         main="Top 10 Differentially Expressed genes"
)

library(grid)
grid.ls(grid.force())

#annotation_legend.3-6-5-6

grid.gedit("GRID.gTree.297",gp = gpar(col="white")) #0:
grid.gedit("GRID.text.1161", gp = gpar(col="white")) #1: removes condition
grid.gedit("GRID.rect.1162", gp = gpar(col="white")) #2: removes box lines for schCTR/shEF1
grid.gedit("GRID.text.1163", gp = gpar(col="white")) #3: removes schCTR/shEF1
grid.gedit("GRID.text.1164", gp = gpar(col="white")) #4: removes run
grid.gedit("GRID.rect.1165", gp = gpar(col="white")) #5: removes box lines for SRR...
grid.gedit("GRID.text.1166", gp = gpar(col="white")) #6:removes SRR...

#legend.3-5-5-5
grid.gedit("GRID.gTree.1194",gp = gpar(col="white")) # : can't tell
grid.gedit("GRID.rect.1192", gp = gpar(col="white")) #7: removes box lines for grdnt
grid.gedit("GRID.text.1193", gp = gpar(col="white")) #8: remmoves grdnt #s




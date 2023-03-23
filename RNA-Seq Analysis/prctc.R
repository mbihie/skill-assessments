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
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
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


#reorder these results by p-value and call summary() on the results object
res = res[order(res$pvalue),]
summary(res)

res$symbol = mapIds(org.Hs.eg.db,
                    keys=row.names(res),
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")

res$entrez = mapIds(org.Hs.eg.db,
                    keys=row.names(res),
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")

res$name =   mapIds(org.Hs.eg.db,
                    keys=row.names(res),
                    column="GENENAME",
                    keytype="ENSEMBL",
                    multiVals="first")

#setup the KEGG data-sets we need.
data(kegg.sets.hs)
data(sigmet.idx.hs)

#remove unnecessary pathway defintions
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]

#creating vector of FCs with Entrez IDs as the names for the gage() function
foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez

#Get the results for the pathway analysis
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)

#join the over-expressed and under-expressed pathways
keggtable <- c(rownames(keggres$greater),rev(rownames(keggres$less)))

#edit the dataframe
keggtable <- as.data.frame(keggtable) 
keggtable <- keggtable %>%
  separate(keggtable,
           sep = 9,
           into=c("pathid","pathway")) %>% 
  select(pathway,pathid) 

#erm <- distinct(keggtable$pathid)

#format p value for downregulated pathways
ls <- as.data.frame(keggres$less) %>%
  rownames_to_column() %>%
  rename(pathway = rowname) %>%
  mutate(pathway = substring(pathway, 10)) %>% map_df(rev)

#format p value for upregulated pathways
grtr <- as.data.frame(keggres$greater) %>%
  rownames_to_column() %>%
  rename(pathway = rowname) %>%
  mutate(pathway = substring(pathway, 10))

#join dataframes
grtrls <- rbind(grtr,ls)
KT <- keggtable %>%
  full_join(grtrls, by="pathway") %>%
  select(pathway, pathid, p.val)

#present the results in a data table
datatable(KT, 
          options=list(pageLength=5),
          caption = htmltools::tags$caption(
            style = 'caption-side: bottom; text-align: center;',
            'Table 2: ',
            width = 3,
            htmltools::em('Top Over-Expressed & Under-Expressed KEGG Pathways.')
          )
)

# we want the log2 fold change 
original_gene_list <- resdf$log2FoldChange

# name the vector
names(original_gene_list) <- rownames(resdf)

# omit any NA values 
original_gene_list <- na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
original_gene_list = sort(original_gene_list, decreasing = TRUE)

# Convert gene IDs for gseKEGG function
#lose some genes here because not all IDs will be converted
ids <- bitr(geneID = names(original_gene_list), 
            fromType = "ENSEMBL", 
            toType = "ENTREZID", 
            OrgDb = org.Hs.eg.db)

# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

df2 = resdf[rownames(resdf) %in% dedup_ids$ENSEMBL,]

df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene universe
kegg_gene_list <- df2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "hsa"
set.seed(1234)
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

out <- kk2@result

datatable(out)

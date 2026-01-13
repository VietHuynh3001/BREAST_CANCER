# SET WORK DIRECTORY
setwd('D:/MSc/demo_TCGA_BRCA/')
# IMPORT LIBRARIES
library(DESeq2)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(limma)



# READ FILE: COUNT DATA, NORMALIZED DATA, META DATA
data_count <- read.csv('D:/MSc/demo_TCGA_BRCA/Final data/count_data.csv', row.names = 1)
data_normalized <- read.csv('D:/MSc/demo_TCGA_BRCA/Final data/normalized_data.csv', row.names = 1)
data_meta <- readRDS('D:/MSc/demo_TCGA_BRCA/Final data/meta_data.rds')
# PREPROCESS CLINICAL DATA
row.names(data_meta)<-gsub('-',".",row.names(data_meta))
# MAKE SURE THAT THE ROW NAMES IN RNA_SEQ DATA MATCH TO THE COLUMN NAMES IN META DATA
all(colnames(data_normalized) %in% rownames(data_meta))
all(colnames(data_count) %in% rownames(data_meta))
# ARE THEY IN THE SAME ORDER
all(colnames(data_normalized) == rownames(data_meta))
all(colnames(data_count) == rownames(data_meta))
# GRAPHIC PRESENTATION: PCA





# CONSTRUCT DESEQDATASET OBJECT
dds <- DESeqDataSetFromMatrix(countData = data_count,
                              colData = data_meta,
                              design = ~ follow_ups_disease_response)
# FILTER GENES >10
keep<-rowSums(counts(dds))>=10
dds<-dds[keep,]
# RE-LEVEL
dds$follow_ups_disease_response <-relevel(x = dds$follow_ups_disease_response,ref = 'TF-Tumor Free')
# RUN DESEQ
dds <- DESeq(dds)
res <- results(dds)
res <-as.data.frame(res)
# SAVE RESULTS
write.csv(res,'D:/MSc/demo_TCGA_BRCA/results/DGE_analysis.csv')
saveRDS(dds,'D:/MSc/demo_TCGA_BRCA/results/dds.rds')




#READ DGE ANALYSIS RESULTS
dds_read <-readRDS('D:/MSc/demo_TCGA_BRCA/results/dss.rds')
df_result<-results(dds_read)


df_result<-read.csv('D:/MSc/demo_TCGA_BRCA/results/DGE_analysis.csv',row.names = 1)
df_result<-na.omit(df_result) # DROP ROWS WHICH CONTAIN NULL VALUE
# FIND TOP UP-REGULATED GENES:
top.up.regulated.gene <- row.names(df_result[(df_result$padj<0.05)&(df_result$log2FoldChange>0.5),])
paste('The number of up-regulated genes:',length(top.up.regulated.gene))
# FIND TOP DOWN-REGULATED GENES:
top.down.regulated.gene <-row.names(df_result[(df_result$padj<0.05)&(df_result$log2FoldChange<(-0.5)),])
paste('The number of down regulated genes:',length(top.down.regulated.gene))
# FIND TOP OVER-EXPRESSED GENES:
selected_df <- df_result[(df_result$padj <0.05)&(abs(df_result$log2FoldChange)>0.5),]
top.over.expressed.genes <-row.names(selected_df)
paste('The number of over-expressed genes:',length(top.over.expressed.genes))

# GRAPHIC PRESENTATION: VOCANO PLOT 
old.pal <-palette(c('#FAAC68','#6DC3BB')) # SET A PALATTE 

mar <- par(mar=c(4,4,2,1),cex.main=1.5) # SET A MARGIN

plot(df_result$log2FoldChange,
     -log10(df_result$padj),
     xlab = 'log2FoldChange',
     ylab = '-log(padj)',
     main = 'UPREGULATED GENES VS DOWNREGULATED GENES',
     cex=0.5,
     pch=20)
with(subset(df_result,abs(log2FoldChange)>=0.5 & padj<0.05),
     points(log2FoldChange,
            -log10(padj),
            pch=20,
            col=(sign(log2FoldChange)+3)/2,
            cex=0.75))
legend('bottomleft', 
       legend = c('up','down'), 
       title = paste('Padj<',0.05,sep = ''),
       pch=20,
       col = 1:2,)

# GRAPHIC PRESENTATION: MA PLOT

# GRAPHIC PRESENTATION: CLUSTERING HEATMAP





# GENE ENRICHMENT ANALYSIS: UP-REGULATED GENES
enrich.analysis.up.genes <- enrichGO(top.up.regulated.gene,
                                     OrgDb = org.Hs.eg.db,
                                     ont = 'ALL',
                                     keyType = 'ENSEMBL',
                                     pvalueCutoff = 0.05,
                                     qvalueCutoff = 0.05,
                                     readable = TRUE)
dotplot(enrich.analysis.up.genes,
        x='GeneRatio',
        color='p.adjust',
        title='Top 15 of GO enrichment for up-regulated genes',
        showCategory = 15,
        label_format =100,
        split='ONTOLOGY') + facet_grid(ONTOLOGY~.,scales='free')
# GENE ENRICHMENT ANALYSIS: DOWN-REGULATED GENES
enrich.analysis.down.genes <- enrichGO(top.down.regulated.gene,
                                       OrgDb = org.Hs.eg.db,
                                       keyType = 'ENSEMBL',
                                       ont = 'ALL',
                                       pvalueCutoff = 0.05,
                                       qvalueCutoff = 0.05,
                                       readable = TRUE)
dotplot(enrich.analysis.down.genes,
        x='GeneRatio',
        title='Top 15 of GO enrichment for down-regulated genes',
        color='p.adjust',
        showCategory=15,
        label_format=100,
        split='ONTOLOGY') +facet_grid(ONTOLOGY~.,scale='free')
enrich.analysis.genes<-enrichGO(gene = top.over.expressed.genes,
                                OrgDb = org.Hs.eg.db,
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05,
                                ont = 'ALL',
                                keyType = 'ENSEMBL',
                                readable = TRUE)
dotplot(enrich.analysis.genes,
        x='GeneRatio',
        title='Top 15 of GO enrichment for genes',
        color='p.adjust',
        showCategory=15,
        label_format=100,
        split='ONTOLOGY')+facet_grid(ONTOLOGY~., scales = 'free')
# SAVE RESULTS OF GO ENRICHMENT AS RDS FILE
saveRDS(enrich.analysis.up.genes,
        file = 'D:/MSc/demo_TCGA_BRCA/results/GO_enrichment_for_up_genes.rds')
saveRDS(enrich.analysis.down.genes,
        file = 'D:/MSc/demo_TCGA_BRCA/results/GO_enrichment_for_down_genes.rds')
saveRDS(enrich.analysis.genes,
        file = 'D:/MSc/demo_TCGA_BRCA/results/GO_enrichment_for_genes.rds')
# LIST BIOLOGICALLY-SIGNIFICANT GENES
bio_sig_genes<-function(enrich.result){
  temp<-as.data.frame(enrich.result)
  ontology=unique(temp$ONTOLOGY)
  list_genes<-c()
  for (i in ontology){
    genes=temp[temp$ONTOLOGY==i,]$geneID
    genes=unique(unlist(strsplit(genes,'/')))
    list_genes=append(list_genes,genes)
  }
  return(unique(list_genes))
}




# FILTER DEG IN DATAFRAME
filter_deg <-function(data.frame,list.deg) {
  data_normalized_deg<-data.frame[row.names(data.frame) %in% list.deg,]
  gene_id <- row.names(data.frame)
  gene_symbol <- mapIds(org.Hs.eg.db,
                        keys = gene_id,
                        column = 'SYMBOL',
                        keytype = 'ENSEMBL',
                        multiVals = 'first')
  gene_symbol <- ifelse(is.na(gene_symbol),gene_id,gene_symbol)
  row.names(data_normalized_deg)<-gene_symbol
  return(data_normalized_deg)
}
# FILTER UP-REGULATED GENES IN NORMALIZED DATA
data_normalized_up_gene<-filter_deg(data_normalized,top.up.regulated.gene)
# FILTER DOWN-REGULATED GENES IN NORMALIZED DATA
data_normalized_down_gene<-filter_deg(data_normalized,top.down.regulated.gene)
# FILTER DIFFENTAILLY-EXPRESSED GENES
data_normalized_DEG <- filter_deg(data_normalized,top.over.expressed.genes)



# SAVE UP-REGULATED GENES, DOWN-REGUALATED, DE GENES:
write.csv(x = data_normalized_up_gene,'D:/MSc/demo_TCGA_BRCA/results/normalized_data_up_genes.csv')
write.csv(x = data_normalized_down_gene,'D:/MSc/demo_TCGA_BRCA/results/normalized_data_down_genes.csv')
write.csv(x = data_normalized_DEG,'D:/MSc/demo_TCGA_BRCA/results/normalized_data_DEG.csv')


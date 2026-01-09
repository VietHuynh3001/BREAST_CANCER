# SET WORKING DIRECTORY
setwd('D:/MSc/demo_TCGA_BRCA/')
getwd()
# LIBRARY
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyr)



# GET PROJECTS
GDC_projects<-getGDCprojects()
GDC_projects
# GET TGCA_BRCA SUMMARY
summary_TGCA_BRCA <-getProjectSummary(project = 'TCGA-BRCA')
summary_TGCA_BRCA
# QUERY TCGA_BRCA
query_TCGA_BRCA <-GDCquery(project = 'TCGA-BRCA',
                           data.category = 'Transcriptome Profiling',
                           access = 'open',
                           experimental.strategy = 'RNA-Seq',
                           workflow.type = 'STAR - Counts')
output_qury_TCGA_BRCA <- getResults(query = query_TCGA_BRCA)
# SELECT 1000 FIRST SAMPLES
set.seed(100)
selected_output <-output_qury_TCGA_BRCA[sample( x = nrow(output_qury_TCGA_BRCA),size = 1000),]
query_TCGA_BRCA$results[[1]] <-selected_output
# DOWLOAD SELECTED PATIENTS
# GDCdownload(query = query_TCGA_BRCA)
# PREPARE DATA
tcga_brca_data <- GDCprepare(query = query_TCGA_BRCA,
                             summarizedExperiment = TRUE,
                             directory = 'D:/MSc/breast cancer data/GDCdata')
tcga_brca_data




# GET RNA SEQ DATA (COUNT DATA)
BRCA_matrix <- assay(tcga_brca_data,'unstranded')
BRCA_matrix
# GET RNA SEQ DATA (NORMALIZED DATA)
# unstranded, stranded_first, stranded_second, tpm_unstrand, fpkm_unstrand, fpkm_uq_unstrand
BRCA_matrix_normalized <- assay(tcga_brca_data,'fpkm_unstrand')
BRCA_matrix_normalized
# GET CLINICAL DATA
clinical_from_data <-colData(tcga_brca_data)
clinical_data <- as.data.frame(clinical_from_data)




# SAVE FILE
write.csv(BRCA_matrix,'D:/MSc/demo_TCGA_BRCA/Final data/count_data.csv')
write.csv(BRCA_matrix_normalized,'D:/MSc/demo_TCGA_BRCA/Final data/normalized_data.csv')
saveRDS(clinical_data,'D:/MSc/demo_TCGA_BRCA/Final data/meta_data.rds')




# READ DATA
BRCA_matrix<-read.csv('D:/MSc/demo_TCGA_BRCA/Final data/count_data.csv',row.names = 1)
BRCA_matrix_normalized<-read.csv('D:/MSc/demo_TCGA_BRCA/Final data/normalized_data.csv',row.names = 1)
clinical_data<-readRDS('D:/MSc/demo_TCGA_BRCA/Final data/meta_data.rds')
#PREPROCESS CLINICAL DATA
clinical_data$follow_ups_disease_response[!duplicated(clinical_data$follow_ups_disease_response)] # CHECK VALUES
clinical_data<-clinical_data[clinical_data$follow_ups_disease_response %in% c("TF-Tumor Free","WT-With Tumor"),]
patient.id <- row.names(clinical_data)

# PREPROCESS READ DATA
gene_id <- gsub("\\..*","",row.names(data_normalized))
# DUPLICATED GENES
non_dup_gene <- !duplicated(gene_id)
# ELIMINATE DUPLICATED GENES
BRCA_matrix<-BRCA_matrix[non_dup_gene,]
BRCA_matrix_normalized<-BRCA_matrix_normalized[non_dup_gene,]
# RENAME ROWS
row.names(BRCA_matrix)<-gene_id[non_dup_gene]
row.names(BRCA_matrix_normalized)<-gene_id[non_dup_gene]
# FILTER PATIENT
BRCA_matrix<-BRCA_matrix[,colnames(BRCA_matrix) %in% patient.id]
colnames(BRCA_matrix_normalized)<-gsub("\\.","-",colnames(BRCA_matrix_normalized))
BRCA_matrix_normalized<-BRCA_matrix_normalized[,colnames(BRCA_matrix_normalized) %in% patient.id]




# RE-CHECK DATA
dim(BRCA_matrix)
dim(BRCA_matrix_normalized)
dim(clinical_data)



# SAVE FILE
write.csv(BRCA_matrix,'D:/MSc/demo_TCGA_BRCA/Final data/count_data.csv')
write.csv(BRCA_matrix_normalized,'D:/MSc/demo_TCGA_BRCA/Final data/normalized_data.csv')
saveRDS(clinical_data,'D:/MSc/demo_TCGA_BRCA/Final data/meta_data.rds')

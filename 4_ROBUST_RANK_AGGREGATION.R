# IMPORT PACKAGES
setwd('D:/MSc/demo_TCGA_BRCA')
library(RobustRankAggreg)
library(pheatmap)

# READ DATA
df=read.csv('D:/MSc/demo_TCGA_BRCA/results/feature_importance.csv',row.names = 1)
total_gene<-nrow(df)
df['TM4SF4',]

# SORT GENE BY DECREASING ORDER OF FEATURE IMPORTANCE
list_model_1<-row.names(df[order(df$XGBoost,decreasing = TRUE),])
list_model_2<-row.names(df[order(df$Random.forest,decreasing = TRUE),])
list_model_3<-row.names(df[order(df$LightGBM,decreasing = TRUE),])
list_model_4<-row.names(df[order(df$Cat.Boost,decreasing = TRUE),])

# Robust Rank Aggregation
input_list<-list(list_model_1,list_model_2,list_model_3,list_model_4)
mat <-rankMatrix(glist = input_list,N=total_gene)
result<-aggregateRanks(glist = input_list,N = total_gene)

# Graphical presentation
top_10_important_genes<-result$Name[1:20]
filtered_mat<-mat[top_10_important_genes,]
pheatmap(filtered_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         main = 'Robust Rank Aggregation heatmap')
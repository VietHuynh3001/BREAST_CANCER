# 1.Introduction
This personal project is conducted to identify potential gene biomarkers for Breast Cancer (BRCA) using the TCGA-BRCA dataset. The main objective is to develop a computational pipeline where multiple Machine Learning (ML) models are combined to find the most significant genes for distinguishing between cancer and normal samples.

In this project, an ensemble-based feature selection strategy is implemented. Four powerful models, including Random Forest, XGBoost, CatBoost, and LightGBM, are utilized to rank genes based on their feature importance. Since features are ranked differently by each model, the Robust Rank Aggregation (RRA) method is applied to integrate these individual lists into a single "consensus" ranking. By using this approach, the bias of any single algorithm is minimized, and the selected biomarkers are ensured to be consistently significant across all models.

# 2.Method
## 1.1. Data Pre-processing and Differential expression analysis
The TCGA-BRCA dataset is utilized as the primary data source. After initial cleaning and normalization, Differential Expression Analysis is performed to compare cancer and normal samples. To select potential features, a threshold of $|\log_2FC| > 0.5$ and $P < 0.05$ is applied. Only genes meeting these criteria are kept for the next steps.
## 1.2. Embedded Feature Selection
Four ensemble-based Machine Learning models are employed to identify the most predictive genes: Random Forest, XGBoost, CatBoost, and LightGBM. These models are trained on the pre-selected genes to perform classification. Embedded Feature Importance is then calculated by each algorithm, and genes are ranked based on their contribution to the model's performance.
## 1.3. Robust Rank Aggregation
The individual ranking lists from the four ML models are integrated using the Robust Rank Aggregation (RRA) method. This statistical step is applied to identify a "consensus" list of genes. A final $P$-value is assigned to each gene, and only those that are consistently ranked at the top across all algorithms are selected as candidate biomarkers.
## 1.4. Oncology Enrichment Analysis
To explore the biological significance of the identified biomarkers, Enrichment Analysis is conducted. Specifically, Gene Ontology (GO) terms are analyzed to determine the biological processes and molecular functions associated with these genes. This step is performed to ensure the selected biomarkers are relevant to cancer biology.

# 3.Results:
## 3.1. Differential Expression Analysis (DEA)

In the first stage of this project, Differential Expression Analysis (DEA) was conducted on the TCGA-BRCA dataset to identify genes with significant expression changes between cancer and normal samples. Based on the statistical threshold of |log₂FC| > 0.5 and P < 0.05, a total of 1,340 differentially expressed genes (DEGs) were identified. Within this set, 456 genes were found to be up-regulated, while 886 genes were classified as down-regulated. The distribution and significance of these candidates are clearly illustrated in the Volcano plot, where specific focus was placed on the up-regulated genes due to their potential roles in diagnosis and treatment.

<p align="center">
  <img src="Vocano_plot.png" width="600"/>
  <br>
  <i>
    Figure 1: Volcano plot representing 1,340 differentially expressed genes, including 456 up-regulated and 886 down-regulated candidates.
  </i>
</p>

Following the identification of these DEGs, the 456 up-regulated genes are selected as the primary input features for the next stage, which is Machine Learning training. These features will be used to train multiple models to evaluate and rank their importance.

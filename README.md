# Melanoma Analysis Pipeline

This R pipeline integrates microarray and RNASeq data to identify prognostic markers in melanoma. The key steps are:

1. **Differentially Expressed CAFs**  
   Filter TCGA and GEO DEGs (FDR < 0.05, |logFC| ≥ 1) and identify common DEGs with CAF markers.

2. **Univariate Cox Regression**  
   Subset expression data with common CAF DEGs and identify significant prognostic genes (p < 0.05).

3. **LASSO Regression**  
   Perform LASSO regression (with cross-validation) to select robust predictors.

4. **Multivariate Cox Regression**  
   Refine candidate genes via multivariate Cox analysis (p < 0.05).

5. **Risk Score Computation**  
   Compute a weighted risk score using gene coefficients, define risk groups (High/Low) based on the median, and merge with clinical data.

6. **GSVA Analysis of Pathways and Immune Cells**  
   Compute GSVA scores for critical pathways involved in cancer and immune cells participating in the tumor microenvironment (TME).

7. **Correlation Analysis**  
   Evaluate the association between identified markers, dysregulated pathways, and immune cells.

8. **Results Saving**  
   Save outputs as CSV files for univariate Cox results, LASSO-selected genes, multivariate Cox results, and clinical data with risk scores.

## File Structure

- **data/input**: Input files (clinical data, gene sets, raw counts, normalized expression matrix, etc.)
- **data/output**: Result CSV files
- **scripts**: Custom functions and main analysis script
- **README.md**: This file

## How to Run

1. **Install Required R Packages**

   You can install the required packages manually. First, install CRAN packages:

   ```r
   install.packages(c("dplyr", "data.table", "glmnet", "plyr", "corrr", "survival"))
   
   if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
   BiocManager::install(c("edgeR", "DESeq2", "MsigDb", "GAVA", "limma", "org.Hs.eg.db", "AnnotationDbi"))


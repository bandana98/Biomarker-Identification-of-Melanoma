setwd("D:/Band_research/Annie/Analysis/Analysis")
source("scripts/functions.R")

#1. Differentially Expressed CAFs

# Load CAF markers
CAFs <- read.csv("data/input/CAFs.csv")  

# Load and filter TCGA DEGs results
load("data/output/tcga_degs_TvsN.RData")  
tcga_degs_TvN_filtered <- subset(tcga_degs_TvN, FDR < 0.05 & abs(logFC) >= 1)

# Load and filter GEO DEGs results
load("data/output/gse15605_degs_TvN.RData")  
geo_degs_TvN_filtered <- subset(gse15605_degs_TvN, adj.P.Val < 0.05 & abs(logFC) >= 1)

# Identify common DEGs among TCGA, GEO, and CAF markers
de_CAFs <- intersect(
  intersect(rownames(tcga_degs_TvN_filtered), rownames(geo_degs_TvN_filtered)),
  CAFs$Gene)

# Report the results
cat("Number of common DEGs among TCGA, GEO, and CAF markers:", length(de_CAFs), "\n")

#2. Univariate Cox Regression

# Load clinical data and normalized expression data (tumor samples)
clinical_data <- read.csv("data/input/clinical_T.csv", row.names = 1)
load("data/input/tcga_exprs_T_norm.RData")  

# Convert overall survival status to binary (0 = Alive, 1 = Dead)
clinical_data <- clinical_data %>%
  mutate(os_status = ifelse(os_status == "Alive", 0, 1))

# Use the common CAF DEGs to subset the expression matrix
de_CAF_exprs <- tcga_exprs_T_norm[de_CAFs, ]

# Run univariate Cox regression 
univ_result <- univariate_cox(clinical_data, de_CAF_exprs)
univ_result <- subset(univ_result, P_Value < 0.05)

# Report the results
cat("Number of significant genes in univariate Cox regression:", nrow(univ_result), "\n")

#3. Perform lasso regression

# Use the common CAF DEGs to subset the expression matrix
univ_CAF_exprs <- tcga_exprs_T_norm[univ_result$Gene,]

# Apply the function
result <- run_lasso(clinical_data, univ_CAF_exprs, response_variable = "sample_type", family_type = "binomial")

#Extract the LASSO model and the selected genes with coefficients
lasso_result <- result$coefficients

# Report the results
cat("Number of significant genes in Lasso regression:", nrow(lasso_result), "\n")


#3. Run multivariate Cox regression 
multi_cox_result <- multivariate_cox(clinical_data, tcga_exprs_T_norm, lasso_result$Gene)
multi_cox_result <- subset(multi_cox_result, P_Value < 0.05)

# Report the results
cat("Number of significant genes in multivariate Cox regression:", nrow(multi_cox_result), "\n")

#4. Compute Risk score
# Load raw expression data or normalized data
load("data/input/tcga_exprs_T.RData")  
clinical_data <- read.csv("data/input/clinical_T.csv", row.names = 1)

# Subset the scaled matrix for genes from the multivariate Cox results and Compute the weighted risk scorefor each sample (using the provided coefficients)
multi_CAF_exprs <- as.data.frame(t(scale(tcga_exprs_T)[multi_cox_result$Gene, ]))
multi_CAF_exprs$risk_score <- as.vector(as.matrix(multi_CAF_exprs) %*% multi_cox_result$Coefficient)

# Define risk groups based on the median risk score
cutoff <- median(multi_CAF_exprs$risk_score)
multi_CAF_exprs$Risk_group <- ifelse(multi_CAF_exprs$risk_score >= cutoff, "High risk", "Low risk")

# Merge Risk Score with Clinical Data 
clinical_data <- cbind(clinical_data, multi_CAF_exprs[, c("risk_score", "Risk_group")])

# Rename column "X" to "barcode"
colnames(clinical_data)[colnames(clinical_data) == "X"] <- "barcode"

# Save the resultsq
write.csv(univ_result, "data/output/univariate_result.csv")   
write.csv(lasso_result, "data/output/lasso_regression_result.csv")
write.csv(multi_cox_result, "data/output/multivariate_result.csv")
write.csv(clinical_data, "data/output/clinical_data_with_RS.csv")

# End of the script























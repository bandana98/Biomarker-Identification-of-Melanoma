# ------------------------------------------------------------------------------
# GSVA Analysis Pipeline
# This script performs GSVA enrichment analysis on critical pathways and immune cell gene sets 
# using clinical data with risk scores and normalized expression data.
# ------------------------------------------------------------------------------

# Load required libraries (install if missing)
required_pkgs <- c("GSVA", "limma", "dplyr")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("Package", pkg, "is required but not installed."))
  }
}
library(GSVA)
library(dplyr)

# Read clinical data containing risk scores
clinical_data <- read.csv("data/output/clinical_data_with_RS.csv", row.names = 1)

# Read normalized expression matrix (loads object 'tcga_exprs_T_norm')
load("data/input/tcga_exprs_T_norm.RData")

#1. Process geneset
# Read the Pathways gene set file and convert it into a list of character vectors
pathway_geneset <- read.csv("data/input/pathway_geneset.csv") %>%
  group_by(gs_name) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  ungroup() %>%
  split(.$gs_name, .$gs_name) %>%
  lapply(function(df) df$gene_symbol)

# Read the immune cell gene set file and convert it into a list of character vectors
ICs_geneset <- read.csv("data/input/ICs_gene_set.csv") %>%
  group_by(cell_type) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  ungroup() %>%
  split(.$cell_type, .$cell_type) %>%  # Split by cell_type (each element is a data frame)
  lapply(function(df) df$gene_symbol)   # Convert each data frame to a vector of gene symbols


#2. GSVA Enrichment Analysis
# Compute GSVA enrichment scores for the Pathways gene sets
gsvapar_PA <- gsvaParam(exprData = tcga_exprs_T_norm, geneSets = pathway_geneset, maxDiff = TRUE)
gsva_score_PA <- gsva(gsvapar_PA)

# Compute GSVA enrichment scores for the immune cell gene sets
gsvapar_ICs <- gsvaParam(exprData = tcga_exprs_T_norm, geneSets = ICs_geneset, maxDiff = TRUE)
gsva_score_ICs <- gsva(gsvapar)

# Save the GSVA enrichment scores   
save(gsva_score_PA, file = "data/output/gsva_score_ICs.RData")
save(gsva_score_ICs, file = "data/output/gsva_score_ICs.RData")

# End of the script
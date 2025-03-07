#1.  Correlation analysis between GSVA scores of critical pathways and risk scores and markers

# Load required libraries
library(dplyr)
library(corrr)

# Load GSVA scores and clinical data with risk scores
load("data/output/gsva_score_ICs.RData")
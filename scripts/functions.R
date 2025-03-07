#1. Function to plot Venn Diagram
VennDiagram <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

#2. Function to perform univariate Cox regression analysis
univariate_cox <- function(clinical_df, exprs_df) {
  # Load required libraries 
  library(survival)
  library(dplyr)

  combined_df <- cbind(clinical_df, t(exprs_df))     # Combine survival data with expression matrix
  genes <- rownames(exprs_df)     # Extract gene names
  
  univ_results <- lapply(genes, function(gene) {
    model <- coxph(as.formula(paste("Surv(os_time, os_status) ~", gene)), data = combined_df)
    summary_model <- summary(model)
    conf_int <- confint.default(model)  # Extract 95% CI
    
    data.frame(
      Gene = gene,
      Coefficient = signif(summary_model$coefficients[gene, "coef"], 2),
      HR = signif(exp(summary_model$coefficients[gene, "coef"]), 2),
      Lower.95 = signif(exp(conf_int[gene, 1]), 2),
      Upper.95 = signif(exp(conf_int[gene, 2]), 2),
      Wald_Test = signif(summary_model$waldtest["test"], 2),
      P_Value = signif(summary_model$waldtest["pvalue"], 2),
      stringsAsFactors = FALSE  # Prevents automatic factor conversion
    )
  })
  
  # Convert list to data frame
  univ_results <- do.call(rbind, univ_results)
  rownames(univ_results) <- NULL  # Remove row names
  
  return(univ_results)
}

#3. Function for multivariate Cox regression analysis
multivariate_cox <- function(clinical_df, exprs_df, genes) {
  
  # Convert overall survival status to binary (0 = Alive, 1 = Dead)
  clinical_df <- clinical_df %>% mutate(os_status = ifelse(os_status == "Alive", 0, 1))

  formula <- as.formula(paste("Surv(os_time, os_status) ~", paste(genes, collapse = " + ")))
  combined_df <- cbind(clinical_df, t(exprs_df[genes, , drop = FALSE]))
  summary_cox <- summary(coxph(formula, data = combined_df))
  
  as.data.frame(cbind(summary_cox$coefficients, summary_cox$conf.int))[, c(1:2, 5, 8:9)] %>%
    setNames(c("Coefficient", "HR", "P_Value", "Lower.95", "Upper.95")) %>% 
    tibble::rownames_to_column("Gene") %>% dplyr::select(Gene, everything())
}

# 4. Function to perform Lasso regression
run_lasso <- function(clinical_data, 
                      norm_exprs, 
                      response_variable, 
                      alpha_value = 1,
                      nfolds = 10, 
                      family_type = "binomial") {
  
  # Load required libraries
  library(data.table)
  library(glmnet)
  library(plyr)
  
  # Create local copies to avoid modifying global objects
  clinical_data_local <- as.data.frame(clinical_data)
  norm_exprs_local   <- as.data.frame(norm_exprs)
  
  # Ensure clinical data has a 'barcode' column; if missing, add it using rownames
  if (!"barcode" %in% colnames(clinical_data_local)) {
    clinical_data_local$barcode <- rownames(clinical_data_local)
  }
  
  # Transpose the expression matrix and add a 'barcode' column
  norm_exprs_transposed <- as.data.frame(t(norm_exprs_local))
  norm_exprs_transposed$barcode <- rownames(norm_exprs_transposed)
  
  # Merge the transposed expression data with clinical data by 'barcode'
  DATA_Model <- plyr::join(norm_exprs_transposed, clinical_data_local, 
                           by = "barcode", type = "left", match = "all")
  
  # Remove clinical columns from the merged data to get the expression matrix for modeling
  DATA_Model_lasso <- DATA_Model[, !(colnames(DATA_Model) %in% colnames(clinical_data_local))]
  
  # Check that the response variable exists in the merged data
  if (!(response_variable %in% colnames(DATA_Model))) {
    stop(paste("Response variable", response_variable, "not found in merged data."))
  }
  
  # Extract and, if necessary, convert the response variable (for binomial, ensure 0/1 numeric)
  y <- DATA_Model[[response_variable]]
  if (family_type == "binomial" && !is.numeric(y)) {
    y_factor <- as.factor(y)
    if (length(levels(y_factor)) != 2) {
      stop("Response variable for binomial family must have exactly 2 levels.")
    }
    y <- as.numeric(y_factor) - 1
  }
  
  # Fit the LASSO model using cross-validation
  cv.lassoModel <- glmnet::cv.glmnet(
    x = data.matrix(DATA_Model_lasso),
    y = y,
    standardize = TRUE,
    alpha = alpha_value,
    nfolds = nfolds,
    family = family_type,
    parallel = FALSE  # ensure reproducibility by disabling parallel processing
  )
  
  # Extract coefficients at lambda.1se (more regularized solution)
  coefs <- coef(cv.lassoModel, s = cv.lassoModel$lambda.1se, exact = TRUE)
  # Identify nonzero predictors (exclude the intercept)
  nonzero_genes <- rownames(coefs)[coefs[, 1] != 0]
  nonzero_genes <- nonzero_genes[nonzero_genes != "(Intercept)"]
  
  coeff <- data.frame(Gene = nonzero_genes, Coefficients = coefs[nonzero_genes, 1])
  
  return(list(lassoModel = cv.lassoModel, coefficients = coeff))
}




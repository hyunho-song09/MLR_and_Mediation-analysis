# 01. Multiple Linear Regression (for multiple features)
# This script performs multiple linear regression (MLR) in parallel over combinations of variables
# from two datasets (X and Y), optionally including covariates (C).
# Progress bars and parallel processing are used to make the computation faster.

# Required Libraries
library(dplyr)
library(doParallel)
library(doSNOW)

# Function to Perform Multiple Linear Regression
meta.MLR.calculator <- function (X, Y, C = NULL, covariate = FALSE) {
  
  print("Running Multiple Linear Regression...")
  
  # Create all combinations of columns between X and Y
  combinations <- expand.grid(X = 1:ncol(X), Y = 1:ncol(Y))
  
  # Register parallel backend
  n_cores <- detectCores() - 6 # Reserve 6 cores for other processes
  cl <- makeCluster(n_cores)
  registerDoSNOW(cl)
  
  # Setup Progress Bar
  pb <- txtProgressBar(max = ncol(X), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  # Process for when covariates are included
  process_combination_cov <- function(comb) {
    set.seed(2023)
    x <- comb$X
    y <- comb$Y
    
    # Combine X, Y, and Covariates, then omit NA rows
    df <- na.omit(cbind(X[, x, drop = FALSE], Y[, y, drop = FALSE], C))
    
    # Create covariate formula
    cov_formula <- paste0("`", colnames(C), "`", collapse = "+")
    
    # Fit the linear model with covariates
    MLR.fit <- lm(paste0("`", colnames(df)[2], "` ~ `", colnames(df)[1], "` + ", cov_formula), data = df)
    
    # Extract MLR results
    MLR.results <- coef(summary(MLR.fit))
    MLR.results_df <- data.frame(t(MLR.results[2,]))
    colnames(MLR.results_df) <- c("coef", "std_ERROR", "tval", "Pval")
    
    # Add feature names
    MLR.results_df$X <- colnames(df)[1]
    MLR.results_df$Y <- colnames(df)[2]
    
    # Reorder columns for output
    MLR.results_df <- MLR.results_df[, c(5, 6, 1, 4, 2, 3)]
    return(MLR.results_df)
  }

  # Process for when covariates are NOT included
  process_combination_no_cov <- function(comb) {
    set.seed(2023)
    x <- comb$X
    y <- comb$Y
    
    # Combine X and Y, then omit NA rows
    df <- na.omit(cbind(X[, x, drop = FALSE], Y[, y, drop = FALSE]))
    
    # Fit the linear model without covariates
    MLR.fit <- lm(paste0("`", colnames(df)[2], "` ~ `", colnames(df)[1], "`"), data = df)
    
    # Extract MLR results
    MLR.results <- coef(summary(MLR.fit))
    MLR.results_df <- data.frame(t(MLR.results[2,]))
    colnames(MLR.results_df) <- c("coef", "std_ERROR", "tval", "Pval")
    
    # Add feature names
    MLR.results_df$X <- colnames(df)[1]
    MLR.results_df$Y <- colnames(df)[2]
    
    # Reorder columns for output
    MLR.results_df <- MLR.results_df[, c(5, 6, 1, 4, 2, 3)]
    return(MLR.results_df)
  }

  # Choose the appropriate process function based on the 'covariate' flag
  process_combination <- if(covariate) process_combination_cov else process_combination_no_cov
  
  # Run the regression in parallel across all combinations
  MLR.results <- foreach(i = 1:nrow(combinations), .combine = rbind, .options.snow = opts) %dopar% {
    process_combination(combinations[i,])
  }

  # Close the progress bar and stop the parallel cluster
  close(pb)
  stopCluster(cl)
  
  return(MLR.results)
}

# 02. Running the Multiple Linear Regression on different datasets
#####

# Example 1: MLR between df.01 and df.02 without covariates
MLR_12 <- meta.MLR.calculator(X = df.01, Y = df.02, covariate = FALSE)
MLR_12$p.adj_BH <- p.adjust(MLR_12$Pval, method = "BH") # Adjust p-values
MLR_12$x_class <- "1"
MLR_12$y_class <- "2"

# Example 2: MLR between df.02 and df.03 without covariates
MLR_23 <- meta.MLR.calculator(X = df.02, Y = df.03, covariate = FALSE)
MLR_23$p.adj_BH <- p.adjust(MLR_23$Pval, method = "BH")
MLR_23$x_class <- "2"
MLR_23$y_class <- "3"

# Example 3: MLR between df.01 and df.03 without covariates
MLR_13 <- meta.MLR.calculator(X = df.01, Y = df.03, covariate = FALSE)
MLR_13$p.adj_BH <- p.adjust(MLR_13$Pval, method = "BH")
MLR_13$x_class <- "1"
MLR_13$y_class <- "3"

# Combine all results into one dataframe and filter significant results (Pval <= 0.05)
asso_results <- rbind(MLR_12, MLR_23, MLR_13)
write.csv(asso_results, paste0(Sys.Date(), "_MLR_results.csv"))
asso_results <- asso_results %>% filter(Pval <= 0.05)
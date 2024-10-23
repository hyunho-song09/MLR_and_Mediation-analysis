rm(list=ls())
setwd("C:/Users/HHSONG/Desktop/programming/R_code/Help/JHJ/241022_Interaction_mediation")

# input
#####

library(xlsx)

raw.meta.df <- read.xlsx("input/241022_Muscle_Serum_example.xlsx", sheetName = "Metadata")
raw.muscle.df <- read.xlsx("input/241022_Muscle_Serum_example.xlsx", sheetName = "muscle_input")
raw.serum.df <- read.xlsx("input/241022_Muscle_Serum_example.xlsx", sheetName = "serum_input")

raw.meta.df$serum_ox[raw.meta.df$serum_ox == "o"] <- 1
raw.meta.df$serum_ox[raw.meta.df$serum_ox == "x"] <- 0
raw.meta.df$muscle_ox[raw.meta.df$muscle_ox == "o"] <- 1
raw.meta.df$muscle_ox[raw.meta.df$muscle_ox == "x"] <- 0

#####


# interaction analysis
#####
meta.col <- colnames(raw.meta.df[,2:9])
muscle.col <- colnames(raw.muscle.df[,2:ncol(raw.muscle.df)])
serum.col <- colnames(raw.serum.df[,2:ncol(raw.serum.df)])

interaction.df <- cbind(raw.meta.df[,2:9],
                        raw.muscle.df[,2:ncol(raw.muscle.df)],
                        raw.serum.df[,2:ncol(raw.serum.df)])

interaction.results <- list()

for (y in 1:length(meta.col)) {
  for (x2 in 1:length(serum.col)) {
    for (x1 in 1:length(muscle.col)) {
      formula <- paste0(y, "~", x1, "+", x2, "+", x1, "*", x2)
      interaction.results[[i]] <- lm(formula, data = interaction.df)
      names(interaction.results)[i] <- paste0(y, "_", x1, "_", x2)
    }
  }
}


#####


# 02-1. Running the Multiple Linear Regression on different datasets
#####

# 02-1. Multiple Linear Regression (for multiple features)
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

#####

# 02-2. mediation analysis
#####

# 02-2. Mediation Analysis with Parallel Processing
# This script performs mediation analysis across multiple combinations of independent (X), mediator (M),
# and dependent (Y) variables, with optional covariates (C). Parallel processing and progress bars 
# are used for computational efficiency.

# Required Libraries
library(dplyr)
library(doParallel)
library(doSNOW)
library(mediation)

# Function for Mediation Analysis
meta.mediate.calculator <- function(X, M, Y, C = NULL, covariate = FALSE) {
  print("Running Mediation Analysis...")
  
  # Create all combinations of columns between X, M, and Y
  combinations <- expand.grid(X = 1:ncol(X), M = 1:ncol(M), Y = 1:ncol(Y))
  
  # Register parallel backend
  n_cores <- detectCores() - 6  # Leave 6 cores free for other tasks
  cl <- makeCluster(n_cores)
  registerDoSNOW(cl)
  
  # Set up progress bar
  pb <- txtProgressBar(max = nrow(combinations), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # Helper function to extract mediation summary
  extract_mediation_summary <- function(x) {
    clp <- 100 * x$conf.level
    isLinear.y <- ((class(x$model.y)[1] %in% c("lm", "rq")) || 
                     (inherits(x$model.y, "glm") && x$model.y$family$family == 
                        "gaussian" && x$model.y$family$link == "identity") || 
                     (inherits(x$model.y, "survreg") && x$model.y$dist == 
                        "gaussian"))
    
    printone <- !x$INT && isLinear.y
    if (printone) {
      smat <- rbind(
        c(x$d1, x$d1.ci, x$d1.p),
        c(x$z0, x$z0.ci, x$z0.p),
        c(x$tau.coef, x$tau.ci, x$tau.p),
        c(x$n0, x$n0.ci, x$n0.p)
      )
      rownames(smat) <- c("ACME", "ADE", "TE", "PROP")
    } else {
      smat <- rbind(
        c(x$d0, x$d0.ci, x$d0.p),
        c(x$d1, x$d1.ci, x$d1.p),
        c(x$z0, x$z0.ci, x$z0.p),
        c(x$z1, x$z1.ci, x$z1.p),
        c(x$tau.coef, x$tau.ci, x$tau.p),
        c(x$n0, x$n0.ci, x$n0.p),
        c(x$n1, x$n1.ci, x$n1.p),
        c(x$d.avg, x$d.avg.ci, x$d.avg.p),
        c(x$z.avg, x$z.avg.ci, x$z.avg.p),
        c(x$n.avg, x$n.avg.ci, x$n.avg.p)
      )
      rownames(smat) <- c("ACME (control)", "ACME (treated)", "ADE (control)", 
                          "ADE (treated)", "TE", "PROP (control)", "PROP (treated)", 
                          "ACME (average)", "ADE (average)", "PROP (average)")
    }
    colnames(smat) <- c("Estimate", paste(clp, "% CI Lower"), paste(clp, "% CI Upper"), "p-value")
    return(smat)
  }
  
  # Main function to process each combination
  process_combination <- function(comb) {
    set.seed(2023)
    x <- comb$X
    m <- comb$M
    y <- comb$Y
    
    # Build the dataset with or without covariates
    if (covariate) {
      df <- na.omit(cbind(X[, x, drop = FALSE], M[, m, drop = FALSE], Y[, y, drop = FALSE], C))
      cov_formula <- paste0("+", "`", colnames(C), "`", collapse = "")
    } else {
      df <- na.omit(cbind(X[, x, drop = FALSE], M[, m, drop = FALSE], Y[, y, drop = FALSE]))
      cov_formula <- ""
    }
    
    # Fit mediation and outcome models
    med.fit <- glm(paste0("`", colnames(df)[2], "` ~ `", colnames(df)[1], "`", cov_formula), data = df)
    out.fit <- glm(paste0("`", colnames(df)[3], "` ~ `", colnames(df)[2], "` + `", colnames(df)[1], "`", cov_formula), data = df)
    
    # Perform mediation analysis
    med.out <- mediate(med.fit, out.fit, treat = colnames(df)[1], mediator = colnames(df)[2], robustSE = FALSE, sims = 1000)
    
    # Extract results from mediation summary
    med.results.tmp.01 <- extract_mediation_summary(summary(med.out))[, c(1, 4, 2, 3)]
    med.results.tmp.02 <- data.frame(cbind(
      t(med.results.tmp.01[1,]), t(med.results.tmp.01[2,]),
      t(med.results.tmp.01[3,]), t(med.results.tmp.01[4,])
    ))
    
    # Naming columns for easier interpretation
    colnames(med.results.tmp.02) <- c(
      paste0(rownames(med.results.tmp.01)[1], "es"), paste0(rownames(med.results.tmp.01)[1], "p"),
      paste0(rownames(med.results.tmp.01)[1], "ci1"), paste0(rownames(med.results.tmp.01)[1], "ci2"),
      paste0(rownames(med.results.tmp.01)[2], "es"), paste0(rownames(med.results.tmp.01)[2], "p"),
      paste0(rownames(med.results.tmp.01)[2], "ci1"), paste0(rownames(med.results.tmp.01)[2], "ci2"),
      paste0(rownames(med.results.tmp.01)[3], "es"), paste0(rownames(med.results.tmp.01)[3], "p"),
      paste0(rownames(med.results.tmp.01)[3], "ci1"), paste0(rownames(med.results.tmp.01)[3], "ci2"),
      paste0(rownames(med.results.tmp.01)[4], "es"), paste0(rownames(med.results.tmp.01)[4], "p"),
      paste0(rownames(med.results.tmp.01)[4], "ci1"), paste0(rownames(med.results.tmp.01)[4], "ci2")
    )
    
    # Adding variable names to the results
    med.results.tmp.02$X <- colnames(df)[1]
    med.results.tmp.02$M <- colnames(df)[2]
    med.results.tmp.02$Y <- colnames(df)[3]
    
    # Reordering columns for the final output
    return(med.results.tmp.02[, c(17:19, 1:16)])
  }
  
  # Parallelize computation across all combinations
  med.results <- foreach(i = 1:nrow(combinations), .combine = rbind, .options.snow = opts) %dopar% {
    tryCatch({
      process_combination(combinations[i, ])
    }, error = function(e) {})
  }
  
  # Close progress bar and stop parallel cluster
  close(pb)
  stopCluster(cl)
  
  return(med.results)
}

# Example: Creating input dataframes for mediation analysis
nm <- unique(c(unlist(asso_results[["X"]]), unlist(asso_results[["Y"]])))
med.df.01 <- df.01[, colnames(df.01) %in% nm]
med.df.02 <- df.02[, colnames(df.02) %in% nm]
med.df.03 <- df.03[, colnames(df.03) %in% nm]

# Direction 1: X -> M -> Y (df1 -> df2 -> df3)
mediation.results.01 <- meta.mediate.calculator(X = med.df.01, M = med.df.02, Y = med.df.03, covariate = FALSE)
mediation.results.01$direction <- 1

# Direction 2: X -> M -> Y (df1 -> df3 -> df2)
mediation.results.02 <- meta.mediate.calculator(X = med.df.01, M = med.df.03, Y = med.df.02, covariate = FALSE)
mediation.results.02$direction <- 2

# Combine both mediation results
mediation.results <- rbind(mediation.results.01, mediation.results.02)

#####
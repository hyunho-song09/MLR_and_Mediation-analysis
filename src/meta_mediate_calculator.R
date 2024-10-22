# Mediation Analysis with Parallel Processing
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
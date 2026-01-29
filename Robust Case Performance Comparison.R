library(glmnet)
library(cluster)
library(ggplot2)
library(dplyr)
library(tidyr)
library(dynamicTreeCut)
library(SIS)
library(MASS)
library(glasso)
library(robustbase)  # For OGK estimator
library(rrcov)       # Alternative robust covariance

set.seed(123)
setwd("/Users/wanru.guo")

# 1. Load CORRUPTED simulation data
simulated_data <- readRDS("100_corrupted_simulations_p1000_varying_cor_robust.rds")[1:100]
n_sim <- length(simulated_data)

############################################################
## Stage 0: Robust group construction
############################################################

# Function for OGK estimator (Option B)
estimate_ogk_covariance <- function(X) {
  # Orthogonalized Gnanadesikan-Kettenring estimator
  tryCatch({
    # Using robustbase::covOGK
    ogk_fit <- robustbase::covOGK(X, sigmamu = robustbase::s_Qn)
    return(ogk_fit$cov)
  }, error = function(e) {
    # Fallback to MCD if OGK fails
    tryCatch({
      mcd_fit <- covMcd(X)
      return(mcd_fit$cov)
    }, error = function(e2) {
      # Final fallback to regular covariance
      warning("Robust covariance estimation failed, using regular covariance")
      return(cov(X))
    })
  })
}

# Function to compute sparse correlation from robust covariance
compute_robust_sparse_correlation <- function(Sigma_robust, rho) {
  tryCatch({
    # Apply graphical lasso to robust covariance
    glasso_fit <- glasso(Sigma_robust, rho = rho, penalize.diagonal = FALSE)
    Theta <- glasso_fit$wi
    
    # Ensure positive definiteness
    if (any(eigen(Theta)$values <= 0)) {
      Theta <- Theta + diag(0.01, nrow(Theta))
    }
    
    # Convert precision to correlation
    D <- diag(1/sqrt(diag(Theta)))
    C_sparse <- D %*% Theta %*% D
    diag(C_sparse) <- 1
    return(C_sparse)
  }, error = function(e) {
    warning("Graphical lasso failed, using robust correlation")
    D <- diag(1/sqrt(diag(Sigma_robust)))
    return(D %*% Sigma_robust %*% D)
  })
}

# Adjusted Rand Index (ARI) function
calculate_ARI <- function(groups1, groups2) {
  if (length(unique(groups1)) == 1 && length(unique(groups2)) == 1) {
    return(1)
  }
  
  tab <- table(groups1, groups2)
  n <- sum(tab)
  sum_comb_a <- sum(choose(rowSums(tab), 2))
  sum_comb_b <- sum(choose(colSums(tab), 2))
  sum_comb <- sum(choose(tab, 2))
  
  expected <- sum_comb_a * sum_comb_b / choose(n, 2)
  max_comb <- (sum_comb_a + sum_comb_b) / 2
  
  if (max_comb - expected == 0) {
    return(1)
  }
  
  (sum_comb - expected) / (max_comb - expected)
}

# Main robust group construction function
construct_robust_groups <- function(X, method = c("spearman", "ogk"),
                                    rho_grid = c(0.01, 0.05, 0.1, 0.2),
                                    deepSplit = 2,
                                    minClusterSize = 3) {
  
  method <- match.arg(method)
  
  if (method == "spearman") {
    # Option A: Spearman correlation
    C <- cor(X, method = "spearman")
    D <- as.dist(1 - abs(C))
    hc <- hclust(D, method = "average")
    
    # Dynamic tree cut
    groups <- cutreeDynamic(
      dendro = hc,
      distM = as.matrix(D),
      deepSplit = deepSplit,
      minClusterSize = minClusterSize,
      pamStage = FALSE,
      pamRespectsDendro = FALSE,
      verbose = 0
    )
    
    return(list(groups = groups, hc = hc, method = "spearman"))
    
  } else {
    # Option B: OGK-based sparse correlation
    # Step 1: Compute robust covariance
    Sigma_ogk <- estimate_ogk_covariance(X)
    
    # Step 2: Try different rho values, compute ARI
    best_ari <- -Inf
    best_rho <- rho_grid[1]
    best_groups <- NULL
    best_hc <- NULL
    
    for (rho in rho_grid) {
      # Compute sparse correlation
      C_sparse <- compute_robust_sparse_correlation(Sigma_ogk, rho)
      
      # Hierarchical clustering
      D_sparse <- as.dist(1 - abs(C_sparse))
      hc_sparse <- hclust(D_sparse, method = "average")
      
      # Dynamic tree cut
      groups_sparse <- cutreeDynamic(
        dendro = hc_sparse,
        distM = as.matrix(D_sparse),
        deepSplit = deepSplit,
        minClusterSize = minClusterSize,
        pamStage = FALSE,
        pamRespectsDendro = FALSE,
        verbose = 0
      )
      
      # For first rho, use as reference
      if (rho == rho_grid[1]) {
        reference_groups <- groups_sparse
        ari <- 1
      } else {
        # Compute ARI with reference clustering
        ari <- calculate_ARI(groups_sparse, reference_groups)
      }
      
      # Update best if better ARI
      if (ari > best_ari) {
        best_ari <- ari
        best_rho <- rho
        best_groups <- groups_sparse
        best_hc <- hc_sparse
      }
    }
    
    return(list(
      groups = best_groups,
      hc = best_hc,
      method = "ogk",
      rho_star = best_rho,
      ARI = best_ari
    ))
  }
}

############################################################
## Stage 1: Robust global group tests
############################################################

perform_robust_group_tests <- function(X, y, groups, alpha1, 
                                       robust_method = c("huber", "bisquare")) {
  
  robust_method <- match.arg(robust_method)
  unique_groups <- unique(groups)
  G1 <- integer(0)
  group_pvals <- numeric(length(unique_groups))
  
  for (k in seq_along(unique_groups)) {
    group_idx <- which(groups == unique_groups[k])
    
    if (length(group_idx) == 1) {
      # Single variable: use robust regression
      fit_rlm <- try(MASS::rlm(y ~ X[, group_idx]), silent = TRUE)
      if (!inherits(fit_rlm, "try-error")) {
        w <- fit_rlm$w
        fit_wls <- lm(y ~ X[, group_idx], weights = w)
        group_pvals[k] <- summary(fit_wls)$coefficients[2, 4]
      } else {
        group_pvals[k] <- 1
      }
      
    } else if (length(group_idx) > nrow(X)) {
      # Group too large
      group_pvals[k] <- 1
      
    } else {
      # Multiple variables: fit robust multivariate model
      Xg <- X[, group_idx, drop = FALSE]
      
      if (robust_method == "huber") {
        fit_rlm <- try(MASS::rlm(y ~ Xg), silent = TRUE)
      } else {
        # bisquare weighting
        fit_rlm <- try(MASS::rlm(y ~ Xg, psi = psi.bisquare), silent = TRUE)
      }
      
      if (!inherits(fit_rlm, "try-error")) {
        # Extract weights from robust fit
        w <- fit_rlm$w
        
        # Weighted least squares
        fit_wls <- try(lm(y ~ Xg, weights = w), silent = TRUE)
        if (!inherits(fit_wls, "try-error")) {
          # F-test for group significance
          null_model <- lm(y ~ 1, weights = w)
          anova_test <- anova(null_model, fit_wls)
          group_pvals[k] <- anova_test$`Pr(>F)`[2]
        } else {
          group_pvals[k] <- 1
        }
      } else {
        group_pvals[k] <- 1
      }
    }
    
    if (group_pvals[k] < alpha1) {
      G1 <- c(G1, group_idx)
    }
  }
  
  return(list(selected_vars = G1, group_pvals = group_pvals))
}

############################################################
## Stage 2: Robust within-group variable tests
############################################################

perform_robust_variable_tests <- function(X, y, groups, selected_group_vars, alpha2,
                                          robust_method = c("huber", "bisquare")) {
  
  robust_method <- match.arg(robust_method)
  unique_groups <- unique(groups)
  S <- integer(0)
  variable_pvals <- list()
  
  for (k in seq_along(unique_groups)) {
    group_idx <- which(groups == unique_groups[k])
    
    # Only test groups that passed Stage 1
    if (any(group_idx %in% selected_group_vars)) {
      pvals <- numeric(length(group_idx))
      
      for (j in seq_along(group_idx)) {
        var_idx <- group_idx[j]
        
        if (robust_method == "huber") {
          fit_rlm <- try(MASS::rlm(y ~ X[, var_idx]), silent = TRUE)
        } else {
          fit_rlm <- try(MASS::rlm(y ~ X[, var_idx], psi = psi.bisquare), silent = TRUE)
        }
        
        if (!inherits(fit_rlm, "try-error")) {
          w <- fit_rlm$w
          fit_wls <- try(lm(y ~ X[, var_idx], weights = w), silent = TRUE)
          if (!inherits(fit_wls, "try-error")) {
            pvals[j] <- summary(fit_wls)$coefficients[2, 4]
          } else {
            pvals[j] <- 1
          }
        } else {
          pvals[j] <- 1
        }
      }
      
      # Select variables within group
      selected_in_group <- group_idx[pvals < alpha2 & !is.na(pvals)]
      S <- c(S, selected_in_group)
      
      variable_pvals[[k]] <- data.frame(
        variable = group_idx,
        p_value = pvals,
        selected = pvals < alpha2
      )
    }
  }
  
  return(list(selected_vars = unique(S), variable_pvals = variable_pvals))
}

############################################################
## Stage 3: Penalized regression
############################################################

perform_penalized_regression <- function(X, y, selected_vars,
                                         method = c("EN", "Adaptive_EN", 
                                                    "LASSO", "Adaptive_LASSO"),
                                         gamma = 1) {
  
  if (length(selected_vars) == 0) {
    return(integer(0))
  }
  
  method <- match.arg(method)
  X_subset <- X[, selected_vars, drop = FALSE]
  
  if (grepl("Adaptive", method)) {
    # Adaptive methods: two-stage procedure
    alpha_val <- ifelse(grepl("EN", method), 0.5, 1)
    
    # First stage: initial estimate
    cv_init <- cv.glmnet(X_subset, y, alpha = alpha_val, nfolds = 5)
    beta_init <- as.numeric(coef(cv_init, s = "lambda.min"))[-1]
    
    # Adaptive weights
    eps <- 1e-6
    w <- 1 / (abs(beta_init)^gamma + eps)
    
    # Second stage: adaptive fit
    cv_adapt <- cv.glmnet(X_subset, y,
                          alpha = alpha_val,
                          penalty.factor = w,
                          nfolds = 5)
    
    beta_final <- as.numeric(coef(cv_adapt, s = "lambda.min"))[-1]
    final_selected <- selected_vars[beta_final != 0]
    
  } else {
    # Non-adaptive methods
    alpha_val <- ifelse(grepl("EN", method), 0.5, 1)
    
    cv_fit <- cv.glmnet(X_subset, y,
                        alpha = alpha_val,
                        nfolds = 5)
    
    beta_final <- as.numeric(coef(cv_fit, s = "lambda.min"))[-1]
    final_selected <- selected_vars[beta_final != 0]
  }
  
  return(final_selected)
}

############################################################
## Complete Robust Dorfman Screening Algorithm
############################################################

robust_dorfman_screening <- function(X, y,
                                     group_method = c("spearman", "ogk"),
                                     alpha1_grid = NULL,
                                     alpha2_grid = NULL,
                                     method = c("EN", "Adaptive_EN"),
                                     tune_cv = TRUE,
                                     nfolds = 5,
                                     deepSplit = 2,
                                     minClusterSize = 3,
                                     robust_method = "huber") {
  
  group_method <- match.arg(group_method)
  method <- match.arg(method)
  
  # Determine if using Vanilla or Adaptive EN
  is_adaptive <- grepl("Adaptive", method)
  
  # Set alpha grids based on method type (as per Algorithm 2)
  if (is.null(alpha1_grid)) {
    if (is_adaptive) {
      alpha1_grid <- c(0.1, 0.2, 0.3)  # Adaptive EN: {0.1, 0.2, 0.3}
    } else {
      alpha1_grid <- c(0.01, 0.05, 0.1, 0.2)  # Vanilla EN: {0.01, 0.05, 0.1, 0.2}
    }
  }
  
  if (is.null(alpha2_grid)) {
    if (is_adaptive) {
      alpha2_grid <- c(0.1, 0.15)  # Adaptive EN: {0.1, 0.15}
    } else {
      alpha2_grid <- c(0.01, 0.05)  # Vanilla EN: {0.01, 0.05}
    }
  }
  
  if (tune_cv) {
    # Tune parameters via 5-fold CV
    tuned_params <- tune_robust_dorfman_params(
      X, y,
      group_method = group_method,
      alpha1_grid = alpha1_grid,
      alpha2_grid = alpha2_grid,
      method = method,
      nfolds = nfolds,
      deepSplit = deepSplit,
      minClusterSize = minClusterSize,
      robust_method = robust_method
    )
    
    alpha1 <- tuned_params$alpha1
    alpha2 <- tuned_params$alpha2
    rho <- ifelse(group_method == "ogk", tuned_params$rho, NULL)
    
  } else {
    # Use default parameters
    alpha1 <- ifelse(is_adaptive, 0.2, 0.05)
    alpha2 <- ifelse(is_adaptive, 0.15, 0.01)
    rho <- 0.1
  }
  
  # Stage 0: Robust group construction
  group_result <- construct_robust_groups(
    X,
    method = group_method,
    deepSplit = deepSplit,
    minClusterSize = minClusterSize
  )
  groups <- group_result$groups
  
  # Stage 1: Robust global group tests
  group_test_result <- perform_robust_group_tests(
    X, y, groups, alpha1,
    robust_method = robust_method
  )
  selected_group_vars <- group_test_result$selected_vars
  
  # Stage 2: Robust within-group variable tests
  variable_test_result <- perform_robust_variable_tests(
    X, y, groups, selected_group_vars, alpha2,
    robust_method = robust_method
  )
  selected_vars <- variable_test_result$selected_vars
  
  # Stage 3: Penalized regression
  final_selected <- perform_penalized_regression(
    X, y, selected_vars,
    method = method
  )
  
  return(list(
    final_selected = final_selected,
    intermediate_selected = selected_vars,
    groups = groups,
    group_pvals = group_test_result$group_pvals,
    variable_pvals = variable_test_result$variable_pvals,
    parameters = list(
      group_method = group_method,
      alpha1 = alpha1,
      alpha2 = alpha2,
      method = method,
      rho = ifelse(group_method == "ogk", group_result$rho_star, NA),
      robust_method = robust_method
    )
  ))
}

############################################################
## Parameter tuning via cross-validation
############################################################

tune_robust_dorfman_params <- function(X, y,
                                       group_method = "spearman",
                                       alpha1_grid,
                                       alpha2_grid,
                                       method = "EN",
                                       nfolds = 5,
                                       deepSplit = 2,
                                       minClusterSize = 3,
                                       robust_method = "huber") {
  
  folds <- sample(rep(1:nfolds, length.out = nrow(X)))
  best_rmse <- Inf
  best_params <- NULL
  
  # Create parameter grid
  param_grid <- expand.grid(
    alpha1 = alpha1_grid,
    alpha2 = alpha2_grid
  )
  
  for (i in 1:nrow(param_grid)) {
    alpha1 <- param_grid$alpha1[i]
    alpha2 <- param_grid$alpha2[i]
    
    fold_rmses <- numeric(nfolds)
    
    for (fold in 1:nfolds) {
      train_idx <- which(folds != fold)
      test_idx <- which(folds == fold)
      
      # Stage 0: Construct groups on training data
      group_result <- construct_robust_groups(
        X[train_idx, ],
        method = group_method,
        deepSplit = deepSplit,
        minClusterSize = minClusterSize
      )
      groups <- group_result$groups
      
      # Stage 1: Group tests
      group_test_result <- perform_robust_group_tests(
        X[train_idx, ], y[train_idx],
        groups, alpha1,
        robust_method = robust_method
      )
      selected_group_vars <- group_test_result$selected_vars
      
      # Stage 2: Variable tests
      variable_test_result <- perform_robust_variable_tests(
        X[train_idx, ], y[train_idx],
        groups, selected_group_vars, alpha2,
        robust_method = robust_method
      )
      selected_vars <- variable_test_result$selected_vars
      
      if (length(selected_vars) > 0) {
        # Stage 3: Penalized regression
        final_selected <- perform_penalized_regression(
          X[train_idx, ], y[train_idx],
          selected_vars,
          method = method
        )
        
        if (length(final_selected) > 0) {
          # Predict on test set
          X_test_selected <- X[test_idx, final_selected, drop = FALSE]
          
          if (length(final_selected) == 1) {
            model <- lm(y[train_idx] ~ X[train_idx, final_selected])
            pred <- predict(model, data.frame(X_test_selected))
          } else {
            model <- lm(y[train_idx] ~ .,
                        data = as.data.frame(X[train_idx, final_selected, drop = FALSE]))
            pred <- predict(model, as.data.frame(X_test_selected))
          }
          
          fold_rmses[fold] <- sqrt(mean((y[test_idx] - pred)^2))
        } else {
          fold_rmses[fold] <- sqrt(mean((y[test_idx] - mean(y[train_idx]))^2))
        }
      } else {
        fold_rmses[fold] <- sqrt(mean((y[test_idx] - mean(y[train_idx]))^2))
      }
    }
    
    avg_rmse <- mean(fold_rmses)
    
    if (avg_rmse < best_rmse) {
      best_rmse <- avg_rmse
      best_params <- list(
        alpha1 = alpha1,
        alpha2 = alpha2,
        rmse = avg_rmse
      )
    }
  }
  
  # For OGK method, also tune rho
  if (group_method == "ogk") {
    best_ari <- -Inf
    best_rho <- 0.1
    
    for (rho in c(0.01, 0.05, 0.1, 0.2)) {
      # Compute ARI across different parameter combinations
      ari_vals <- numeric(nrow(param_grid))
      
      for (i in 1:nrow(param_grid)) {
        alpha1 <- param_grid$alpha1[i]
        
        # Construct groups with current rho
        Sigma_ogk <- estimate_ogk_covariance(X)
        C_sparse <- compute_robust_sparse_correlation(Sigma_ogk, rho)
        D_sparse <- as.dist(1 - abs(C_sparse))
        hc_sparse <- hclust(D_sparse, method = "average")
        
        groups_sparse <- cutreeDynamic(
          dendro = hc_sparse,
          distM = as.matrix(D_sparse),
          deepSplit = deepSplit,
          minClusterSize = minClusterSize,
          pamStage = FALSE,
          pamRespectsDendro = FALSE,
          verbose = 0
        )
        
        # Compare with reference (first parameter combination)
        if (i == 1) {
          reference_groups <- groups_sparse
          ari_vals[i] <- 1
        } else {
          ari_vals[i] <- calculate_ARI(groups_sparse, reference_groups)
        }
      }
      
      avg_ari <- mean(ari_vals, na.rm = TRUE)
      
      if (avg_ari > best_ari) {
        best_ari <- avg_ari
        best_rho <- rho
      }
    }
    
    best_params$rho <- best_rho
  }
  
  return(best_params)
}

############################################################
## Specific Robust Dorfman Method Implementations
############################################################

# 1. Robust Dorfman with Spearman correlation + Vanilla EN
run_robust_dorfman_spearman_en <- function(X, y) {
  result <- robust_dorfman_screening(
    X, y,
    group_method = "spearman",
    method = "EN"
  )
  return(result$final_selected)
}

# 2. Robust Dorfman with Spearman correlation + Adaptive EN
run_robust_dorfman_spearman_adaptive_en <- function(X, y) {
  result <- robust_dorfman_screening(
    X, y,
    group_method = "spearman",
    method = "Adaptive_EN"
  )
  return(result$final_selected)
}

# 3. Robust Dorfman with OGK + Vanilla EN
run_robust_dorfman_ogk_en <- function(X, y) {
  result <- robust_dorfman_screening(
    X, y,
    group_method = "ogk",
    method = "EN"
  )
  return(result$final_selected)
}

# 4. Robust Dorfman with OGK + Adaptive EN
run_robust_dorfman_ogk_adaptive_en <- function(X, y) {
  result <- robust_dorfman_screening(
    X, y,
    group_method = "ogk",
    method = "Adaptive_EN"
  )
  return(result$final_selected)
}

############################################################
## Comparison methods (for baseline)
############################################################

run_en <- function(X, y) {
  cv_fit <- cv.glmnet(X, y, alpha = 0.5, nfolds = 5)
  which(coef(cv_fit, s = "lambda.min")[-1] != 0)
}

run_adaptive_en <- function(X, y, gamma = 1) {
  init <- cv.glmnet(X, y, alpha = 0.5, nfolds = 5)
  beta_init <- as.numeric(coef(init, s = "lambda.min"))[-1]
  w <- 1 / (abs(beta_init) + 1e-6)^gamma
  
  cv_adapt <- cv.glmnet(X, y,
                        alpha = 0.5,
                        penalty.factor = w,
                        nfolds = 5)
  which(coef(cv_adapt, s = "lambda.min")[-1] != 0)
}

run_lasso <- function(X, y) {
  cv_fit <- cv.glmnet(X, y, alpha = 1, nfolds = 5)
  which(coef(cv_fit, s = "lambda.min")[-1] != 0)
}

run_adaptive_lasso <- function(X, y, gamma = 1) {
  init <- cv.glmnet(X, y, alpha = 1, nfolds = 5)
  beta_init <- as.numeric(coef(init, s = "lambda.min"))[-1]
  w <- 1 / (abs(beta_init) + 1e-6)^gamma
  
  cv_adapt <- cv.glmnet(X, y,
                        alpha = 1,
                        penalty.factor = w,
                        nfolds = 5)
  which(coef(cv_adapt, s = "lambda.min")[-1] != 0)
}

run_sis_lasso <- function(X, y) {
  sis_fit <- SIS(X, y, family = "gaussian", tune = "bic")
  sis_fit$sis.ix0
}

############################################################
## Metrics calculation
############################################################

calculate_metrics <- function(selected_vars, true_vars,
                              X_train, y_train, X_test, y_test) {
  true_pos  <- sum(selected_vars %in% true_vars)
  false_pos <- sum(!(selected_vars %in% true_vars))
  false_neg <- sum(!(true_vars %in% selected_vars))
  
  TPR <- ifelse((true_pos + false_neg) > 0,
                true_pos / (true_pos + false_neg), 0)
  FDR <- ifelse((true_pos + false_pos) > 0,
                false_pos / (true_pos + false_pos), 0)
  PPV <- ifelse((true_pos + false_pos) > 0,
                true_pos / (true_pos + false_pos), 1)
  F1  <- ifelse((PPV + TPR) > 0,
                2 * (PPV * TPR) / (PPV + TPR), 0)
  
  if (length(selected_vars) > 0) {
    model <- lm(y_train ~ .,
                data.frame(X_train[, selected_vars, drop = FALSE]))
    pred  <- predict(model,
                     data.frame(X_test[, selected_vars, drop = FALSE]))
    rmse  <- sqrt(mean((y_test - pred)^2))
  } else {
    pred   <- mean(y_train)
    rmse   <- sqrt(mean((y_test - pred)^2))
  }
  
  data.frame(
    Vars  = length(selected_vars),
    RMSE  = rmse,
    TPR   = TPR,
    FDR   = FDR,
    F1    = F1
  )
}

############################################################
## Main simulation loop for Robust Scenario
############################################################

methods <- c(
  # Robust Dorfman methods (Vanilla EN)
  "Robust_Dorfman_Spearman_EN",
  "Robust_Dorfman_OGK_EN",
  
  # Robust Dorfman methods (Adaptive EN)
  "Robust_Dorfman_Spearman_Adaptive_EN",
  "Robust_Dorfman_OGK_Adaptive_EN",
  
  # Baseline methods
  "Adaptive_LASSO",
  "Adaptive_EN",
  "LASSO",
  "EN",
  "SIS_LASSO"
)

results <- list()
timing_results <- list()

for (sim in 1:n_sim) {
  dat  <- simulated_data[[sim]]
  start_time <- Sys.time()
  
  X    <- dat$X
  y    <- dat$y
  true_vars <- which(dat$beta != 0)
  
  # Split data (60% train, 40% test)
  train_idx <- sample(nrow(X), 0.6 * nrow(X))
  X_train   <- X[train_idx, ]
  y_train   <- y[train_idx]
  X_test    <- X[-train_idx, ]
  y_test    <- y[-train_idx]
  
  # Run all methods
  selected_methods <- list(
    # Robust Dorfman Vanilla EN
    tryCatch(run_robust_dorfman_spearman_en(X_train, y_train),
             error = function(e) integer(0)),
    tryCatch(run_robust_dorfman_ogk_en(X_train, y_train),
             error = function(e) integer(0)),
    
    # Robust Dorfman Adaptive EN
    tryCatch(run_robust_dorfman_spearman_adaptive_en(X_train, y_train),
             error = function(e) integer(0)),
    tryCatch(run_robust_dorfman_ogk_adaptive_en(X_train, y_train),
             error = function(e) integer(0)),
    
    # Baseline methods
    tryCatch(run_adaptive_lasso(X_train, y_train),
             error = function(e) integer(0)),
    tryCatch(run_adaptive_en(X_train, y_train),
             error = function(e) integer(0)),
    tryCatch(run_lasso(X_train, y_train),
             error = function(e) integer(0)),
    tryCatch(run_en(X_train, y_train),
             error = function(e) integer(0)),
    tryCatch(run_sis_lasso(X_train, y_train),
             error = function(e) integer(0))
  )
  
  # Calculate metrics
  metrics <- lapply(selected_methods, function(sel) {
    calculate_metrics(sel, true_vars, X_train, y_train, X_test, y_test)
  })
  
  # Store results
  for (i in seq_along(methods)) {
    results[[paste0("Sim", sim, "_", methods[i])]] <-
      metrics[[i]] %>% mutate(Simulation = sim, Method = methods[i])
  }
  
  end_time <- Sys.time()
  timing_results[[sim]] <- end_time - start_time
  
  cat(sprintf("Simulation %2d completed\n", sim))
}

# Process results
results_df <- bind_rows(results)

# Performance summary
performance_summary <- results_df %>%
  group_by(Method) %>%
  summarise(
    N_Sims = n(),
    Mean_Vars = round(mean(Vars, na.rm = TRUE), 1),
    SD_Vars = round(sd(Vars, na.rm = TRUE), 2),
    Mean_RMSE = round(mean(RMSE, na.rm = TRUE), 3),
    SD_RMSE = round(sd(RMSE, na.rm = TRUE), 3),
    Mean_TPR = round(mean(TPR, na.rm = TRUE), 3),
    SD_TPR = round(sd(TPR, na.rm = TRUE), 3),
    Mean_FDR = round(mean(FDR, na.rm = TRUE), 3),
    SD_FDR = round(sd(FDR, na.rm = TRUE), 3),
    Mean_F1 = round(mean(F1, na.rm = TRUE), 3),
    SD_F1 = round(sd(F1, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_F1))

cat("\n=== ROBUST SCENARIO PERFORMANCE SUMMARY ===\n")
print(performance_summary)

# Visualization
ggplot(results_df, aes(x = Method, y = F1, fill = Method)) +
  geom_boxplot() +
  labs(title = "F1 Score Comparison: Robust Scenario (Corrupted Data)",
       subtitle = "Methods using Spearman/OGK correlation + Huber robust regression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(results_df, aes(x = Method, y = RMSE, fill = Method)) +
  geom_boxplot() +
  labs(title = "RMSE Comparison: Robust Scenario (Corrupted Data)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Compare Vanilla vs Adaptive EN within Robust Dorfman
robust_comparison <- results_df %>%
  filter(grepl("Robust_Dorfman", Method)) %>%
  mutate(
    group_method = ifelse(grepl("Spearman", Method), "Spearman", "OGK"),
    penalty_type = ifelse(grepl("Adaptive", Method), "Adaptive_EN", "Vanilla_EN")
  )

ggplot(robust_comparison, aes(x = group_method, y = F1, fill = penalty_type)) +
  geom_boxplot() +
  labs(title = "Robust Dorfman: Group Method vs Penalty Type",
       x = "Group Construction Method",
       y = "F1 Score",
       fill = "Penalty Type") +
  theme_minimal()

# Final table for paper
final_table <- performance_summary %>%
  mutate(
    Variables = sprintf("%.1f (%.2f)", Mean_Vars, SD_Vars),
    TPR = sprintf("%.3f (%.3f)", Mean_TPR, SD_TPR),
    FDR = sprintf("%.3f (%.3f)", Mean_FDR, SD_FDR),
    F1 = sprintf("%.3f (%.3f)", Mean_F1, SD_F1),
    RMSE = sprintf("%.3f (%.3f)", Mean_RMSE, SD_RMSE)
  ) %>%
  select(Method, Variables, TPR, FDR, F1, RMSE) %>%
  arrange(desc(gsub(".*\\(|\\)", "", F1)))  # Sort by F1 score

cat("\n=== FINAL PERFORMANCE TABLE FOR PAPER ===\n")
print(final_table)

# Timing analysis
timing_df <- data.frame(
  Simulation = 1:n_sim,
  Time_seconds = sapply(timing_results, as.numeric)
)

cat("\n=== TIMING ANALYSIS ===\n")
cat("Average time per simulation:", round(mean(timing_df$Time_seconds), 2), "seconds\n")
cat("Total time for all simulations:", round(sum(timing_df$Time_seconds), 2), "seconds\n")



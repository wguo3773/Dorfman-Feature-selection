library(glmnet)
library(cluster)
library(ggplot2)
library(dplyr)
library(tidyr)
library(dynamicTreeCut)
library(SIS)
library(glasso)  # For graphical lasso
library(spcov)   # Alternative for sparse covariance estimation
library(caret)   # For cross-validation

set.seed(123)
setwd("/Users/wanru.guo")

# 1. Load simulation data
simulated_data <- readRDS("100_simulations_p1000_noise1.rds")[1:100]  
n_sim <- length(simulated_data)

# 2. Tuning parameters according to algorithm specifications
# h: dendrogram cut height (tuned via CV)
# rho: graphical lasso penalty (selected from grid)
# alpha1: group significance threshold (Normal EN: {0.05, 0.1, 0.2}, Adaptive EN: {0.2, 0.3, 0.4})
# alpha2: variable significance threshold (Normal EN: 0.05, Adaptive EN: {0.15, 0.2})

############################################################
## Helper: Graphical Lasso for sparse correlation
############################################################
compute_sparse_correlation <- function(X, rho, use_glasso = TRUE) {
  # Compute Pearson covariance matrix
  S <- cov(X)
  
  if (use_glasso) {
    # Apply graphical lasso
    tryCatch({
      glasso_fit <- glasso(S, rho = rho, penalize.diagonal = FALSE)
      # Convert precision matrix to correlation matrix
      Theta <- glasso_fit$wi
      # Ensure Theta is positive definite
      if (any(eigen(Theta)$values <= 0)) {
        Theta <- Theta + diag(0.01, nrow(Theta))
      }
      # Convert precision to correlation
      D <- diag(1/sqrt(diag(Theta)))
      C_sparse <- D %*% Theta %*% D
      diag(C_sparse) <- 1
      return(C_sparse)
    }, error = function(e) {
      # Fallback to regular correlation if glasso fails
      warning("Graphical lasso failed, using regular correlation")
      return(cor(X))
    })
  } else {
    # Alternative: Threshold regular correlation
    C <- cor(X)
    C[abs(C) < rho] <- 0
    diag(C) <- 1
    return(C)
  }
}

############################################################
## Helper: ARI calculation for tuning rho
############################################################
calculate_ARI <- function(groups1, groups2) {
  # Adjusted Rand Index for comparing two clusterings
  if (length(unique(groups1)) == 1 && length(unique(groups2)) == 1) {
    return(1)  # Both have only one cluster
  }
  
  # Create contingency table
  tab <- table(groups1, groups2)
  
  # Calculate ARI
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

############################################################
## Stage 0: Group construction with Option A/B
############################################################
construct_groups <- function(X, h, rho = NULL, use_sparse = FALSE, return_ARI = FALSE) {
  if (use_sparse && !is.null(rho)) {
    # Option B: Sparse correlation refinement with graphical lasso
    C_sparse <- compute_sparse_correlation(X, rho)
    D <- as.dist(1 - abs(C_sparse))
  } else {
    # Option A: Pearson correlation
    C <- cor(X)
    D <- as.dist(1 - abs(C))
  }
  
  # Hierarchical clustering with average linkage
  hc <- hclust(D, method = "average")
  
  # Cut dendrogram at height h
  groups <- cutree(hc, h = h)
  
  if (return_ARI) {
    # For tuning rho, compute ARI between sparse and regular correlation clustering
    D_regular <- as.dist(1 - abs(cor(X)))
    hc_regular <- hclust(D_regular, method = "average")
    groups_regular <- cutree(hc_regular, h = h)
    
    ari <- calculate_ARI(groups, groups_regular)
    return(list(groups = groups, ARI = ari))
  }
  
  return(list(groups = groups, hc = hc))
}

############################################################
## Stage 1: Global group tests
############################################################
perform_group_tests <- function(X, y, groups, alpha1) {
  unique_groups <- unique(groups)
  G1 <- integer(0)  # Selected groups
  group_pvals <- numeric(length(unique_groups))
  
  for (k in seq_along(unique_groups)) {
    group_idx <- which(groups == unique_groups[k])
    
    if (length(group_idx) == 1) {
      # Single variable group
      model <- lm(y ~ X[, group_idx])
      group_pvals[k] <- summary(model)$coefficients[2, 4]
    } else if (length(group_idx) > nrow(X)) {
      # Group too large for regression
      group_pvals[k] <- 1  # Conservative: don't select
    } else {
      # Fit group model
      tryCatch({
        group_data <- data.frame(y = y, X = X[, group_idx, drop = FALSE])
        model <- lm(y ~ ., data = group_data)
        
        # F-test for group significance
        if (length(group_idx) > 1) {
          null_model <- lm(y ~ 1)
          anova_test <- anova(null_model, model)
          group_pvals[k] <- anova_test$`Pr(>F)`[2]
        } else {
          group_pvals[k] <- summary(model)$coefficients[2, 4]
        }
      }, error = function(e) {
        group_pvals[k] <- 1  # Error in fitting
      })
    }
    
    if (group_pvals[k] < alpha1) {
      G1 <- c(G1, group_idx)
    }
  }
  
  return(list(selected_vars = G1, group_pvals = group_pvals))
}

############################################################
## Stage 2: Within-group variable tests
############################################################
perform_variable_tests <- function(X, y, groups, selected_group_vars, alpha2) {
  unique_groups <- unique(groups)
  S <- integer(0)  # Final selected variables
  variable_pvals <- list()
  
  for (k in seq_along(unique_groups)) {
    group_idx <- which(groups == unique_groups[k])
    
    # Only test within groups that had variables selected in Stage 1
    if (any(group_idx %in% selected_group_vars)) {
      pvals <- numeric(length(group_idx))
      
      for (j in seq_along(group_idx)) {
        var_idx <- group_idx[j]
        tryCatch({
          model <- lm(y ~ X[, var_idx])
          pvals[j] <- summary(model)$coefficients[2, 4]
        }, error = function(e) {
          pvals[j] <- 1  # Error in fitting
        })
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
    # Adaptive methods: need initial estimate
    alpha_val <- ifelse(grepl("EN", method), 0.5, 1)
    
    # First stage: get initial estimate
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
## Complete Dorfman Screening Algorithm
############################################################
dorfman_screening <- function(X, y, 
                              h = NULL, 
                              rho_grid = c(0.01, 0.05, 0.1, 0.2),
                              alpha1_grid = c(0.05, 0.1, 0.2),
                              alpha2 = 0.05,
                              method = "EN",
                              tune_cv = TRUE,
                              nfolds = 5) {
  
  # Determine if using Normal or Adaptive version
  is_adaptive <- grepl("Adaptive", method)
  
  # Set alpha grids based on method type
  if (is_adaptive) {
    if (missing(alpha1_grid)) {
      alpha1_grid <- c(0.2, 0.3, 0.4)
    }
    if (missing(alpha2)) {
      alpha2 <- c(0.15, 0.2)
    }
  }
  
  if (tune_cv) {
    # Tune parameters via cross-validation
    tuned_params <- tune_dorfman_params(X, y, 
                                        rho_grid = rho_grid,
                                        alpha1_grid = alpha1_grid,
                                        alpha2_grid = alpha2,
                                        method = method,
                                        nfolds = nfolds)
    
    h <- tuned_params$h
    rho <- tuned_params$rho
    alpha1 <- tuned_params$alpha1
    alpha2 <- tuned_params$alpha2
    use_sparse <- tuned_params$use_sparse
  } else {
    # Use provided parameters
    use_sparse <- !is.null(rho)
    if (is.null(h)) h <- 0.5
    if (use_sparse) rho <- rho_grid[1]
    if (is.null(alpha1)) alpha1 <- ifelse(is_adaptive, 0.3, 0.1)
    if (is.null(alpha2)) alpha2 <- ifelse(is_adaptive, 0.15, 0.05)
  }
  
  # Stage 0: Group construction
  group_result <- construct_groups(X, h = h, rho = rho, use_sparse = use_sparse)
  groups <- group_result$groups
  
  # Stage 1: Global group tests
  group_test_result <- perform_group_tests(X, y, groups, alpha1)
  selected_group_vars <- group_test_result$selected_vars
  
  # Stage 2: Within-group variable tests
  variable_test_result <- perform_variable_tests(X, y, groups, selected_group_vars, alpha2)
  selected_vars <- variable_test_result$selected_vars
  
  # Stage 3: Penalized regression
  final_selected <- perform_penalized_regression(X, y, selected_vars, method = method)
  
  return(list(
    final_selected = final_selected,
    intermediate_selected = selected_vars,
    groups = groups,
    group_pvals = group_test_result$group_pvals,
    variable_pvals = variable_test_result$variable_pvals,
    parameters = list(h = h, rho = ifelse(use_sparse, rho, NA), 
                      alpha1 = alpha1, alpha2 = alpha2,
                      method = method)
  ))
}

############################################################
## Parameter tuning via cross-validation
############################################################
tune_dorfman_params <- function(X, y, 
                                rho_grid = c(0.01, 0.05, 0.1, 0.2),
                                alpha1_grid = c(0.05, 0.1, 0.2),
                                alpha2_grid = 0.05,
                                method = "EN",
                                nfolds = 5) {
  
  is_adaptive <- grepl("Adaptive", method)
  
  # Create folds
  folds <- sample(rep(1:nfolds, length.out = nrow(X)))
  
  # Grid for tuning
  h_grid <- seq(0.1, 0.9, by = 0.2)
  param_grid <- expand.grid(
    h = h_grid,
    rho = c(NA, rho_grid),  # NA means use Option A (regular correlation)
    alpha1 = alpha1_grid,
    alpha2 = alpha2_grid
  )
  
  best_rmse <- Inf
  best_params <- NULL
  
  for (i in 1:nrow(param_grid)) {
    h <- param_grid$h[i]
    rho <- param_grid$rho[i]
    alpha1 <- param_grid$alpha1[i]
    alpha2 <- param_grid$alpha2[i]
    
    fold_rmses <- numeric(nfolds)
    
    for (fold in 1:nfolds) {
      train_idx <- which(folds != fold)
      test_idx <- which(folds == fold)
      
      # Use sparse correlation if rho is not NA
      use_sparse <- !is.na(rho)
      
      # Stage 0: Group construction
      group_result <- construct_groups(X[train_idx, ], h = h, 
                                       rho = ifelse(use_sparse, rho, NULL), 
                                       use_sparse = use_sparse)
      groups <- group_result$groups
      
      # Stage 1: Group tests
      group_test_result <- perform_group_tests(X[train_idx, ], y[train_idx], 
                                               groups, alpha1)
      selected_group_vars <- group_test_result$selected_vars
      
      # Stage 2: Variable tests
      variable_test_result <- perform_variable_tests(X[train_idx, ], y[train_idx], 
                                                     groups, selected_group_vars, alpha2)
      selected_vars <- variable_test_result$selected_vars
      
      if (length(selected_vars) > 0) {
        # Stage 3: Penalized regression
        final_selected <- perform_penalized_regression(X[train_idx, ], y[train_idx], 
                                                       selected_vars, method = method)
        
        if (length(final_selected) > 0) {
          # Predict on test set
          X_test_selected <- X[test_idx, final_selected, drop = FALSE]
          if (length(final_selected) == 1) {
            pred <- predict(lm(y[train_idx] ~ X[train_idx, final_selected]),
                            data.frame(X_test_selected))
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
        h = h,
        rho = rho,
        alpha1 = alpha1,
        alpha2 = alpha2,
        use_sparse = use_sparse,
        rmse = avg_rmse
      )
    }
  }
  
  # If sparse correlation wasn't best, try tuning rho separately
  if (is.na(best_params$rho)) {
    # Tune rho using ARI criterion
    best_ari <- -Inf
    best_rho <- rho_grid[1]
    
    for (rho_val in rho_grid) {
      ari_vals <- numeric(length(h_grid))
      for (j in seq_along(h_grid)) {
        group_result <- construct_groups(X, h = h_grid[j], 
                                         rho = rho_val, 
                                         use_sparse = TRUE,
                                         return_ARI = TRUE)
        ari_vals[j] <- group_result$ARI
      }
      avg_ari <- mean(ari_vals, na.rm = TRUE)
      
      if (avg_ari > best_ari) {
        best_ari <- avg_ari
        best_rho <- rho_val
      }
    }
    
    best_params$rho <- best_rho
    best_params$use_sparse <- TRUE
  }
  
  return(best_params)
}

############################################################
## Specific method implementations
############################################################
run_dorfman_EN <- function(X, y) {
  dorfman_screening(X, y, method = "EN")
}

run_dorfman_Adaptive_EN <- function(X, y) {
  dorfman_screening(X, y, method = "Adaptive_EN")
}

run_dorfman_sparse_EN <- function(X, y) {
  # Force use of sparse correlation
  result <- dorfman_screening(X, y, method = "EN", tune_cv = FALSE)
  # Use the best rho from grid based on ARI
  return(result)
}

run_dorfman_sparse_Adaptive_EN <- function(X, y) {
  dorfman_screening(X, y, method = "Adaptive_EN", tune_cv = FALSE)
}

############################################################
## Other comparison methods (same as before)
############################################################
run_group_ar2 <- function(X, y) {
  # Implement Group_AR2_gp_LASSO
  cor_mat <- cor(X)
  dist_mat <- as.dist(1 - abs(cor_mat))
  hc <- hclust(dist_mat, method = "average")
  
  # Try different cut heights
  best_rmse <- Inf
  best_groups <- NULL
  
  for (h in seq(0.1, 0.9, by = 0.2)) {
    groups <- cutree(hc, h = h)
    
    # Calculate group scores
    unique_groups <- unique(groups)
    group_scores <- numeric(length(unique_groups))
    
    for (k in seq_along(unique_groups)) {
      group_idx <- which(groups == unique_groups[k])
      if (length(group_idx) > 1) {
        # Fit group model
        tryCatch({
          model <- lm(y ~ X[, group_idx])
          group_scores[k] <- summary(model)$adj.r.squared
        }, error = function(e) {
          group_scores[k] <- -Inf
        })
      }
    }
    
    # Select top groups
    top_groups <- unique_groups[order(group_scores, decreasing = TRUE)[1:min(10, length(unique_groups))]]
    selected_vars <- which(groups %in% top_groups)
    
    if (length(selected_vars) > 0) {
      # Fit LASSO on selected variables
      cv_fit <- cv.glmnet(X[, selected_vars], y, alpha = 1)
      final_selected <- which(coef(cv_fit, s = "lambda.min")[-1] != 0)
      final_vars <- selected_vars[final_selected]
      
      # Calculate RMSE via CV
      folds <- sample(rep(1:5, length.out = nrow(X)))
      rmse_vals <- numeric(5)
      for (fold in 1:5) {
        train_idx <- which(folds != fold)
        test_idx <- which(folds == fold)
        
        if (length(final_vars) > 0) {
          model <- lm(y[train_idx] ~ X[train_idx, final_vars])
          pred <- predict(model, data.frame(X[test_idx, final_vars]))
        } else {
          pred <- mean(y[train_idx])
        }
        rmse_vals[fold] <- sqrt(mean((y[test_idx] - pred)^2))
      }
      
      avg_rmse <- mean(rmse_vals)
      if (avg_rmse < best_rmse) {
        best_rmse <- avg_rmse
        best_groups <- final_vars
      }
    }
  }
  
  return(best_groups)
}

run_lasso <- function(X, y) {
  cv_fit <- cv.glmnet(X, y, alpha = 1)
  which(coef(cv_fit, s = "lambda.min")[-1] != 0)
}

run_en <- function(X, y) {
  cv_fit <- cv.glmnet(X, y, alpha = 0.5)
  which(coef(cv_fit, s = "lambda.min")[-1] != 0)
}

run_adaptive_lasso <- function(X, y) {
  # First stage: get initial estimate
  cv_init <- cv.glmnet(X, y, alpha = 1)
  beta_init <- as.numeric(coef(cv_init, s = "lambda.min"))[-1]
  
  # Adaptive weights
  w <- 1 / (abs(beta_init) + 1e-6)
  
  # Second stage: adaptive fit
  cv_adapt <- cv.glmnet(X, y, alpha = 1, penalty.factor = w)
  which(coef(cv_adapt, s = "lambda.min")[-1] != 0)
}

run_adaptive_en <- function(X, y) {
  # First stage: get initial estimate
  cv_init <- cv.glmnet(X, y, alpha = 0.5)
  beta_init <- as.numeric(coef(cv_init, s = "lambda.min"))[-1]
  
  # Adaptive weights
  w <- 1 / (abs(beta_init) + 1e-6)
  
  # Second stage: adaptive fit
  cv_adapt <- cv.glmnet(X, y, alpha = 0.5, penalty.factor = w)
  which(coef(cv_adapt, s = "lambda.min")[-1] != 0)
}

run_sis_lasso <- function(X, y) {
  sis_fit <- SIS(X, y, family = "gaussian", tune = "bic")
  sis_fit$sis.ix0
}

############################################################
## Metrics calculation (same as before)
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
  
  n_test <- length(y_test)
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
## Main simulation loop
############################################################
methods <- c(
  "Dorfman_EN",
  "Dorfman_Adaptive_EN",
  "Dorfman_sparse_EN",
  "Dorfman_sparse_Adaptive_EN",
  "Adaptive_LASSO",
  "Adaptive_EN",
  "LASSO",
  "EN",
  "SIS_LASSO",
  "Group_AR2_gp_LASSO"
)

results <- list()
timing_results <- list()

for (sim in 1:n_sim) {
  dat  <- simulated_data[[sim]]
  start_time <- Sys.time()
  
  X    <- dat$X
  y    <- dat$y
  true_vars <- which(dat$beta != 0)
  
  # Split data
  train_idx <- sample(nrow(X), 0.6 * nrow(X))
  X_train   <- X[train_idx, ]
  y_train   <- y[train_idx]
  X_test    <- X[-train_idx, ]
  y_test    <- y[-train_idx]
  
  # Run all methods
  selected_methods <- list(
    tryCatch(run_dorfman_EN(X_train, y_train)$final_selected, 
             error = function(e) integer(0)),
    tryCatch(run_dorfman_Adaptive_EN(X_train, y_train)$final_selected, 
             error = function(e) integer(0)),
    tryCatch(run_dorfman_sparse_EN(X_train, y_train)$final_selected, 
             error = function(e) integer(0)),
    tryCatch(run_dorfman_sparse_Adaptive_EN(X_train, y_train)$final_selected, 
             error = function(e) integer(0)),
    tryCatch(run_adaptive_lasso(X_train, y_train), 
             error = function(e) integer(0)),
    tryCatch(run_adaptive_en(X_train, y_train), 
             error = function(e) integer(0)),
    tryCatch(run_lasso(X_train, y_train), 
             error = function(e) integer(0)),
    tryCatch(run_en(X_train, y_train), 
             error = function(e) integer(0)),
    tryCatch(run_sis_lasso(X_train, y_train), 
             error = function(e) integer(0)),
    tryCatch(run_group_ar2(X_train, y_train), 
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
  
  cat(sprintf("Sim %2d completed\n", sim))
}

# Process results
results_df <- bind_rows(results)

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

print(performance_summary)

# Create visualization
ggplot(results_df, aes(x = Method, y = F1, fill = Method)) +
  geom_boxplot() +
  labs(title = "F1 Score Comparison of Methods",
       subtitle = "Normal Case (p=1000)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(results_df, aes(x = Method, y = RMSE, fill = Method)) +
  geom_boxplot() +
  labs(title = "RMSE Comparison of Methods",
       subtitle = "Normal Case (p=1000)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Summary table for paper
final_table <- performance_summary %>%
  mutate(
    Variables = sprintf("%.1f (%.2f)", Mean_Vars, SD_Vars),
    TPR = sprintf("%.3f (%.3f)", Mean_TPR, SD_TPR),
    FDR = sprintf("%.3f (%.3f)", Mean_FDR, SD_FDR),
    F1 = sprintf("%.3f (%.3f)", Mean_F1, SD_F1),
    RMSE = sprintf("%.3f (%.3f)", Mean_RMSE, SD_RMSE)
  ) %>%
  select(Method, Variables, TPR, FDR, F1, RMSE)

print(final_table)
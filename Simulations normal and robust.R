
p = 500 
Normal simulation 
```{r}
library(MASS)
set.seed(123)

# Parameters
n_obs <- 200          # Number of observations
n_groups <- 100       # Number of variable groups
n_per_group <- 5      # Variables per group
p <- n_groups * n_per_group  # Total variables (500)
n_simulations <- 100  # Number of simulations 

# Correlation structure (20 groups at each level)
cor_levels <- rep(c(0.1, 0.3, 0.5, 0.7, 0.9), each = 20)

# Active groups (ALL variables in these groups are active)
active_groups <- data.frame(
  group_id = c(1, 21, 41, 61, 81),  # ONE group per correlation level
  cor_level = c(0.1, 0.3, 0.5, 0.7, 0.9),
  beta = c(1, 0.9, 0.7, 0.5, 0.3)
)

# Initialize list to store all simulations
all_simulations <- vector("list", n_simulations)

# Function to generate one simulation
generate_simulation <- function() {
  # Initialize design matrix X (200 x 500)
  X <- matrix(rnorm(n_obs * p), nrow = n_obs, ncol = p)
  
  # Initialize beta vector (all zeros initially)
  beta <- numeric(p)
  
  # Create response variable y
  y <- numeric(n_obs)
  
  # Generate data with perfect within-group correlation
  for (g in 1:n_groups) {
    group_cols <- ((g-1)*n_per_group + 1):(g*n_per_group)
    rho <- cor_levels[g]
    
    # Create covariance matrix
    Sigma <- matrix(rho, n_per_group, n_per_group)
    diag(Sigma) <- 1
    
    # Generate correlated variables
    L <- chol(Sigma)
    Z <- matrix(rnorm(n_obs * n_per_group), n_obs)
    X[, group_cols] <- Z %*% L
    
    # Set ALL variables as active if this is an active group
    if (g %in% active_groups$group_id) {
      effect_size <- active_groups$beta[active_groups$group_id == g]
      # ALL 5 variables in this group are active with the same beta
      beta[group_cols] <- effect_size
      y <- y + effect_size * rowSums(X[, group_cols])  # Sum of all 5 active vars
    }
  }
  
  # Add noise to response
  y <- y + rnorm(n_obs)
  
  # Create variable group info
  var_info <- data.frame(
    var_id = 1:p,
    group_id = rep(1:n_groups, each = n_per_group),
    cor_level = rep(cor_levels, each = n_per_group),
    is_active = rep(1:n_groups %in% active_groups$group_id, each = n_per_group),
    is_signal_var = (1:p) %in% unlist(lapply(active_groups$group_id, function(g) {
      ((g-1)*n_per_group + 1):(g*n_per_group)
    }))
  )
  
  return(list(
    X = X,
    y = y,
    var_info = var_info,
    beta = beta,
    active_groups = active_groups
  ))
}

# Generate all simulations
for (i in 1:n_simulations) {
  all_simulations[[i]] <- generate_simulation()
  cat("Generated simulation", i, "of", n_simulations, "\n")
}

# Save all simulations
saveRDS(all_simulations, "100_simulations_noise1.rds") 

# Verification function
verify_group_correlations <- function(X, var_info) {
  results <- data.frame()
  for (g in unique(var_info$group_id)) {
    group_vars <- var_info$var_id[var_info$group_id == g]
    cor_mat <- cor(X[, group_vars])
    avg_cor <- mean(cor_mat[upper.tri(cor_mat)])
    
    results <- rbind(results, data.frame(
      group_id = g,
      specified_cor = var_info$cor_level[var_info$group_id == g][1],
      empirical_cor = avg_cor,
      is_active = var_info$is_active[var_info$group_id == g][1]
    ))
  }
  return(results)
}

# Verification on first simulation
cat("\nVerification for first simulation:\n")
first_sim <- all_simulations[[1]]
cor_results <- verify_group_correlations(first_sim$X, first_sim$var_info)

library(dplyr)
cor_results %>%
  group_by(specified_cor) %>%
  summarise(mean_empirical = mean(empirical_cor),
            sd_empirical = sd(empirical_cor)) %>%
  print()

# Check active variables
cat("\nActive Variables (25 total):\n")
active_vars <- which(first_sim$beta != 0)
active_vars_info <- data.frame(
  var_id = active_vars,
  group_id = first_sim$var_info$group_id[active_vars],
  cor_level = first_sim$var_info$cor_level[active_vars],
  beta = first_sim$beta[active_vars]
)
print(active_vars_info)

# Count active variables by correlation level
cat("\nActive Variables by Correlation Level:\n")
active_vars_info %>%
  group_by(cor_level, beta) %>%
  summarise(count = n()) %>%
  print()
```


```{r}
library(MASS)
library(moments)
library(Matrix)
library(dplyr)

set.seed(123)

# Parameters
n_obs        <- 200
n_groups     <- 100
n_per_group  <- 5
p            <- n_groups * n_per_group
n_simulations <- 100

# Correlation structure (20 groups at each level)
cor_levels <- rep(c(0.1, 0.3, 0.5, 0.7, 0.9), each = 20)

# Active groups (ALL variables in these groups are active)
active_groups <- data.frame(
  group_id   = c(1, 21, 41, 61, 81),
  cor_level  = c(0.1, 0.3, 0.5, 0.7, 0.9),
  beta       = c(2.0, 1.8, 1.5, 1.2, 1.0)  # stronger signal than before
)

inactive_groups <- setdiff(1:n_groups, active_groups$group_id)

all_simulations_corrupted <- vector("list", n_simulations)

generate_corrupted_simulation <- function() {
  X_clean <- matrix(NA_real_, nrow = n_obs, ncol = p)
  beta    <- numeric(p)
  y_clean <- numeric(n_obs)
  
  # Clean X and signal
  for (g in 1:n_groups) {
    group_cols <- ((g - 1) * n_per_group + 1):(g * n_per_group)
    rho        <- cor_levels[g]
    Sigma      <- matrix(rho, n_per_group, n_per_group)
    diag(Sigma) <- 1
    L <- chol(Sigma)
    Z <- matrix(rnorm(n_obs * n_per_group), n_obs)
    X_clean[, group_cols] <- Z %*% L
    
    if (g %in% active_groups$group_id) {
      eff <- active_groups$beta[active_groups$group_id == g]
      beta[group_cols] <- eff
      y_clean <- y_clean + eff * rowSums(X_clean[, group_cols, drop = FALSE])
    }
  }
  
  y_clean <- y_clean + rnorm(n_obs, 0, 0.3)
  
  X <- X_clean
  
  # Mild global noise
  X <- X + matrix(rnorm(n_obs * p, 0, 0.1), n_obs, p)
  
  # Aligned X+y outliers in INACTIVE variables (small fraction of rows)
  n_aligned <- 8
  aligned_idx <- sample(1:n_obs, size = n_aligned, replace = FALSE)
  vars_inactive <- unlist(lapply(inactive_groups, function(g) {
    ((g - 1) * n_per_group + 1):(g * n_per_group)
  }))
  
  for (i in aligned_idx) {
    bad_vars <- sample(vars_inactive, size = 40, replace = FALSE)
    base     <- rnorm(1, mean = 8, sd = 1)
    X[i, bad_vars] <- base + rnorm(length(bad_vars), 0, 0.3)
  }
  
  # False cross-group correlations in inactive vars (moderate)
  n_pairs <- 20
  false_pair_idx <- sample(setdiff(1:n_obs, aligned_idx),
                           size = 20, replace = FALSE)
  for (i in false_pair_idx) {
    vset <- sample(vars_inactive, size = 2 * n_pairs, replace = FALSE)
    for (k in seq(1, 2 * n_pairs - 1, by = 2)) {
      v1 <- vset[k]
      v2 <- vset[k + 1]
      X[i, v2] <- X[i, v1] + rnorm(1, 0, 0.1)
    }
  }
  
  # Moderate leverage mostly in inactive variables
  leverage_idx <- sample(setdiff(1:n_obs, c(aligned_idx, false_pair_idx)),
                         size = 12, replace = FALSE)
  for (i in leverage_idx) {
    lev_inactive <- sample(vars_inactive, size = 35, replace = FALSE)
    rem          <- setdiff(1:p, lev_inactive)
    lev_any      <- sample(rem, size = 10, replace = FALSE)
    lev_vars     <- c(lev_inactive, lev_any)
    X[i, lev_vars] <- X[i, lev_vars] * 5 + rnorm(length(lev_vars), 0, 2)
  }
  
  # Batch effects – large only in inactive groups, mild in active
  batch1_idx <- sample(1:n_obs, round(0.15 * n_obs))
  batch2_idx <- sample(setdiff(1:n_obs, batch1_idx), round(0.15 * n_obs))
  
  batch1_inactive <- sample(inactive_groups, 12)
  for (g in batch1_inactive) {
    group_cols <- ((g - 1) * n_per_group + 1):(g * n_per_group)
    X[batch1_idx, group_cols] <- X[batch1_idx, group_cols] + 3
  }
  
  batch2_inactive <- sample(setdiff(inactive_groups, batch1_inactive), 12)
  for (g in batch2_inactive) {
    group_cols <- ((g - 1) * n_per_group + 1):(g * n_per_group)
    X[batch2_idx, group_cols] <- X[batch2_idx, group_cols] - 3
  }
  
  # Very mild symmetric shifts in active groups (does not hurt non-robust much)
  for (g in active_groups$group_id) {
    group_cols <- ((g - 1) * n_per_group + 1):(g * n_per_group)
    X[batch1_idx, group_cols] <- X[batch1_idx, group_cols] + rnorm(length(batch1_idx), 0, 0.3)
    X[batch2_idx, group_cols] <- X[batch2_idx, group_cols] + rnorm(length(batch2_idx), 0, 0.3)
  }
  
  # Response: strong linear signal + moderate heavy-tailed + hetero + aligned spikes
  y_signal <- as.numeric(X_clean %*% beta)
  
  eps_heavy <- rt(n_obs, df = 5) * 0.8
  
  g0      <- active_groups$group_id[1]
  g0_cols <- ((g0 - 1) * n_per_group + 1):(g0 * n_per_group)
  x_ref   <- rowMeans(X_clean[, g0_cols, drop = FALSE])
  hetero_scale <- 0.7 + 1.0 * abs(x_ref)
  eps_hetero   <- rnorm(n_obs, 0, hetero_scale)
  
  eps_spike <- numeric(n_obs)
  eps_spike[aligned_idx] <- rnorm(length(aligned_idx), mean = 8, sd = 1.5)
  
  epsilon <- eps_heavy + eps_hetero + eps_spike
  
  y <- y_signal + epsilon
  
  var_info <- data.frame(
    var_id        = 1:p,
    group_id      = rep(1:n_groups, each = n_per_group),
    cor_level     = rep(cor_levels, each = n_per_group),
    is_active     = rep(1:n_groups %in% active_groups$group_id, each = n_per_group),
    is_signal_var = (1:p) %in% unlist(lapply(active_groups$group_id, function(g) {
      ((g - 1) * n_per_group + 1):(g * n_per_group)
    }))
  )
  
  list(
    X              = X,
    X_clean        = X_clean,
    y              = y,
    y_clean        = y_clean,
    var_info       = var_info,
    beta           = beta,
    active_groups  = active_groups,
    corruption_info = list(
      aligned_obs      = aligned_idx,
      false_pair_obs   = false_pair_idx,
      leverage_obs     = leverage_idx,
      batch1_obs       = batch1_idx,
      batch2_obs       = batch2_idx,
      corruption_types = c(
        "mild global noise",
        "aligned X+y outliers in inactive vars",
        "false cross-group correlations in inactive vars",
        "moderate leverage mostly in inactive vars",
        "batch shifts mainly in inactive groups",
        "heavy-tailed + heteroscedastic + moderate spikes in y"
      )
    )
  )
}

for (i in 1:n_simulations) {
  all_simulations_corrupted[[i]] <- generate_corrupted_simulation()
  if (i %% 10 == 0) cat("Generated robust-favoring simulation", i, "of", n_simulations, "\n")
}

saveRDS(all_simulations_corrupted, "100_corrupted_simulations_varying_cor_robust.rds")

cat("\nFinished generating 100 robust-favoring simulations and saved to:\n")
cat("  100_corrupted_simulations_varying_cor_robust.rds\n")
```


p = 1000 
"100_simulations_p1000_noise1.rds" 
```{r}
library(MASS)
set.seed(123)

# Parameters
n_obs <- 200          # Number of observations
n_groups <- 200       # Number of variable groups (increased from 100)
n_per_group <- 5      # Variables per group
p <- n_groups * n_per_group  # Total variables (1000)
n_simulations <- 100  # Number of simulations 

# Correlation structure (40 groups at each level)
cor_levels <- rep(c(0.1, 0.3, 0.5, 0.7, 0.9), each = 40)

# Active groups (ALL variables in these groups are active)
# One group per correlation level
active_groups <- data.frame(
  group_id = c(1, 41, 81, 121, 161),  # ONE group per correlation level
  cor_level = c(0.1, 0.3, 0.5, 0.7, 0.9),
  beta = c(1, 0.9, 0.7, 0.5, 0.3)
)

# Initialize list to store all simulations
all_simulations <- vector("list", n_simulations)

# Function to generate one simulation
generate_simulation <- function() {
  # Initialize design matrix X (200 x 1000)
  X <- matrix(rnorm(n_obs * p), nrow = n_obs, ncol = p)
  
  # Initialize beta vector (all zeros initially)
  beta <- numeric(p)
  
  # Create response variable y
  y <- numeric(n_obs)
  
  # Generate data with perfect within-group correlation
  for (g in 1:n_groups) {
    group_cols <- ((g-1)*n_per_group + 1):(g*n_per_group)
    rho <- cor_levels[g]
    
    # Create covariance matrix
    Sigma <- matrix(rho, n_per_group, n_per_group)
    diag(Sigma) <- 1
    
    # Generate correlated variables
    L <- chol(Sigma)
    Z <- matrix(rnorm(n_obs * n_per_group), n_obs)
    X[, group_cols] <- Z %*% L
    
    # Set ALL variables as active if this is an active group
    if (g %in% active_groups$group_id) {
      effect_size <- active_groups$beta[active_groups$group_id == g]
      # ALL 5 variables in this group are active with the same beta
      beta[group_cols] <- effect_size
      y <- y + effect_size * rowSums(X[, group_cols])  # Sum of all 5 active vars
    }
  }
  
  # Add noise to response
  y <- y + rnorm(n_obs)
  
  # Create variable group info
  var_info <- data.frame(
    var_id = 1:p,
    group_id = rep(1:n_groups, each = n_per_group),
    cor_level = rep(cor_levels, each = n_per_group),
    is_active = rep(1:n_groups %in% active_groups$group_id, each = n_per_group),
    is_signal_var = (1:p) %in% unlist(lapply(active_groups$group_id, function(g) {
      ((g-1)*n_per_group + 1):(g*n_per_group)
    }))
  )
  
  return(list(
    X = X,
    y = y,
    var_info = var_info,
    beta = beta,
    active_groups = active_groups
  ))
}

# Generate all simulations
for (i in 1:n_simulations) {
  all_simulations[[i]] <- generate_simulation()
  cat("Generated simulation", i, "of", n_simulations, "\n")
}

# Save all simulations
saveRDS(all_simulations, "100_simulations_p1000_noise1.rds") 

# Verification function
verify_group_correlations <- function(X, var_info) {
  results <- data.frame()
  for (g in unique(var_info$group_id)) {
    group_vars <- var_info$var_id[var_info$group_id == g]
    cor_mat <- cor(X[, group_vars])
    avg_cor <- mean(cor_mat[upper.tri(cor_mat)])
    
    results <- rbind(results, data.frame(
      group_id = g,
      specified_cor = var_info$cor_level[var_info$group_id == g][1],
      empirical_cor = avg_cor,
      is_active = var_info$is_active[var_info$group_id == g][1]
    ))
  }
  return(results)
}

# Verification on first simulation
cat("\nVerification for first simulation:\n")
first_sim <- all_simulations[[1]]
cor_results <- verify_group_correlations(first_sim$X, first_sim$var_info)

library(dplyr)
cor_results %>%
  group_by(specified_cor) %>%
  summarise(mean_empirical = mean(empirical_cor),
            sd_empirical = sd(empirical_cor),
            n_groups = n()) %>%
  print()

# Check active variables
cat("\nActive Variables (25 total):\n")
active_vars <- which(first_sim$beta != 0)
active_vars_info <- data.frame(
  var_id = active_vars,
  group_id = first_sim$var_info$group_id[active_vars],
  cor_level = first_sim$var_info$cor_level[active_vars],
  beta = first_sim$beta[active_vars]
)
print(active_vars_info)

# Count active variables by correlation level
cat("\nActive Variables by Correlation Level:\n") 
active_vars_info %>%
  group_by(cor_level, beta) %>%
  summarise(count = n()) %>%
  print()

# Summary statistics
cat("\nSummary Statistics:\n")
cat("Total variables:", p, "\n")
cat("Total groups:", n_groups, "\n")
cat("Groups per correlation level:", 40, "\n")
cat("Active groups:", nrow(active_groups), "\n")
cat("Active variables:", length(active_vars), "\n")
cat("Signal-to-noise ratio variables:", 25, "/", 1000, "=", 25/1000*100, "%\n")

# Check group structure
cat("\nActive Groups Details:\n")
for (i in 1:nrow(active_groups)) {
  g <- active_groups$group_id[i]
  group_vars <- ((g-1)*n_per_group + 1):(g*n_per_group)
  cat(sprintf("Group %d (correlation=%.1f, beta=%.1f): Variables %d-%d\n",
              g, active_groups$cor_level[i], active_groups$beta[i],
              min(group_vars), max(group_vars)))
}
```

corrupted p=1000 
100_corrupted_simulations_p1000_varying_cor_robust.rds 
```{r}
library(MASS)
library(moments)
library(Matrix)
library(dplyr)

set.seed(123)

# Parameters - Modified for p=1000
n_obs        <- 200
n_groups     <- 200       # Increased from 100
n_per_group  <- 5
p            <- n_groups * n_per_group  # 1000
n_simulations <- 100

# Correlation structure (40 groups at each level)
cor_levels <- rep(c(0.1, 0.3, 0.5, 0.7, 0.9), each = 40)

# Active groups (ALL variables in these groups are active)
active_groups <- data.frame(
  group_id   = c(1, 41, 81, 121, 161),  # Adjusted for 40 groups per level
  cor_level  = c(0.1, 0.3, 0.5, 0.7, 0.9),
  beta       = c(2.0, 1.8, 1.5, 1.2, 1.0)  # stronger signal than before
)

inactive_groups <- setdiff(1:n_groups, active_groups$group_id)

all_simulations_corrupted <- vector("list", n_simulations)

generate_corrupted_simulation <- function() {
  X_clean <- matrix(NA_real_, nrow = n_obs, ncol = p)
  beta    <- numeric(p)
  y_clean <- numeric(n_obs)
  
  # Clean X and signal
  for (g in 1:n_groups) {
    group_cols <- ((g - 1) * n_per_group + 1):(g * n_per_group)
    rho        <- cor_levels[g]
    Sigma      <- matrix(rho, n_per_group, n_per_group)
    diag(Sigma) <- 1
    L <- chol(Sigma)
    Z <- matrix(rnorm(n_obs * n_per_group), n_obs)
    X_clean[, group_cols] <- Z %*% L
    
    if (g %in% active_groups$group_id) {
      eff <- active_groups$beta[active_groups$group_id == g]
      beta[group_cols] <- eff
      y_clean <- y_clean + eff * rowSums(X_clean[, group_cols, drop = FALSE])
    }
  }
  
  y_clean <- y_clean + rnorm(n_obs, 0, 0.3)
  
  X <- X_clean 
  
  # Mild global noise
  X <- X + matrix(rnorm(n_obs * p, 0, 0.1), n_obs, p)
  
  # Aligned X+y outliers in INACTIVE variables (small fraction of rows)
  n_aligned <- 8
  aligned_idx <- sample(1:n_obs, size = n_aligned, replace = FALSE)
  vars_inactive <- unlist(lapply(inactive_groups, function(g) {
    ((g - 1) * n_per_group + 1):(g * n_per_group) 
  }))
  
  for (i in aligned_idx) {
    bad_vars <- sample(vars_inactive, size = 80, replace = FALSE)  # Increased from 40
    base     <- rnorm(1, mean = 8, sd = 1)
    X[i, bad_vars] <- base + rnorm(length(bad_vars), 0, 0.3)
  }
  
  # False cross-group correlations in inactive vars (moderate)
  n_pairs <- 40  # Increased from 20
  false_pair_idx <- sample(setdiff(1:n_obs, aligned_idx),
                           size = 20, replace = FALSE)
  for (i in false_pair_idx) {
    vset <- sample(vars_inactive, size = 2 * n_pairs, replace = FALSE)
    for (k in seq(1, 2 * n_pairs - 1, by = 2)) {
      v1 <- vset[k]
      v2 <- vset[k + 1]
      X[i, v2] <- X[i, v1] + rnorm(1, 0, 0.1)
    }
  }
  
  # Moderate leverage mostly in inactive variables
  leverage_idx <- sample(setdiff(1:n_obs, c(aligned_idx, false_pair_idx)),
                         size = 12, replace = FALSE)
  for (i in leverage_idx) {
    lev_inactive <- sample(vars_inactive, size = 70, replace = FALSE)  # Increased from 35
    rem          <- setdiff(1:p, lev_inactive)
    lev_any      <- sample(rem, size = 20, replace = FALSE)  # Increased from 10
    lev_vars     <- c(lev_inactive, lev_any)
    X[i, lev_vars] <- X[i, lev_vars] * 5 + rnorm(length(lev_vars), 0, 2)
  }
  
  # Batch effects – large only in inactive groups, mild in active
  batch1_idx <- sample(1:n_obs, round(0.15 * n_obs))
  batch2_idx <- sample(setdiff(1:n_obs, batch1_idx), round(0.15 * n_obs))
  
  # Increase number of inactive groups affected by batch effects
  batch1_inactive <- sample(inactive_groups, 24)  # Increased from 12
  for (g in batch1_inactive) {
    group_cols <- ((g - 1) * n_per_group + 1):(g * n_per_group)
    X[batch1_idx, group_cols] <- X[batch1_idx, group_cols] + 3
  }
  
  batch2_inactive <- sample(setdiff(inactive_groups, batch1_inactive), 24)  # Increased from 12
  for (g in batch2_inactive) {
    group_cols <- ((g - 1) * n_per_group + 1):(g * n_per_group)
    X[batch2_idx, group_cols] <- X[batch2_idx, group_cols] - 3
  }
  
  # Very mild symmetric shifts in active groups (does not hurt non-robust much)
  for (g in active_groups$group_id) {
    group_cols <- ((g - 1) * n_per_group + 1):(g * n_per_group)
    X[batch1_idx, group_cols] <- X[batch1_idx, group_cols] + rnorm(length(batch1_idx), 0, 0.3)
    X[batch2_idx, group_cols] <- X[batch2_idx, group_cols] + rnorm(length(batch2_idx), 0, 0.3)
  }
  
  # Response: strong linear signal + moderate heavy-tailed + hetero + aligned spikes
  y_signal <- as.numeric(X_clean %*% beta)
  
  eps_heavy <- rt(n_obs, df = 5) * 0.8
  
  g0      <- active_groups$group_id[1]
  g0_cols <- ((g0 - 1) * n_per_group + 1):(g0 * n_per_group)
  x_ref   <- rowMeans(X_clean[, g0_cols, drop = FALSE])
  hetero_scale <- 0.7 + 1.0 * abs(x_ref)
  eps_hetero   <- rnorm(n_obs, 0, hetero_scale)
  
  eps_spike <- numeric(n_obs)
  eps_spike[aligned_idx] <- rnorm(length(aligned_idx), mean = 8, sd = 1.5)
  
  epsilon <- eps_heavy + eps_hetero + eps_spike
  
  y <- y_signal + epsilon
  
  var_info <- data.frame(
    var_id        = 1:p,
    group_id      = rep(1:n_groups, each = n_per_group),
    cor_level     = rep(cor_levels, each = n_per_group),
    is_active     = rep(1:n_groups %in% active_groups$group_id, each = n_per_group),
    is_signal_var = (1:p) %in% unlist(lapply(active_groups$group_id, function(g) {
      ((g - 1) * n_per_group + 1):(g * n_per_group)
    }))
  )
  
  list(
    X              = X,
    X_clean        = X_clean,
    y              = y,
    y_clean        = y_clean,
    var_info       = var_info,
    beta           = beta,
    active_groups  = active_groups,
    corruption_info = list(
      aligned_obs      = aligned_idx,
      false_pair_obs   = false_pair_idx,
      leverage_obs     = leverage_idx,
      batch1_obs       = batch1_idx,
      batch2_obs       = batch2_idx,
      corruption_types = c(
        "mild global noise",
        "aligned X+y outliers in inactive vars",
        "false cross-group correlations in inactive vars",
        "moderate leverage mostly in inactive vars",
        "batch shifts mainly in inactive groups",
        "heavy-tailed + heteroscedastic + moderate spikes in y"
      )
    )
  )
}

for (i in 1:n_simulations) {
  all_simulations_corrupted[[i]] <- generate_corrupted_simulation()
  if (i %% 10 == 0) cat("Generated robust-favoring simulation", i, "of", n_simulations, "\n")
}

saveRDS(all_simulations_corrupted, "100_corrupted_simulations_p1000_varying_cor_robust.rds")

# Verification
cat("\n=== VERIFICATION FOR FIRST SIMULATION (p=1000) ===\n")
first_sim <- all_simulations_corrupted[[1]]

cat("\n1. Dimensionality:\n")
cat("   n_obs:", n_obs, "\n")
cat("   p:", p, "\n")
cat("   n_groups:", n_groups, "\n")
cat("   Variables per group:", n_per_group, "\n")
cat("   Total active variables:", sum(first_sim$beta != 0), "\n")

cat("\n2. Active Groups:\n")
for (i in 1:nrow(first_sim$active_groups)) {
  g <- first_sim$active_groups$group_id[i]
  cat(sprintf("   Group %d: correlation=%.1f, beta=%.1f\n", 
              g, first_sim$active_groups$cor_level[i], first_sim$active_groups$beta[i]))
}

cat("\n3. Corruption Information:\n")
cat("   Types of corruption:\n")
for (type in first_sim$corruption_info$corruption_types) {
  cat("   -", type, "\n")
}
cat("   Number of aligned outliers:", length(first_sim$corruption_info$aligned_obs), "\n")
cat("   Number of false-pair observations:", length(first_sim$corruption_info$false_pair_obs), "\n")
cat("   Number of leverage points:", length(first_sim$corruption_info$leverage_obs), "\n")
cat("   Total corrupted observations (approx):", 
    length(unique(c(first_sim$corruption_info$aligned_obs,
                    first_sim$corruption_info$false_pair_obs,
                    first_sim$corruption_info$leverage_obs,
                    first_sim$corruption_info$batch1_obs,
                    first_sim$corruption_info$batch2_obs))), "\n")

cat("\n4. Signal vs Noise:\n")
cat("   Signal variables:", sum(first_sim$var_info$is_signal_var), "/", p, 
    sprintf("(%.1f%%)", sum(first_sim$var_info$is_signal_var)/p*100), "\n")
cat("   Inactive variables:", sum(!first_sim$var_info$is_signal_var), "/", p,
    sprintf("(%.1f%%)", sum(!first_sim$var_info$is_signal_var)/p*100), "\n")

cat("\n5. Response Distribution:\n")
cat("   Mean(y):", mean(first_sim$y), "\n")
cat("   SD(y):", sd(first_sim$y), "\n")
cat("   Skewness(y):", skewness(first_sim$y), "\n")
cat("   Kurtosis(y):", kurtosis(first_sim$y), "\n")

# Check within-group correlations
cat("\n6. Sample Within-Group Correlations (first 5 groups):\n")
for (g in 1:5) {
  group_cols <- ((g-1)*n_per_group + 1):(g*n_per_group)
  cor_mat <- cor(first_sim$X[, group_cols])
  avg_cor <- mean(cor_mat[upper.tri(cor_mat)])
  cat(sprintf("   Group %d (corr=%.1f): empirical avg = %.3f\n", 
              g, cor_levels[g], avg_cor))
}

cat("\nFinished generating 100 robust-favoring simulations with p=1000 and saved to:\n")
cat("  100_corrupted_simulations_p1000_varying_cor_robust.rds\n")
```




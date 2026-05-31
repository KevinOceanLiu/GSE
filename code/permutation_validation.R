# =========================================================
# 05_permutation_validation.R
# Permutation validation for GSE
# =========================================================

# ---- 1. Load packages ----
library(FNN)
library(dplyr)

# ---- 2. Set paths ----
residual_dir <- "result"
output_dir <- "result"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ---- 3. User parameters ----
model_list <- c(
  "lm", "rf", "xgbTree", "gbm", "knn",
  "gam", "gwr", "svmRadial", "cubist"
)

# Model-specific optimal K values selected from GSE-K curves
optimal_k <- c(
  lm = 19,
  rf = 25,
  xgbTree = 23,
  gbm = 23,
  knn = 25,
  gam = 21,
  gwr = 27,
  svmRadial = 23,
  cubist = 25
)

n_perm <- 999
set.seed(123)

# ---- 4. Standardize residual column name ----
standardize_residual_column <- function(df) {
  
  if ("residuals" %in% names(df)) {
    df$residual <- df$residuals
  }
  
  if (!"residual" %in% names(df)) {
    stop("Residual column not found. Expected 'residual' or 'residuals'.")
  }
  
  return(df)
}

# ---- 5. Calculate k-NN GSE ----
calculate_gse_knn <- function(residual_data, k) {
  
  if (k < 2) {
    stop("k must be >= 2.")
  }
  
  if (nrow(residual_data) < k) {
    stop("Not enough data points for this k.")
  }
  
  required_cols <- c("x_proj", "y_proj", "residual")
  missing_cols <- setdiff(required_cols, names(residual_data))
  
  if (length(missing_cols) > 0) {
    stop(paste(
      "Missing columns for GSE calculation:",
      paste(missing_cols, collapse = ", ")
    ))
  }
  
  coords <- residual_data[, c("x_proj", "y_proj")]
  nn_indices <- get.knnx(coords, coords, k = k)$nn.index
  
  n <- nrow(residual_data)
  local_rmse <- numeric(n)
  local_variance <- numeric(n)
  
  for (i in seq_len(n)) {
    
    local_indices <- nn_indices[i, ]
    local_errors <- residual_data$residual[local_indices]
    
    local_rmse[i] <- sqrt(mean(local_errors^2))
    local_variance[i] <- var(local_errors) * ((k - 1) / k)
  }
  
  variance_sum <- sum(local_variance, na.rm = TRUE)
  
  if (variance_sum == 0) {
    weights <- rep(1 / n, n)
  } else {
    weights <- local_variance / variance_sum
  }
  
  gse <- sqrt(sum(weights * local_rmse^2, na.rm = TRUE))
  
  return(gse)
}

# ---- 6. Run permutation validation ----
summary_list <- list()

for (model_name in model_list) {
  
  message("Permutation validation for model: ", model_name)
  
  residual_file <- file.path(
    residual_dir,
    paste0(model_name, "_test_residuals.csv")
  )
  
  if (!file.exists(residual_file)) {
    warning(paste("Residual file not found:", residual_file))
    next
  }
  
  residual_data <- read.csv(residual_file, fileEncoding = "UTF-8")
  residual_data <- standardize_residual_column(residual_data)
  
  required_cols <- c("ID", "longitude", "latitude", "x_proj", "y_proj", "residual")
  missing_cols <- setdiff(required_cols, names(residual_data))
  
  if (length(missing_cols) > 0) {
    stop(paste(
      "Missing columns in", model_name, ":",
      paste(missing_cols, collapse = ", ")
    ))
  }
  
  residual_data <- residual_data[order(residual_data$ID), ]
  
  k_use <- optimal_k[[model_name]]
  
  if (is.null(k_use) || is.na(k_use)) {
    stop(paste("Optimal K is missing for model:", model_name))
  }
  
  # Observed RMSE and GSE
  rmse_obs <- sqrt(mean(residual_data$residual^2))
  gse_obs <- calculate_gse_knn(
    residual_data = residual_data,
    k = k_use
  )
  
  # Permutation test:
  # coordinates are fixed, residuals are randomly reassigned
  perm_gse <- numeric(n_perm)
  
  for (b in seq_len(n_perm)) {
    
    perm_data <- residual_data
    perm_data$residual <- sample(residual_data$residual, replace = FALSE)
    
    perm_gse[b] <- calculate_gse_knn(
      residual_data = perm_data,
      k = k_use
    )
  }
  
  mean_perm_gse <- mean(perm_gse)
  sd_perm_gse <- sd(perm_gse)
  delta_gse_percent <- ((gse_obs - mean_perm_gse) / mean_perm_gse) * 100
  
  # One-sided p-value: probability that permuted GSE >= observed GSE
  p_value <- (sum(perm_gse >= gse_obs) + 1) / (n_perm + 1)
  
  summary_list[[model_name]] <- data.frame(
    Model = model_name,
    Optimal_K = k_use,
    RMSE = rmse_obs,
    Observed_GSE = gse_obs,
    Mean_permuted_GSE = mean_perm_gse,
    SD_permuted_GSE = sd_perm_gse,
    Delta_GSE_percent = delta_gse_percent,
    Permutation_p_value = p_value
  )
  
  message("  Optimal K: ", k_use)
  message("  RMSE: ", round(rmse_obs, 4))
  message("  Observed GSE: ", round(gse_obs, 4))
  message("  Mean permuted GSE: ", round(mean_perm_gse, 4))
  message("  Delta GSE (%): ", round(delta_gse_percent, 3))
  message("  Permutation p-value: ", round(p_value, 4))
}

# ---- 7. Export summary ----
if (length(summary_list) == 0) {
  stop("No permutation results were generated. Please check residual files.")
}

permutation_summary <- bind_rows(summary_list) %>%
  mutate(
    RMSE = round(RMSE, 4),
    Observed_GSE = round(Observed_GSE, 4),
    Mean_permuted_GSE = round(Mean_permuted_GSE, 4),
    SD_permuted_GSE = round(SD_permuted_GSE, 4),
    Delta_GSE_percent = round(Delta_GSE_percent, 3),
    Permutation_p_value = round(Permutation_p_value, 4)
  )

output_file <- file.path(output_dir, "permutation_results.csv")

write.csv(
  permutation_summary,
  file = output_file,
  row.names = FALSE
)

print(permutation_summary)

message("Permutation validation results saved to: ", output_file)
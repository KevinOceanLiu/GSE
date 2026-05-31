# =========================================================
# 09_runtime_comparison.R
# Runtime comparison among RMSE, fixed-K GSE, and GSEopt
# =========================================================

# ---- 1. Load packages ----
library(FNN)
library(microbenchmark)
library(dplyr)

# ---- 2. Set paths ----
residual_file <- "result/rf_test_residuals.csv"
output_dir <- "result"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (!file.exists(residual_file)) {
  stop(paste("Residual file not found:", residual_file))
}

# ---- 3. Read residual data ----
residual_data <- read.csv(
  residual_file,
  fileEncoding = "UTF-8",
  check.names = FALSE
)

# Support both old and new residual column names
if ("residual" %in% names(residual_data) && !"residuals" %in% names(residual_data)) {
  residual_data$residuals <- residual_data$residual
}

required_cols <- c("x_proj", "y_proj", "residuals")
missing_cols <- setdiff(required_cols, names(residual_data))

if (length(missing_cols) > 0) {
  stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
}

residual_data <- residual_data %>%
  select(x_proj, y_proj, residuals) %>%
  na.omit()

# ---- 4. Candidate K values ----
k_values <- c(
  3, 5, 7, 9, 11, 13, 15, 17, 19, 21,
  23, 25, 27, 29, 31, 33, 35, 37, 39, 41,
  43, 45, 47, 49, 51, 53, 55, 57, 59, 61,
  63, 65, 67, 69, 71, 73, 75, 77, 79, 81,
  83, 85, 87, 89, 91, 93, 95, 97, 99, 101
)

fixed_k <- 25

# ---- 5. RMSE ----
calculate_rmse <- function(residual_data) {
  sqrt(mean(residual_data$residuals^2))
}

# ---- 6. Fixed-K GSE ----
calculate_gse_fixed_k <- function(residual_data, k) {
  
  if (k < 2) {
    stop("k must be >= 2.")
  }
  
  if (nrow(residual_data) < k) {
    stop("Not enough data points for this k.")
  }
  
  coords <- residual_data[, c("x_proj", "y_proj")]
  nn_indices <- get.knnx(coords, coords, k = k)$nn.index
  
  n <- nrow(residual_data)
  local_rmse <- numeric(n)
  local_variance <- numeric(n)
  
  for (i in seq_len(n)) {
    
    local_indices <- nn_indices[i, ]
    local_errors <- residual_data$residuals[local_indices]
    
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

# ---- 7. GSEopt ----
calculate_gseopt <- function(residual_data, k_values) {
  
  coords <- residual_data[, c("x_proj", "y_proj")]
  errors <- residual_data$residuals
  n <- nrow(residual_data)
  
  max_k <- max(k_values)
  
  if (n < max_k) {
    stop("Not enough data points for max K.")
  }
  
  # Search k-NN once and truncate the neighbor matrix for smaller K values
  nn_indices_max <- get.knnx(coords, coords, k = max_k)$nn.index
  
  gse_values <- numeric(length(k_values))
  
  for (j in seq_along(k_values)) {
    
    k <- k_values[j]
    nn_indices <- nn_indices_max[, 1:k, drop = FALSE]
    
    local_rmse <- numeric(n)
    local_variance <- numeric(n)
    
    for (i in seq_len(n)) {
      
      local_errors <- errors[nn_indices[i, ]]
      
      local_rmse[i] <- sqrt(mean(local_errors^2))
      local_variance[i] <- var(local_errors) * ((k - 1) / k)
    }
    
    variance_sum <- sum(local_variance, na.rm = TRUE)
    
    if (variance_sum == 0) {
      weights <- rep(1 / n, n)
    } else {
      weights <- local_variance / variance_sum
    }
    
    gse_values[j] <- sqrt(sum(weights * local_rmse^2, na.rm = TRUE))
  }
  
  gse_df <- data.frame(
    K_Value = k_values,
    GSE = gse_values
  ) %>%
    na.omit() %>%
    arrange(K_Value)
  
  loess_model <- loess(GSE ~ K_Value, data = gse_df, span = 0.75)
  
  y_fit <- predict(
    loess_model,
    newdata = data.frame(K_Value = gse_df$K_Value)
  )
  
  fit_df <- data.frame(
    K = gse_df$K_Value,
    y_fit = y_fit
  ) %>%
    na.omit()
  
  y_prev <- fit_df$y_fit[-nrow(fit_df)]
  y_next <- fit_df$y_fit[-1]
  k_next <- fit_df$K[-1]
  
  marginal_rate <- ((y_next - y_prev) / y_prev) * 100
  marginal_rate[!is.finite(marginal_rate)] <- NA
  
  rate_df <- data.frame(
    K_next = k_next,
    marginal_rate = marginal_rate
  ) %>%
    na.omit()
  
  threshold_value <- -1
  optimal_k <- NA_real_
  
  if (nrow(rate_df) >= 3) {
    
    mr <- rate_df$marginal_rate
    
    for (i in seq_len(length(mr) - 2)) {
      
      if (
        mr[i] > threshold_value &&
        mr[i + 1] > threshold_value &&
        mr[i + 2] > threshold_value
      ) {
        optimal_k <- rate_df$K_next[i]
        break
      }
    }
  }
  
  if (is.na(optimal_k)) {
    optimal_k <- fit_df$K[which.min(fit_df$y_fit)]
  }
  
  gse_opt <- gse_df$GSE[gse_df$K_Value == optimal_k][1]
  
  return(list(
    optimal_k = optimal_k,
    gse_opt = gse_opt,
    gse_table = gse_df
  ))
}

# ---- 8. Check outputs ----
rmse_value <- calculate_rmse(residual_data)
gse_fixed_value <- calculate_gse_fixed_k(residual_data, k = fixed_k)
gseopt_result <- calculate_gseopt(residual_data, k_values)

message("RMSE: ", round(rmse_value, 6))
message("Fixed-K GSE: ", round(gse_fixed_value, 6))
message("GSEopt K: ", gseopt_result$optimal_k)
message("GSEopt: ", round(gseopt_result$gse_opt, 6))

# ---- 9. Runtime comparison ----
time_compare <- microbenchmark(
  RMSE = calculate_rmse(residual_data),
  GSE_fixed_K = calculate_gse_fixed_k(residual_data, k = fixed_k),
  GSEopt = calculate_gseopt(residual_data, k_values),
  times = 10
)

print(time_compare)

# ---- 10. Summarize runtime ----
time_summary <- summary(time_compare)

runtime_summary <- data.frame(
  Metric = as.character(time_summary$expr),
  Median_seconds = time_summary$median / 1e9,
  Mean_seconds = time_summary$mean / 1e9
)

rmse_time <- runtime_summary$Median_seconds[
  runtime_summary$Metric == "RMSE"
]

gse_fixed_time <- runtime_summary$Median_seconds[
  runtime_summary$Metric == "GSE_fixed_K"
]

gseopt_time <- runtime_summary$Median_seconds[
  runtime_summary$Metric == "GSEopt"
]

runtime_ratios <- data.frame(
  Ratio = c(
    "GSE_fixed_K / RMSE",
    "GSEopt / RMSE",
    "GSEopt / GSE_fixed_K"
  ),
  Value = c(
    gse_fixed_time / rmse_time,
    gseopt_time / rmse_time,
    gseopt_time / gse_fixed_time
  )
)

runtime_ratios$Value <- round(runtime_ratios$Value, 2)

# ---- 11. Export results ----
write.csv(
  runtime_summary,
  file = file.path(output_dir, "runtime_comparison_results.csv"),
  row.names = FALSE
)

write.csv(
  runtime_ratios,
  file = file.path(output_dir, "runtime_comparison_ratios.csv"),
  row.names = FALSE
)

print(runtime_summary)
print(runtime_ratios)

message("Runtime comparison results saved to: ", output_dir)
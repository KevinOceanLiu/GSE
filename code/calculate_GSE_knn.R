# =========================================================
# 03_calculate_GSE_knn.R
# Calculate k-NN-based GSE for all models
# =========================================================

# ---- 1. Load packages ----
library(FNN)
library(sf)

# ---- 2. Set paths ----
data_dir <- "result"
output_dir <- "result"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ---- 3. User parameters ----
model_list <- c(
  "lm", "rf", "xgbTree", "gbm", "knn",
  "gam", "gwr", "svmRadial", "cubist"
)

k_values <- c(
  3, 5, 7, 9, 11, 13, 15, 17, 19, 21,
  23, 25, 27, 29, 31, 33, 35, 37, 39, 41,
  43, 45, 47, 49, 51, 53, 55, 57, 59, 61,
  63, 65, 67, 69, 71, 73, 75, 77, 79, 81,
  83, 85, 87, 89, 91, 93, 95, 97, 99, 101
)

# ---- 4. Project coordinates if needed ----
add_projected_coordinates <- function(df) {
  
  if (all(c("x_proj", "y_proj") %in% names(df))) {
    return(df)
  }
  
  required_coord_cols <- c("longitude", "latitude")
  missing_coord_cols <- setdiff(required_coord_cols, names(df))
  
  if (length(missing_coord_cols) > 0) {
    stop(paste(
      "Missing coordinate columns:",
      paste(missing_coord_cols, collapse = ", ")
    ))
  }
  
  df_sf <- st_as_sf(
    df,
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  )
  
  df_sf_3577 <- st_transform(df_sf, crs = 3577)
  xy <- st_coordinates(df_sf_3577)
  
  df$x_proj <- xy[, 1]
  df$y_proj <- xy[, 2]
  
  return(df)
}

# ---- 5. Standardize residual column name ----
standardize_residual_column <- function(df) {
  
  if ("residuals" %in% names(df)) {
    df$residual <- df$residuals
  }
  
  if (!"residual" %in% names(df)) {
    stop("Residual column not found. Expected 'residual' or 'residuals'.")
  }
  
  return(df)
}

# ---- 6. Calculate GSE for one K ----
calculate_gse_knn <- function(residual_data, k) {
  
  if (k < 2) {
    return(NA_real_)
  }
  
  if (nrow(residual_data) < k) {
    return(NA_real_)
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

# ---- 7. Calculate GSE for all models ----
gse_results_list <- list()

for (model_name in model_list) {
  
  message("Processing model: ", model_name)
  
  residual_file <- file.path(
    data_dir,
    paste0(model_name, "_test_residuals.csv")
  )
  
  if (!file.exists(residual_file)) {
    warning(paste("Residual file not found:", residual_file))
    next
  }
  
  residual_data <- read.csv(residual_file, fileEncoding = "UTF-8")
  residual_data <- standardize_residual_column(residual_data)
  residual_data <- add_projected_coordinates(residual_data)
  
  gse_values <- numeric(length(k_values))
  
  for (i in seq_along(k_values)) {
    
    k <- k_values[i]
    
    gse_values[i] <- calculate_gse_knn(
      residual_data = residual_data,
      k = k
    )
  }
  
  gse_results_list[[model_name]] <- gse_values
}

# ---- 8. Export results ----
if (length(gse_results_list) == 0) {
  stop("No GSE results were generated. Please check residual files.")
}

gse_results <- data.frame(
  K_Value = k_values,
  as.data.frame(gse_results_list)
)

output_file <- file.path(output_dir, "GSE_results_kNN2.csv")

write.csv(
  gse_results,
  file = output_file,
  row.names = FALSE
)

print(gse_results)

message("GSE results saved to: ", output_file)
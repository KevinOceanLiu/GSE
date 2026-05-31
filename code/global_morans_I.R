# =========================================================
# 06_global_morans_I.R
# Calculate Global Moran's I for test-set residuals
# =========================================================

# ---- 1. Load packages ----
library(spdep)
library(dplyr)

# ---- 2. Set paths ----
residual_dir <- "result"
output_dir <- "result"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ---- 3. User parameters ----
model_list <- c(
  "rf", "svmRadial", "gam", "lm", "xgbTree",
  "cubist", "gbm", "knn", "gwr"
)

# Model-specific optimal K values selected from GSE-K curves
optimal_k <- c(
  rf = 25,
  svmRadial = 23,
  gam = 21,
  lm = 19,
  xgbTree = 23,
  cubist = 25,
  gbm = 23,
  knn = 25,
  gwr = 27
)

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

# ---- 5. Calculate Global Moran's I ----
moran_results <- list()

for (model_name in model_list) {
  
  message("Calculating Moran's I for model: ", model_name)
  
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
  
  required_cols <- c("ID", "x_proj", "y_proj", "residual")
  missing_cols <- setdiff(required_cols, names(residual_data))
  
  if (length(missing_cols) > 0) {
    stop(paste(
      "Missing columns in", model_name, ":",
      paste(missing_cols, collapse = ", ")
    ))
  }
  
  residual_data <- residual_data[order(residual_data$ID), ]
  
  coords <- as.matrix(residual_data[, c("x_proj", "y_proj")])
  k_use <- optimal_k[[model_name]]
  
  if (is.null(k_use) || is.na(k_use)) {
    stop(paste("Optimal K is missing for model:", model_name))
  }
  
  # Build k-NN spatial weights using the same K as GSE
  knn_obj <- knearneigh(coords, k = k_use)
  nb_obj <- knn2nb(knn_obj)
  listw_obj <- nb2listw(nb_obj, style = "W", zero.policy = TRUE)
  
  # One-sided Moran's I test for positive spatial autocorrelation
  moran_test <- moran.test(
    residual_data$residual,
    listw_obj,
    zero.policy = TRUE,
    alternative = "greater"
  )
  
  moran_results[[model_name]] <- data.frame(
    Model = model_name,
    Optimal_K = k_use,
    Moran_I = as.numeric(moran_test$estimate["Moran I statistic"]),
    Expected_I = as.numeric(moran_test$estimate["Expectation"]),
    Variance = as.numeric(moran_test$estimate["Variance"]),
    p_value = moran_test$p.value
  )
}

# ---- 6. Export summary ----
if (length(moran_results) == 0) {
  stop("No Moran's I results were generated. Please check residual files.")
}

moran_summary <- bind_rows(moran_results) %>%
  mutate(
    Moran_I = round(Moran_I, 4),
    Expected_I = round(Expected_I, 4),
    Variance = round(Variance, 6),
    p_value = round(p_value, 4)
  )

output_file <- file.path(output_dir, "global_morans_I_results.csv")

write.csv(
  moran_summary,
  file = output_file,
  row.names = FALSE
)

print(moran_summary)

message("Global Moran's I results saved to: ", output_file)
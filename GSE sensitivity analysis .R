# =========================================================
# Figure 10 plotting data calculation
# Shared-version script
# Computational method unchanged
#
# Function:
# Refit the nine models using the same fixed train/test split,
# then recalculate RMSE and GSEopt for area-, richness-, and
# elevation-based subsets on the common test set.
# =========================================================

# -------------------------
# 0. Load packages
# -------------------------
library(caret)
library(dplyr)
library(sf)
library(FNN)

library(randomForest)
library(xgboost)
library(gbm)
library(mgcv)
library(pls)
library(e1071)
library(Cubist)
library(kknn)

# -------------------------
# 1. User settings
# -------------------------
data_path <- "data.csv"
output_dir <- "GSE_outputs"
output_path <- file.path(output_dir, "figure10_plot_data.csv")

MODEL_LIST <- c(
  "lm", "rf", "xgbTree", "gbm", "knn",
  "gam", "pls", "svmRadial", "cubist"
)

# K range for GSE calculation
K_VALUES <- c(
  3, 5, 7, 9, 11, 13, 15, 17, 19, 21,
  23, 25, 27, 29, 31, 33, 35, 37, 39, 41,
  43, 45, 47, 49, 51, 53, 55, 57, 59, 61,
  63, 65, 67, 69, 71, 73, 75, 77, 79, 81,
  83, 85, 87, 89, 91, 93, 95, 97, 99, 101
)

CRS_CODE <- 9473   # GDA2020 / Australian Albers

# Current GSEopt rule
# first consecutive 3 marginal rate > -1
threshold_value <- -1
n_consecutive <- 3

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -------------------------
# 2. Read data
# -------------------------
data <- read.csv(
  data_path,
  fileEncoding = "UTF-8",
  check.names = FALSE
)

# Clean column names
names(data) <- trimws(names(data))

# Compatibility for dot-style names
names(data)[names(data) == "Elevation.Division"] <- "Elevation Division"
names(data)[names(data) == "richness.Division"]  <- "richness Division"
names(data)[names(data) == "area.Division"]      <- "area Division"

cat("Current column names:\n")
print(names(data))

# -------------------------
# 3. Check required columns
# -------------------------
required_cols <- c(
  "ID", "richness", "longitude", "latitude",
  "Short.wave", "Wind.speed", "Soil.depth", "Soil.total.nitrogen",
  "Soil.pH", "Soil.CEC", "Elevation", "Slope",
  "sinAspect", "cosAspect",
  "Distance.to.artificial.land", "Distance.to.water",
  "Elevation Division", "richness Division", "area Division"
)

missing_cols <- setdiff(required_cols, names(data))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

# -------------------------
# 4. Remove missing values and sort by ID
# -------------------------
data <- na.omit(data)
data <- data[order(data$ID), ]

# -------------------------
# 5. Fixed train/test split based on ID
# Same logic as the main modeling script
# -------------------------
set.seed(123)

unique_ids <- sort(unique(data$ID))
train_id_index <- createDataPartition(unique_ids, p = 0.7, list = FALSE)
train_ids <- unique_ids[train_id_index]
test_ids  <- setdiff(unique_ids, train_ids)

trainData <- data[data$ID %in% train_ids, ]
testData  <- data[data$ID %in% test_ids, ]

trainData <- trainData[order(trainData$ID), ]
testData  <- testData[order(testData$ID), ]

cat("Train size:", nrow(trainData), "\n")
cat("Test size:", nrow(testData), "\n")

# -------------------------
# 6. Fixed 5-fold CV on training data
# -------------------------
set.seed(123)
folds <- createFolds(trainData$richness, k = 5, returnTrain = TRUE)

# -------------------------
# 7. Fix caret seeds
# -------------------------
tune_length <- 3

set.seed(123)
seeds <- vector(mode = "list", length = 5 + 1)
for (i in 1:5) {
  seeds[[i]] <- sample.int(100000, tune_length)
}
seeds[[6]] <- sample.int(100000, 1)

ctrl <- trainControl(
  method = "cv",
  number = 5,
  index = folds,
  allowParallel = FALSE,
  seeds = seeds
)

# -------------------------
# 8. Predictor formula
# IMPORTANT:
# Division columns are used only for stratification,
# not as predictors in model fitting
# -------------------------
predictor_formula <- richness ~ . - ID - longitude - latitude -
  `Elevation Division` - `richness Division` - `area Division`

# -------------------------
# 9. Helper functions
# -------------------------
prepare_projected_coords <- function(df, crs_code = 9473) {
  pts_sf <- st_as_sf(
    df,
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  )
  pts_proj <- st_transform(pts_sf, crs = crs_code)
  st_coordinates(pts_proj)
}

calculate_gse_knn <- function(residual_data, coords_matrix, k) {
  
  if (k < 2) return(NA_real_)
  if (nrow(residual_data) < k) return(NA_real_)
  
  nn_indices <- get.knnx(coords_matrix, coords_matrix, k = k)$nn.index
  
  num_points <- nrow(residual_data)
  local_rmse_k <- numeric(num_points)
  local_variance_k <- numeric(num_points)
  
  for (i in seq_len(num_points)) {
    current_indices <- nn_indices[i, ]
    local_errors <- residual_data$residuals[current_indices]
    
    local_rmse_k[i] <- sqrt(mean(local_errors^2))
    
    if (length(local_errors) > 1) {
      local_variance_k[i] <- var(local_errors) * ((k - 1) / k)
    } else {
      local_variance_k[i] <- 0
    }
  }
  
  total_variance_sum <- sum(local_variance_k)
  
  if (total_variance_sum == 0) {
    weights_k <- rep(1 / num_points, num_points)
  } else {
    weights_k <- local_variance_k / total_variance_sum
  }
  
  sqrt(sum(weights_k * local_rmse_k^2))
}

find_gseopt <- function(gse_df,
                        threshold_value = -1,
                        n_consecutive = 3) {
  
  gse_df <- gse_df %>%
    filter(!is.na(GSE)) %>%
    arrange(K_Value)
  
  if (nrow(gse_df) < (n_consecutive + 1)) {
    return(list(optimal_k = NA_real_, gse_opt = NA_real_))
  }
  
  loess_model <- loess(GSE ~ K_Value, data = gse_df, span = 0.75)
  fitted_vals <- predict(loess_model, newdata = data.frame(K_Value = gse_df$K_Value))
  
  fit_df <- data.frame(
    K = gse_df$K_Value,
    GSE_fit = fitted_vals
  ) %>%
    na.omit()
  
  if (nrow(fit_df) < (n_consecutive + 1)) {
    return(list(optimal_k = NA_real_, gse_opt = NA_real_))
  }
  
  y_prev <- fit_df$GSE_fit[-nrow(fit_df)]
  y_next <- fit_df$GSE_fit[-1]
  k_next <- fit_df$K[-1]
  
  marginal_rate <- ((y_next - y_prev) / y_prev) * 100
  marginal_rate[!is.finite(marginal_rate)] <- NA
  
  mr_df <- data.frame(
    K_next = k_next,
    marginal_rate = marginal_rate
  ) %>%
    na.omit()
  
  optimal_k <- NA_real_
  
  if (nrow(mr_df) >= n_consecutive) {
    mr <- mr_df$marginal_rate
    
    for (i in 1:(length(mr) - n_consecutive + 1)) {
      if (all(mr[i:(i + n_consecutive - 1)] > threshold_value)) {
        optimal_k <- mr_df$K_next[i]
        break
      }
    }
  }
  
  if (is.na(optimal_k)) {
    optimal_k <- fit_df$K[which.min(fit_df$GSE_fit)]
  }
  
  gse_opt <- gse_df$GSE[gse_df$K_Value == optimal_k][1]
  
  list(optimal_k = optimal_k, gse_opt = gse_opt)
}

calc_rmse <- function(obs, pred) {
  sqrt(mean((obs - pred)^2))
}

train_one_model <- function(model_name, trainData, ctrl, tune_length = 3) {
  
  if (model_name == "xgbTree") {
    tune_grid <- expand.grid(
      nrounds = 100,
      max_depth = 3,
      eta = 0.1,
      gamma = 0,
      colsample_bytree = 0.8,
      min_child_weight = 1,
      subsample = 0.8
    )
  } else {
    tune_grid <- NULL
  }
  
  set.seed(123)
  
  if (model_name == "xgbTree") {
    train(
      predictor_formula,
      data = trainData,
      method = model_name,
      trControl = ctrl,
      tuneGrid = tune_grid,
      verbose = FALSE
    )
  } else if (model_name == "rf") {
    train(
      predictor_formula,
      data = trainData,
      method = model_name,
      trControl = ctrl,
      tuneLength = tune_length,
      ntree = 500
    )
  } else if (model_name == "gbm") {
    train(
      predictor_formula,
      data = trainData,
      method = model_name,
      trControl = ctrl,
      tuneLength = tune_length,
      verbose = FALSE
    )
  } else if (model_name == "cubist") {
    train(
      predictor_formula,
      data = trainData,
      method = model_name,
      trControl = ctrl,
      tuneLength = tune_length
    )
  } else {
    train(
      predictor_formula,
      data = trainData,
      method = model_name,
      trControl = ctrl,
      tuneLength = tune_length
    )
  }
}

process_one_subset <- function(df_subset, model_name, strat_type, subset_name) {
  
  if (nrow(df_subset) < 10) {
    return(NULL)
  }
  
  rmse_value <- calc_rmse(df_subset$richness_true, df_subset$richness_pred)
  
  k_use <- K_VALUES[K_VALUES <= nrow(df_subset)]
  
  if (length(k_use) < 4) {
    return(data.frame(
      Model = model_name,
      Stratification_Type = strat_type,
      Subset = as.character(subset_name),
      n_points = nrow(df_subset),
      RMSE = rmse_value,
      GSEopt = NA_real_,
      Optimal_K = NA_real_
    ))
  }
  
  coords_matrix <- prepare_projected_coords(df_subset, crs_code = CRS_CODE)
  
  gse_vals <- sapply(k_use, function(k) {
    calculate_gse_knn(
      residual_data = df_subset,
      coords_matrix = coords_matrix,
      k = k
    )
  })
  
  gse_df <- data.frame(
    K_Value = k_use,
    GSE = gse_vals
  )
  
  opt_res <- find_gseopt(
    gse_df = gse_df,
    threshold_value = threshold_value,
    n_consecutive = n_consecutive
  )
  
  data.frame(
    Model = model_name,
    Stratification_Type = strat_type,
    Subset = as.character(subset_name),
    n_points = nrow(df_subset),
    RMSE = rmse_value,
    GSEopt = opt_res$gse_opt,
    Optimal_K = opt_res$optimal_k
  )
}

# -------------------------
# 10. Main loop
# -------------------------
all_results <- list()

for (model_name in MODEL_LIST) {
  
  cat("\n=============================\n")
  cat("Running model:", model_name, "\n")
  cat("=============================\n")
  
  model <- tryCatch({
    train_one_model(
      model_name = model_name,
      trainData = trainData,
      ctrl = ctrl,
      tune_length = tune_length
    )
  }, error = function(e) {
    cat("Model failed:", model_name, "\n")
    cat("Reason:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(model)) next
  
  test_pred <- predict(model, newdata = testData)
  
  test_result_df <- testData %>%
    mutate(
      richness_true = richness,
      richness_pred = test_pred,
      residuals = richness - test_pred
    )
  
  # area subsets
  area_levels <- sort(unique(test_result_df$`area Division`))
  for (lv in area_levels) {
    df_sub <- test_result_df %>% filter(`area Division` == lv)
    res <- process_one_subset(df_sub, model_name, "area", lv)
    if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
  }
  
  # richness subsets
  richness_levels <- sort(unique(test_result_df$`richness Division`))
  for (lv in richness_levels) {
    df_sub <- test_result_df %>% filter(`richness Division` == lv)
    res <- process_one_subset(df_sub, model_name, "richness", lv)
    if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
  }
  
  # elevation subsets
  elev_levels <- sort(unique(test_result_df$`Elevation Division`))
  for (lv in elev_levels) {
    df_sub <- test_result_df %>% filter(`Elevation Division` == lv)
    res <- process_one_subset(df_sub, model_name, "elevation", lv)
    if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
  }
}

# -------------------------
# 11. Save output
# -------------------------
if (length(all_results) == 0) {
  stop("No results were generated.")
}

figure10_plot_data <- bind_rows(all_results)

write.csv(figure10_plot_data, output_path, row.names = FALSE)

cat("\nDone.\n")
cat("Figure 10 plotting data saved to:\n", output_path, "\n")
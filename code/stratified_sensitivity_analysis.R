# =========================================================
# 07_stratified_sensitivity_analysis.R
# Calculate stratified RMSE and GSEopt for Figure 10
# =========================================================

# ---- 1. Load packages ----
library(caret)
library(dplyr)
library(sf)
library(FNN)

library(randomForest)
library(xgboost)
library(gbm)
library(mgcv)
library(e1071)
library(Cubist)
library(kknn)

library(sp)
library(GWmodel)

# ---- 2. Set paths ----
data_path <- "data/data3_partitions.csv"
output_dir <- "result"
output_file <- file.path(output_dir, "figure10_plot_data.csv")

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

# GDA2020 / Australian Albers
crs_code <- 9473

# Fixed adaptive GWR bandwidth selected from prior CV search
gwr_bw <- 227

# GSEopt rule: first K at which three consecutive marginal rates are greater than -1%
threshold_value <- -1
n_consecutive <- 3

set.seed(123)

# ---- 4. Read data ----
data <- read.csv(
  data_path,
  fileEncoding = "UTF-8",
  check.names = FALSE
)

names(data) <- trimws(names(data))

# Allow both original names and R-converted names
names(data)[names(data) == "Elevation.Division"] <- "Elevation Division"
names(data)[names(data) == "richness.Division"] <- "richness Division"
names(data)[names(data) == "area.Division"] <- "area Division"

# ---- 5. Check required columns ----
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
  stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
}

# ---- 6. Prepare data ----
data <- na.omit(data)
data <- data[order(data$ID), ]

# ---- 7. Fixed train/test split by ID ----
set.seed(123)

unique_ids <- sort(unique(data$ID))
train_id_index <- createDataPartition(unique_ids, p = 0.7, list = FALSE)

train_ids <- unique_ids[train_id_index]
test_ids <- setdiff(unique_ids, train_ids)

trainData <- data[data$ID %in% train_ids, ]
testData <- data[data$ID %in% test_ids, ]

trainData <- trainData[order(trainData$ID), ]
testData <- testData[order(testData$ID), ]

message("Train size: ", nrow(trainData))
message("Test size: ", nrow(testData))

# ---- 8. Fixed 5-fold cross-validation ----
set.seed(123)

folds <- createFolds(trainData$richness, k = 5, returnTrain = TRUE)

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

# ---- 9. Model formula ----
# Division columns are excluded from model predictors
predictor_formula <- richness ~ . - ID - longitude - latitude -
  `Elevation Division` - `richness Division` - `area Division`

predictors <- c(
  "Short.wave", "Wind.speed", "Soil.depth", "Soil.total.nitrogen",
  "Soil.pH", "Soil.CEC", "Elevation", "Slope",
  "sinAspect", "cosAspect",
  "Distance.to.artificial.land", "Distance.to.water"
)

# ---- 10. Helper functions ----
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

add_projected_coords <- function(df, crs_code = 9473) {
  
  xy <- prepare_projected_coords(df, crs_code)
  
  df$x_proj <- xy[, 1]
  df$y_proj <- xy[, 2]
  
  df
}

calculate_gse_knn <- function(residual_data, coords_matrix, k) {
  
  if (k < 2) {
    return(NA_real_)
  }
  
  if (nrow(residual_data) < k) {
    return(NA_real_)
  }
  
  nn_indices <- get.knnx(coords_matrix, coords_matrix, k = k)$nn.index
  
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

find_gseopt <- function(
    gse_df,
    threshold_value = -1,
    n_consecutive = 3
) {
  
  gse_df <- gse_df %>%
    filter(!is.na(GSE)) %>%
    arrange(K_Value)
  
  if (nrow(gse_df) < (n_consecutive + 1)) {
    return(list(optimal_k = NA_real_, gse_opt = NA_real_))
  }
  
  loess_model <- loess(GSE ~ K_Value, data = gse_df, span = 0.75)
  
  fitted_vals <- predict(
    loess_model,
    newdata = data.frame(K_Value = gse_df$K_Value)
  )
  
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
    
    for (i in seq_len(length(mr) - n_consecutive + 1)) {
      
      current_window <- mr[i:(i + n_consecutive - 1)]
      
      if (all(current_window > threshold_value)) {
        optimal_k <- mr_df$K_next[i]
        break
      }
    }
  }
  
  if (is.na(optimal_k)) {
    optimal_k <- fit_df$K[which.min(fit_df$GSE_fit)]
  }
  
  gse_opt <- gse_df$GSE[gse_df$K_Value == optimal_k][1]
  
  return(list(optimal_k = optimal_k, gse_opt = gse_opt))
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

predict_gwr <- function(
    trainData,
    testData,
    predictors,
    gwr_bw = 227,
    crs_code = 9473
) {
  
  trainData_gwr <- add_projected_coords(trainData, crs_code)
  testData_gwr <- add_projected_coords(testData, crs_code)
  
  gwr_formula <- as.formula(
    paste("richness ~", paste(predictors, collapse = " + "))
  )
  
  train_sp <- trainData_gwr
  coordinates(train_sp) <- ~ x_proj + y_proj
  proj4string(train_sp) <- CRS(SRS_string = paste0("EPSG:", crs_code))
  
  test_sp <- testData_gwr
  coordinates(test_sp) <- ~ x_proj + y_proj
  proj4string(test_sp) <- CRS(SRS_string = paste0("EPSG:", crs_code))
  
  gwr_test <- gwr.basic(
    formula = gwr_formula,
    data = train_sp,
    regression.points = test_sp,
    bw = gwr_bw,
    kernel = "bisquare",
    adaptive = TRUE,
    longlat = FALSE
  )
  
  required_coef_cols <- c("Intercept", predictors)
  missing_coef <- setdiff(required_coef_cols, names(gwr_test$SDF))
  
  if (length(missing_coef) > 0) {
    stop(paste(
      "Missing GWR coefficient columns:",
      paste(missing_coef, collapse = ", ")
    ))
  }
  
  test_coef <- as.data.frame(gwr_test$SDF)[, required_coef_cols]
  test_x <- testData_gwr[, predictors]
  
  test_pred <- test_coef$Intercept +
    rowSums(as.matrix(test_coef[, predictors]) * as.matrix(test_x))
  
  return(test_pred)
}

process_one_subset <- function(
    df_subset,
    model_name,
    strat_type,
    subset_name
) {
  
  if (nrow(df_subset) < 10) {
    return(NULL)
  }
  
  rmse_value <- calc_rmse(
    obs = df_subset$richness_true,
    pred = df_subset$richness_pred
  )
  
  k_use <- k_values[k_values <= nrow(df_subset)]
  
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
  
  coords_matrix <- prepare_projected_coords(
    df_subset,
    crs_code = crs_code
  )
  
  gse_values <- sapply(k_use, function(k) {
    
    calculate_gse_knn(
      residual_data = df_subset,
      coords_matrix = coords_matrix,
      k = k
    )
  })
  
  gse_df <- data.frame(
    K_Value = k_use,
    GSE = gse_values
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

# ---- 11. Main analysis ----
all_results <- list()

for (model_name in model_list) {
  
  message("Running model: ", model_name)
  
  if (model_name == "gwr") {
    
    test_pred <- tryCatch(
      predict_gwr(
        trainData = trainData,
        testData = testData,
        predictors = predictors,
        gwr_bw = gwr_bw,
        crs_code = crs_code
      ),
      error = function(e) {
        warning(paste("GWR failed:", e$message))
        return(NULL)
      }
    )
    
    if (is.null(test_pred)) {
      next
    }
    
  } else {
    
    model <- tryCatch(
      train_one_model(
        model_name = model_name,
        trainData = trainData,
        ctrl = ctrl,
        tune_length = tune_length
      ),
      error = function(e) {
        warning(paste("Model failed:", model_name, e$message))
        return(NULL)
      }
    )
    
    if (is.null(model)) {
      next
    }
    
    test_pred <- predict(model, newdata = testData)
  }
  
  test_result_df <- testData %>%
    mutate(
      richness_true = richness,
      richness_pred = as.numeric(test_pred),
      residuals = richness - as.numeric(test_pred)
    )
  
  # Area-based stratification
  area_levels <- sort(unique(test_result_df$`area Division`))
  
  for (lv in area_levels) {
    
    df_sub <- test_result_df %>%
      filter(`area Division` == lv)
    
    res <- process_one_subset(
      df_subset = df_sub,
      model_name = model_name,
      strat_type = "area",
      subset_name = lv
    )
    
    if (!is.null(res)) {
      all_results[[length(all_results) + 1]] <- res
    }
  }
  
  # Richness-based stratification
  richness_levels <- sort(unique(test_result_df$`richness Division`))
  
  for (lv in richness_levels) {
    
    df_sub <- test_result_df %>%
      filter(`richness Division` == lv)
    
    res <- process_one_subset(
      df_subset = df_sub,
      model_name = model_name,
      strat_type = "richness",
      subset_name = lv
    )
    
    if (!is.null(res)) {
      all_results[[length(all_results) + 1]] <- res
    }
  }
  
  # Elevation-based stratification
  elevation_levels <- sort(unique(test_result_df$`Elevation Division`))
  
  for (lv in elevation_levels) {
    
    df_sub <- test_result_df %>%
      filter(`Elevation Division` == lv)
    
    res <- process_one_subset(
      df_subset = df_sub,
      model_name = model_name,
      strat_type = "elevation",
      subset_name = lv
    )
    
    if (!is.null(res)) {
      all_results[[length(all_results) + 1]] <- res
    }
  }
}

# ---- 12. Export results ----
if (length(all_results) == 0) {
  stop("No stratified sensitivity results were generated.")
}

figure10_plot_data <- bind_rows(all_results)

write.csv(
  figure10_plot_data,
  file = output_file,
  row.names = FALSE
)

print(figure10_plot_data)

message("Figure 10 plotting data saved to: ", output_file)
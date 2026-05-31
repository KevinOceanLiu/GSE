# =========================================================
# 01_model_training_and_residuals.R
# Train non-spatial models and export test-set residuals
# =========================================================

# ---- 1. Load packages ----
library(caret)
library(readr)
library(dplyr)
library(ggplot2)
library(grid)
library(sf)

library(randomForest)
library(xgboost)
library(gbm)
library(Cubist)
library(kknn)
library(mgcv)
library(kernlab)

# ---- 2. Set paths ----
data_path <- "data/data.csv"
output_dir <- "result"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ---- 3. Read data ----
data <- read.csv(data_path, fileEncoding = "UTF-8")

# ---- 4. Check required columns ----
required_cols <- c(
  "ID", "richness", "longitude", "latitude",
  "Short.wave", "Wind.speed", "Soil.depth", "Soil.total.nitrogen",
  "Soil.pH", "Soil.CEC", "Elevation", "Slope",
  "sinAspect", "cosAspect",
  "Distance.to.artificial.land", "Distance.to.water"
)

missing_cols <- setdiff(required_cols, names(data))

if (length(missing_cols) > 0) {
  stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
}

# ---- 5. Prepare data ----
data <- na.omit(data)
data <- data[order(data$ID), ]

# ---- 6. Add projected coordinates ----
# EPSG:4326 = longitude/latitude
# EPSG:3577 = GDA94 / Australian Albers

data_sf <- st_as_sf(
  data,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

data_sf_3577 <- st_transform(data_sf, crs = 3577)
data_xy <- st_coordinates(data_sf_3577)

data$x_proj <- data_xy[, 1]
data$y_proj <- data_xy[, 2]

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

# ---- 8. Fixed 5-fold cross-validation ----
set.seed(123)

folds <- createFolds(trainData$richness, k = 5, returnTrain = TRUE)

ctrl <- trainControl(
  method = "cv",
  number = 5,
  index = folds,
  allowParallel = FALSE
)

# ---- 9. Model list ----
model_names <- c(
  "lm",
  "rf",
  "xgbTree",
  "gbm",
  "knn",
  "gam",
  "svmRadial",
  "cubist"
)

tune_length <- 3

# ---- 10. Train models and export results ----
performance_list <- list()

for (model_name in model_names) {
  
  message("Running model: ", model_name)
  
  set.seed(123)
  
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
    
    model <- train(
      richness ~ . - ID - longitude - latitude - x_proj - y_proj,
      data = trainData,
      method = model_name,
      trControl = ctrl,
      tuneGrid = tune_grid,
      verbose = FALSE
    )
    
  } else if (model_name == "rf") {
    
    model <- train(
      richness ~ . - ID - longitude - latitude - x_proj - y_proj,
      data = trainData,
      method = model_name,
      trControl = ctrl,
      tuneLength = tune_length,
      ntree = 500
    )
    
  } else if (model_name == "gbm") {
    
    model <- train(
      richness ~ . - ID - longitude - latitude - x_proj - y_proj,
      data = trainData,
      method = model_name,
      trControl = ctrl,
      tuneLength = tune_length,
      verbose = FALSE
    )
    
  } else {
    
    model <- train(
      richness ~ . - ID - longitude - latitude - x_proj - y_proj,
      data = trainData,
      method = model_name,
      trControl = ctrl,
      tuneLength = tune_length
    )
  }
  
  # ---- 11. Predict ----
  train_pred <- predict(model, newdata = trainData)
  test_pred <- predict(model, newdata = testData)
  
  train_residual <- trainData$richness - as.numeric(train_pred)
  test_residual <- testData$richness - as.numeric(test_pred)
  
  train_results <- postResample(pred = train_pred, obs = trainData$richness)
  test_results <- postResample(pred = test_pred, obs = testData$richness)
  
  train_mae <- mean(abs(train_residual))
  test_mae <- mean(abs(test_residual))
  
  performance_list[[model_name]] <- data.frame(
    Model = model_name,
    Train_R2 = as.numeric(train_results["Rsquared"]),
    Train_RMSE = as.numeric(train_results["RMSE"]),
    Train_MAE = train_mae,
    Test_R2 = as.numeric(test_results["Rsquared"]),
    Test_RMSE = as.numeric(test_results["RMSE"]),
    Test_MAE = test_mae
  )
  
  # ---- 12. Export test-set residuals ----
  test_residuals <- data.frame(
    ID = testData$ID,
    longitude = testData$longitude,
    latitude = testData$latitude,
    x_proj = testData$x_proj,
    y_proj = testData$y_proj,
    richness_true = testData$richness,
    richness_pred = as.numeric(test_pred),
    observed = testData$richness,
    predicted = as.numeric(test_pred),
    residual = test_residual,
    residuals = test_residual
  )
  
  write.csv(
    test_residuals,
    file = file.path(output_dir, paste0(model_name, "_test_residuals.csv")),
    row.names = FALSE
  )
  
  # ---- 13. Export scatterplot ----
  test_r2 <- round(as.numeric(test_results["Rsquared"]), 3)
  test_rmse <- round(as.numeric(test_results["RMSE"]), 3)
  test_mae_round <- round(test_mae, 3)
  
  test_plot <- ggplot(
    data = data.frame(
      Actual = testData$richness,
      Predicted = as.numeric(test_pred)
    ),
    aes(x = Actual, y = Predicted)
  ) +
    geom_point(color = "#eec9b9", alpha = 0.5) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "black",
      linewidth = 1
    ) +
    labs(
      title = model_name,
      x = "Actual richness",
      y = "Predicted richness"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 20),
      panel.background = element_rect(fill = "white", color = NA),
      axis.text = element_text(size = 20, color = "black"),
      axis.title = element_text(size = 20, color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      axis.ticks = element_line(color = "black", linewidth = 1),
      axis.line.x.bottom = element_line(color = "black", linewidth = 1),
      axis.line.y.left = element_line(color = "black", linewidth = 1)
    ) +
    scale_x_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, by = 20)
    ) +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, by = 20)
    ) +
    annotate(
      "text",
      x = 0,
      y = 100,
      label = paste("R² =", test_r2),
      size = 7,
      hjust = 0,
      color = "black"
    ) +
    annotate(
      "text",
      x = 0,
      y = 90,
      label = paste("RMSE =", test_rmse),
      size = 7,
      hjust = 0,
      color = "black"
    ) +
    annotate(
      "text",
      x = 0,
      y = 80,
      label = paste("MAE =", test_mae_round),
      size = 7,
      hjust = 0,
      color = "black"
    )
  
  ggsave(
    filename = file.path(output_dir, paste0(model_name, "_scatterplot.png")),
    plot = test_plot,
    width = 5,
    height = 5,
    dpi = 300,
    bg = "white"
  )
}

# ---- 14. Export model performance ----
performance_table <- bind_rows(performance_list)

write.csv(
  performance_table,
  file = file.path(output_dir, "model_performance.csv"),
  row.names = FALSE
)

print(performance_table)
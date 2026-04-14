# =========================================================
# Stable modeling code for GSE study
# Shared-version script
# Computational method unchanged
# =========================================================

# -------------------------
# 0. Load packages
# -------------------------
library(caret)
library(readr)
library(Cubist)
library(dplyr)
library(rpart)
library(GD)
library(randomForest)
library(e1071)
library(xgboost)
library(class)
library(mgcv)
library(pls)
library(gbm)
library(ggplot2)
library(grid)

# -------------------------
# 1. User settings
# -------------------------
data_path <- "data.csv"
output_dir <- "GSE_outputs"
residual_dir <- file.path(output_dir, "residuals")

# Models used in the paper
MODEL_LIST <- c(
  "lm", "rf", "xgbTree", "gbm", "knn",
  "gam", "pls", "svmRadial", "cubist"
)

# -------------------------
# 2. Create output folders
# -------------------------
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(residual_dir)) dir.create(residual_dir, recursive = TRUE)

# -------------------------
# 3. Read data
# -------------------------
data <- read.csv(
  data_path,
  fileEncoding = "UTF-8",
  check.names = FALSE
)

# Clean column names
names(data) <- trimws(names(data))

# Compatibility for possible dot-version names
names(data)[names(data) == "Elevation.Division"] <- "Elevation Division"
names(data)[names(data) == "richness.Division"]  <- "richness Division"
names(data)[names(data) == "area.Division"]      <- "area Division"

# -------------------------
# 4. Check required columns
# -------------------------
required_cols <- c(
  "ID", "richness", "longitude", "latitude",
  "Short.wave", "Wind.speed", "Soil.depth", "Soil.total.nitrogen",
  "Soil.pH", "Soil.CEC", "Elevation", "Slope",
  "sinAspect", "cosAspect",
  "Distance.to.artificial.land", "Distance.to.water"
)

missing_cols <- setdiff(required_cols, names(data))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Optional stratification columns used only for sensitivity analysis
division_cols <- c("Elevation Division", "richness Division", "area Division")
existing_division_cols <- intersect(division_cols, names(data))

# -------------------------
# 5. Remove missing values and sort by ID
# -------------------------
data <- na.omit(data)
data <- data[order(data$ID), ]

# -------------------------
# 6. Fixed train/test split based on ID
# This keeps results stable even if row order changes
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

# Save split IDs for reproducibility
write.csv(
  data.frame(ID = train_ids),
  file.path(output_dir, "train_ids.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(ID = test_ids),
  file.path(output_dir, "test_ids.csv"),
  row.names = FALSE
)

# -------------------------
# 7. Fixed 5-fold CV on training data
# -------------------------
set.seed(123)
folds <- createFolds(trainData$richness, k = 5, returnTrain = TRUE)

# -------------------------
# 8. Fix caret seeds for reproducibility
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
# 9. Predictor formula
# IMPORTANT:
# Division columns are excluded from model fitting
# so they are not used as predictors
# -------------------------
if (length(existing_division_cols) > 0) {
  predictor_formula <- as.formula(
    paste(
      "richness ~ . - ID - longitude - latitude -",
      paste(sprintf("`%s`", existing_division_cols), collapse = " - ")
    )
  )
} else {
  predictor_formula <- richness ~ . - ID - longitude - latitude
}

cat("Model formula:\n")
print(predictor_formula)

# -------------------------
# 10. Helper function: train one model
# -------------------------
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
    model <- train(
      predictor_formula,
      data = trainData,
      method = model_name,
      trControl = ctrl,
      tuneGrid = tune_grid,
      verbose = FALSE
    )
  } else if (model_name == "rf") {
    model <- train(
      predictor_formula,
      data = trainData,
      method = model_name,
      trControl = ctrl,
      tuneLength = tune_length,
      ntree = 500
    )
  } else if (model_name == "gbm") {
    model <- train(
      predictor_formula,
      data = trainData,
      method = model_name,
      trControl = ctrl,
      tuneLength = tune_length,
      verbose = FALSE
    )
  } else if (model_name == "cubist") {
    model <- train(
      predictor_formula,
      data = trainData,
      method = model_name,
      trControl = ctrl,
      tuneLength = tune_length
    )
  } else {
    model <- train(
      predictor_formula,
      data = trainData,
      method = model_name,
      trControl = ctrl,
      tuneLength = tune_length
    )
  }

  return(model)
}

# -------------------------
# 11. Main loop: run all models
# -------------------------
all_metrics <- list()

for (model_name in MODEL_LIST) {

  cat("\n=============================\n")
  cat("Running model:", model_name, "\n")
  cat("=============================\n")

  model <- tryCatch(
    train_one_model(
      model_name = model_name,
      trainData = trainData,
      ctrl = ctrl,
      tune_length = tune_length
    ),
    error = function(e) {
      cat("Model failed:", model_name, "\n")
      cat("Reason:", e$message, "\n")
      return(NULL)
    }
  )

  if (is.null(model)) next

  # Predict
  train_pred <- predict(model, newdata = trainData)
  test_pred  <- predict(model, newdata = testData)

  # Residuals
  train_residuals <- trainData$richness - train_pred
  test_residuals  <- testData$richness - test_pred

  # Metrics
  train_results <- postResample(pred = train_pred, obs = trainData$richness)
  test_results  <- postResample(pred = test_pred, obs = testData$richness)

  # Save residual file for test set
  residual_df <- data.frame(
    ID = testData$ID,
    longitude = testData$longitude,
    latitude = testData$latitude,
    richness_true = testData$richness,
    richness_pred = test_pred,
    residuals = test_residuals
  )

  write.csv(
    residual_df,
    file.path(residual_dir, paste0(model_name, "_test_residuals.csv")),
    row.names = FALSE
  )

  # Collect metrics
  all_metrics[[length(all_metrics) + 1]] <- data.frame(
    Model = model_name,
    Train_R2 = round(train_results["Rsquared"], 3),
    Train_RMSE = round(train_results["RMSE"], 3),
    Train_MAE = round(train_results["MAE"], 3),
    Test_R2 = round(test_results["Rsquared"], 3),
    Test_RMSE = round(test_results["RMSE"], 3),
    Test_MAE = round(test_results["MAE"], 3)
  )

  # Optional: save test scatter plot
  test_plot <- ggplot(
    data = data.frame(
      Actual = testData$richness,
      Predicted = test_pred
    ),
    aes(x = Actual, y = Predicted)
  ) +
    geom_point(color = "#eec9b9", alpha = 0.5) +
    geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed",
      color = "black",
      linewidth = 0.8
    ) +
    labs(
      title = model_name,
      x = "Actual richness",
      y = "Predicted richness"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 16),
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 12, color = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      axis.line = element_line(color = "black", linewidth = 0.6)
    ) +
    annotate(
      "text", x = min(testData$richness, na.rm = TRUE), 
      y = max(test_pred, na.rm = TRUE),
      label = paste0("R² = ", round(test_results["Rsquared"], 3)),
      hjust = 0, size = 4.5
    ) +
    annotate(
      "text", x = min(testData$richness, na.rm = TRUE), 
      y = max(test_pred, na.rm = TRUE) * 0.92,
      label = paste0("RMSE = ", round(test_results["RMSE"], 3)),
      hjust = 0, size = 4.5
    ) +
    annotate(
      "text", x = min(testData$richness, na.rm = TRUE), 
      y = max(test_pred, na.rm = TRUE) * 0.84,
      label = paste0("MAE = ", round(test_results["MAE"], 3)),
      hjust = 0, size = 4.5
    )

  ggsave(
    filename = file.path(output_dir, paste0(model_name, "_test_scatter.png")),
    plot = test_plot,
    width = 4.7,
    height = 4.7,
    dpi = 300
  )
}

# -------------------------
# 12. Save all model metrics
# -------------------------
if (length(all_metrics) == 0) {
  stop("No model results were generated.")
}

model_metrics_df <- bind_rows(all_metrics)

write.csv(
  model_metrics_df,
  file.path(output_dir, "model_test_metrics.csv"),
  row.names = FALSE
)

cat("\nDone.\n")
cat("Outputs saved in:\n", output_dir, "\n")
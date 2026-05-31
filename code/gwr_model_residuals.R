# =========================================================
# 02_gwr_model_residuals.R
# Fit GWR and export test-set residuals
# =========================================================

# ---- 1. Load packages ----
library(caret)
library(dplyr)
library(ggplot2)
library(grid)
library(sf)
library(sp)
library(GWmodel)

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

# ---- 6. Fixed train/test split by ID ----
set.seed(123)

unique_ids <- sort(unique(data$ID))
train_id_index <- createDataPartition(unique_ids, p = 0.7, list = FALSE)

train_ids <- unique_ids[train_id_index]
test_ids <- setdiff(unique_ids, train_ids)

# ---- 7. Project coordinates for GWR ----
# EPSG:4326 = longitude/latitude
# EPSG:3577 = GDA94 / Australian Albers

all_sf <- st_as_sf(
  data,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

all_sf_3577 <- st_transform(all_sf, crs = 3577)
all_xy <- st_coordinates(all_sf_3577)

data$x_proj <- all_xy[, 1]
data$y_proj <- all_xy[, 2]

# ---- 8. Recreate train/test data after projection ----
trainData <- data[data$ID %in% train_ids, ]
testData <- data[data$ID %in% test_ids, ]

trainData <- trainData[order(trainData$ID), ]
testData <- testData[order(testData$ID), ]

# ---- 9. Define GWR formula ----
predictors <- c(
  "Short.wave", "Wind.speed", "Soil.depth", "Soil.total.nitrogen",
  "Soil.pH", "Soil.CEC", "Elevation", "Slope",
  "sinAspect", "cosAspect",
  "Distance.to.artificial.land", "Distance.to.water"
)

gwr_formula <- as.formula(
  paste("richness ~", paste(predictors, collapse = " + "))
)

# ---- 10. Convert data to SpatialPointsDataFrame ----
train_sp <- trainData
coordinates(train_sp) <- ~ x_proj + y_proj
proj4string(train_sp) <- CRS(SRS_string = "EPSG:3577")

test_sp <- testData
coordinates(test_sp) <- ~ x_proj + y_proj
proj4string(test_sp) <- CRS(SRS_string = "EPSG:3577")

# ---- 11. Fit GWR ----
# Fixed adaptive bandwidth selected by prior CV search
gwr_bw <- 227

cat("Using fixed adaptive GWR bandwidth:", gwr_bw, "\n")

gwr_train <- gwr.basic(
  formula = gwr_formula,
  data = train_sp,
  regression.points = train_sp,
  bw = gwr_bw,
  kernel = "bisquare",
  adaptive = TRUE,
  longlat = FALSE
)

gwr_test <- gwr.basic(
  formula = gwr_formula,
  data = train_sp,
  regression.points = test_sp,
  bw = gwr_bw,
  kernel = "bisquare",
  adaptive = TRUE,
  longlat = FALSE
)

# ---- 12. Calculate predictions from local coefficients ----
# yhat = Intercept + beta1*x1 + beta2*x2 + ...

required_coef_cols <- c("Intercept", predictors)

missing_train_coef <- setdiff(required_coef_cols, names(gwr_train$SDF))
missing_test_coef <- setdiff(required_coef_cols, names(gwr_test$SDF))

if (length(missing_train_coef) > 0) {
  stop(paste(
    "Missing GWR training coefficient columns:",
    paste(missing_train_coef, collapse = ", ")
  ))
}

if (length(missing_test_coef) > 0) {
  stop(paste(
    "Missing GWR test coefficient columns:",
    paste(missing_test_coef, collapse = ", ")
  ))
}

train_coef <- as.data.frame(gwr_train$SDF)[, required_coef_cols]
test_coef <- as.data.frame(gwr_test$SDF)[, required_coef_cols]

train_x <- trainData[, predictors]
test_x <- testData[, predictors]

train_pred <- train_coef$Intercept +
  rowSums(as.matrix(train_coef[, predictors]) * as.matrix(train_x))

test_pred <- test_coef$Intercept +
  rowSums(as.matrix(test_coef[, predictors]) * as.matrix(test_x))

# ---- 13. Calculate residuals and metrics ----
train_residuals <- trainData$richness - train_pred
test_residuals <- testData$richness - test_pred

train_results <- postResample(pred = train_pred, obs = trainData$richness)
test_results <- postResample(pred = test_pred, obs = testData$richness)

train_mae <- mean(abs(train_residuals))
test_mae <- mean(abs(test_residuals))

performance <- data.frame(
  Model = "gwr",
  Bandwidth = gwr_bw,
  Train_R2 = as.numeric(train_results["Rsquared"]),
  Train_RMSE = as.numeric(train_results["RMSE"]),
  Train_MAE = train_mae,
  Test_R2 = as.numeric(test_results["Rsquared"]),
  Test_RMSE = as.numeric(test_results["RMSE"]),
  Test_MAE = test_mae
)

print(performance)

# ---- 14. Export test-set residuals ----
test_residuals_output <- data.frame(
  ID = testData$ID,
  longitude = testData$longitude,
  latitude = testData$latitude,
  observed = testData$richness,
  predicted = as.numeric(test_pred),
  residual = as.numeric(test_residuals)
)

write.csv(
  test_residuals_output,
  file = file.path(output_dir, "gwr_test_residuals.csv"),
  row.names = FALSE
)

write.csv(
  performance,
  file = file.path(output_dir, "gwr_model_performance.csv"),
  row.names = FALSE
)

# ---- 15. Export scatterplot ----
test_r2 <- round(performance$Test_R2, 3)
test_rmse <- round(performance$Test_RMSE, 3)
test_mae_round <- round(performance$Test_MAE, 3)

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
    title = "gwr",
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
  filename = file.path(output_dir, "gwr_scatterplot.png"),
  plot = test_plot,
  width = 1410 / 300,
  height = 1410 / 300,
  dpi = 300,
  units = "in",
  bg = "white"
)

print(test_plot)
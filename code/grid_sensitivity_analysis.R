# =========================================================
# 08_grid_GSE_analysis.R
# Calculate square-grid and hexagonal-grid GSE for all models
# =========================================================

# ---- 1. Load packages ----
library(sf)
library(dplyr)

# ---- 2. Set paths ----
boundary_path <- "data/boundary/boundary.shp"
residual_dir <- "result"
output_dir <- "result/grid_GSE_outputs"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (!file.exists(boundary_path)) {
  stop(paste("Boundary shapefile not found:", boundary_path))
}

if (!dir.exists(residual_dir)) {
  stop(paste("Residual directory not found:", residual_dir))
}

# ---- 3. User parameters ----
crs_code <- 9473  # GDA2020 / Australian Albers

grid_sizes <- round(
  seq(50000, 2000000, length.out = 100),
  0
)

target_models <- c(
  "lm", "rf", "xgbTree", "gbm", "knn",
  "gam", "gwr", "svmRadial", "cubist"
)

# ---- 4. Read study-area boundary ----
study_area <- st_read(
  boundary_path,
  options = "ENCODING=UTF-8",
  quiet = TRUE
)

study_area <- st_transform(study_area, crs = crs_code)
study_area <- st_make_valid(study_area)

boundary_union <- st_union(study_area)

boundary_union <- st_sf(
  id = 1,
  geometry = st_sfc(boundary_union, crs = crs_code)
)

boundary_union <- st_make_valid(boundary_union)

# ---- 5. Read residual file list ----
residual_files <- list.files(
  residual_dir,
  pattern = "_test_residuals\\.csv$",
  full.names = TRUE
)

if (length(residual_files) == 0) {
  stop(paste("No residual files found in:", residual_dir))
}

model_names <- gsub(
  "_test_residuals\\.csv$",
  "",
  basename(residual_files)
)

keep_idx <- model_names %in% target_models

residual_files <- residual_files[keep_idx]
model_names <- model_names[keep_idx]

order_idx <- match(target_models, model_names)
valid_order_idx <- order_idx[!is.na(order_idx)]

residual_files <- residual_files[valid_order_idx]
model_names <- model_names[valid_order_idx]

if (length(residual_files) == 0) {
  stop("No residual files matched target model names.")
}

missing_models <- setdiff(target_models, model_names)

if (length(missing_models) > 0) {
  warning(paste(
    "Missing residual files for:",
    paste(missing_models, collapse = ", ")
  ))
}

message("Models to process:")
print(model_names)

# ---- 6. Read residual points ----
read_residual_points <- function(file_path, crs_code = 9473) {
  
  df <- read.csv(
    file_path,
    fileEncoding = "UTF-8",
    check.names = FALSE
  )
  
  # Support both old and new residual column names
  if ("residual" %in% names(df) && !"residuals" %in% names(df)) {
    df$residuals <- df$residual
  }
  
  required_cols <- c("longitude", "latitude", "residuals")
  missing_cols <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0) {
    stop(paste(
      "Residual file missing columns:",
      paste(missing_cols, collapse = ", ")
    ))
  }
  
  if (!"ID" %in% names(df)) {
    df$ID <- seq_len(nrow(df))
  }
  
  pts_sf <- st_as_sf(
    df,
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  )
  
  pts_proj <- st_transform(pts_sf, crs = crs_code)
  pts_proj <- st_make_valid(pts_proj)
  
  return(pts_proj)
}

# ---- 7. Build populated grid ----
build_populated_grid <- function(base_map, points, cellsize, square = TRUE) {
  
  grid_geom <- st_make_grid(
    base_map,
    cellsize = cellsize,
    square = square
  )
  
  grid_sf <- st_sf(
    grid_id = seq_along(grid_geom),
    geometry = grid_geom
  )
  
  grid_sf <- st_make_valid(grid_sf)
  
  intersects_list <- st_intersects(
    grid_sf,
    points,
    sparse = TRUE
  )
  
  populated_idx <- which(lengths(intersects_list) > 0)
  
  if (length(populated_idx) == 0) {
    return(NULL)
  }
  
  populated_grid <- grid_sf[populated_idx, ]
  
  populated_grid <- st_intersection(
    populated_grid,
    boundary_union
  )
  
  populated_grid <- st_make_valid(populated_grid)
  populated_grid <- st_collection_extract(
    populated_grid,
    "POLYGON",
    warn = FALSE
  )
  
  populated_grid <- populated_grid[!st_is_empty(populated_grid), ]
  
  if (nrow(populated_grid) == 0) {
    return(NULL)
  }
  
  return(populated_grid)
}

# ---- 8. Calculate GSE for one grid system ----
calculate_gse_grid <- function(populated_grid, points_proj) {
  
  empty_result <- list(
    GSE = NA_real_,
    n_populated_cells = 0,
    n_points_used = 0,
    mean_points_per_cell = NA_real_
  )
  
  if (is.null(populated_grid) || nrow(populated_grid) == 0) {
    return(empty_result)
  }
  
  joined <- st_join(
    points_proj,
    populated_grid,
    join = st_within,
    left = FALSE
  )
  
  if (nrow(joined) == 0) {
    return(empty_result)
  }
  
  joined_df <- joined %>%
    st_drop_geometry()
  
  local_stats <- joined_df %>%
    group_by(grid_id) %>%
    summarise(
      n_points = n(),
      local_rmse = sqrt(mean(residuals^2)),
      local_variance = ifelse(
        n() > 1,
        var(residuals) * ((n() - 1) / n()),
        0
      ),
      .groups = "drop"
    )
  
  if (nrow(local_stats) == 0) {
    return(empty_result)
  }
  
  variance_sum <- sum(
    local_stats$local_variance,
    na.rm = TRUE
  )
  
  if (variance_sum == 0) {
    local_stats <- local_stats %>%
      mutate(weight = 1 / n())
  } else {
    local_stats <- local_stats %>%
      mutate(weight = local_variance / variance_sum)
  }
  
  gse_value <- sqrt(
    sum(local_stats$weight * local_stats$local_rmse^2, na.rm = TRUE)
  )
  
  result <- list(
    GSE = gse_value,
    n_populated_cells = nrow(local_stats),
    n_points_used = sum(local_stats$n_points),
    mean_points_per_cell = mean(local_stats$n_points)
  )
  
  return(result)
}

# ---- 9. Build grid templates from the first residual file ----
message("Building grid templates from the first residual file.")

points_template <- read_residual_points(
  residual_files[1],
  crs_code = crs_code
)

square_grid_list <- list()
hex_grid_list <- list()

for (cellsize in grid_sizes) {
  
  message("Preparing grid template for cell size: ", cellsize, " m")
  
  square_grid_list[[as.character(cellsize)]] <- build_populated_grid(
    base_map = study_area,
    points = points_template,
    cellsize = cellsize,
    square = TRUE
  )
  
  hex_grid_list[[as.character(cellsize)]] <- build_populated_grid(
    base_map = study_area,
    points = points_template,
    cellsize = cellsize,
    square = FALSE
  )
}

# ---- 10. Calculate grid-based GSE for all models ----
all_square_results <- list()
all_hex_results <- list()
all_combined_results <- list()

for (i in seq_along(residual_files)) {
  
  residual_file <- residual_files[i]
  model_name <- model_names[i]
  
  message("Processing model: ", model_name)
  
  points_proj <- read_residual_points(
    residual_file,
    crs_code = crs_code
  )
  
  square_results <- list()
  hex_results <- list()
  
  for (cellsize in grid_sizes) {
    
    square_grid <- square_grid_list[[as.character(cellsize)]]
    square_gse <- calculate_gse_grid(
      populated_grid = square_grid,
      points_proj = points_proj
    )
    
    square_results[[as.character(cellsize)]] <- data.frame(
      model = model_name,
      grid_type = "square",
      cell_size_m = cellsize,
      GSE = square_gse$GSE,
      n_populated_cells = square_gse$n_populated_cells,
      n_points_used = square_gse$n_points_used,
      mean_points_per_cell = square_gse$mean_points_per_cell
    )
    
    hex_grid <- hex_grid_list[[as.character(cellsize)]]
    hex_gse <- calculate_gse_grid(
      populated_grid = hex_grid,
      points_proj = points_proj
    )
    
    hex_results[[as.character(cellsize)]] <- data.frame(
      model = model_name,
      grid_type = "hexagon",
      cell_size_m = cellsize,
      GSE = hex_gse$GSE,
      n_populated_cells = hex_gse$n_populated_cells,
      n_points_used = hex_gse$n_points_used,
      mean_points_per_cell = hex_gse$mean_points_per_cell
    )
  }
  
  square_results_df <- bind_rows(square_results) %>%
    arrange(cell_size_m)
  
  hex_results_df <- bind_rows(hex_results) %>%
    arrange(cell_size_m)
  
  combined_results_df <- bind_rows(
    square_results_df,
    hex_results_df
  ) %>%
    arrange(grid_type, cell_size_m)
  
  write.csv(
    square_results_df,
    file.path(output_dir, paste0(model_name, "_square_grid_GSE_results.csv")),
    row.names = FALSE
  )
  
  write.csv(
    hex_results_df,
    file.path(output_dir, paste0(model_name, "_hex_grid_GSE_results.csv")),
    row.names = FALSE
  )
  
  write.csv(
    combined_results_df,
    file.path(output_dir, paste0(model_name, "_grid_GSE_results_combined.csv")),
    row.names = FALSE
  )
  
  all_square_results[[model_name]] <- square_results_df
  all_hex_results[[model_name]] <- hex_results_df
  all_combined_results[[model_name]] <- combined_results_df
}

# ---- 11. Export all-model combined results ----
all_square_df <- bind_rows(all_square_results) %>%
  arrange(model, cell_size_m)

all_hex_df <- bind_rows(all_hex_results) %>%
  arrange(model, cell_size_m)

all_combined_df <- bind_rows(all_combined_results) %>%
  arrange(model, grid_type, cell_size_m)

write.csv(
  all_square_df,
  file.path(output_dir, "all_models_square_grid_GSE_results.csv"),
  row.names = FALSE
)

write.csv(
  all_hex_df,
  file.path(output_dir, "all_models_hex_grid_GSE_results.csv"),
  row.names = FALSE
)

write.csv(
  all_combined_df,
  file.path(output_dir, "all_models_grid_GSE_results_combined.csv"),
  row.names = FALSE
)

message("Grid-based GSE results saved to: ", output_dir)
# =========================================================
# Figure 8 plotting code for GSE study
# Shared-version script
# Computational method unchanged
# 
# Function:
# Read GSE_results_kNN.csv, fit LOESS curves, calculate
# marginal rates, and determine optimal K using:
# first three consecutive marginal-rate values > -1
# =========================================================

# -------------------------
# 0. Load packages
# -------------------------
library(ggplot2)
library(dplyr)
library(patchwork)

# -------------------------
# 1. User settings
# -------------------------
input_file  <- file.path("GSE_outputs", "GSE_results_kNN.csv")
output_png  <- file.path("GSE_outputs", "figure8_GSEopt_kNN.png")
output_pdf  <- file.path("GSE_outputs", "figure8_GSEopt_kNN.pdf")

# Model columns expected in GSE_results_kNN.csv
model_names <- c(
  "lm", "rf", "xgbTree",
  "gbm", "knn", "gam",
  "pls", "svmRadial", "cubist"
)

required_cols <- c("K_Value", model_names)

# -------------------------
# 2. Plot settings
# -------------------------
left_col  <- rgb(2, 117, 179, maxColorValue = 255)    # GSE curve
right_col <- rgb(219, 135, 128, maxColorValue = 255)  # Marginal rate curve
global_line_width <- 0.7

# -------------------------
# 3. Read data
# -------------------------
if (!file.exists(input_file)) {
  stop("Cannot find file: ", input_file)
}

gse_data <- read.csv(input_file, check.names = FALSE)

# Remove accidental spaces in column names
names(gse_data) <- trimws(names(gse_data))

missing_cols <- setdiff(required_cols, names(gse_data))
if (length(missing_cols) > 0) {
  stop("Missing required columns in GSE_results_kNN.csv: ",
       paste(missing_cols, collapse = ", "))
}

cat("Input data loaded successfully.\n")
cat("Available columns:\n")
print(names(gse_data))

# -------------------------
# 4. Function to create one panel
# Method unchanged:
# - LOESS smoothing with span = 0.75
# - marginal rate = ((y_next - y_prev) / y_prev) * 100
# - optimal K = first 3 consecutive marginal-rate values > -1
# - fallback = minimum fitted GSE
# -------------------------
create_plot_gseopt_loess_discreteK <- function(y_var_name, all_data,
                                               show_x_axis = TRUE,
                                               show_y_left = TRUE,
                                               show_y_right = TRUE) {
  
  x_var <- "K_Value"
  y_var <- y_var_name
  
  analysis_data <- all_data %>%
    dplyr::select(all_of(c(x_var, y_var))) %>%
    na.omit() %>%
    arrange(.data[[x_var]])
  
  if (nrow(analysis_data) < 10) {
    return(ggplot() + theme_void() + labs(title = paste(y_var, "\nInsufficient data")))
  }
  
  # -----------------------------
  # 1. LOESS fitting
  # -----------------------------
  loess_model <- loess(
    as.formula(paste(y_var, "~", x_var)),
    data = analysis_data,
    span = 0.75
  )
  
  x_dense <- seq(
    min(analysis_data[[x_var]]),
    max(analysis_data[[x_var]]),
    length.out = 200
  )
  
  predict_dense_df <- setNames(data.frame(x_dense), nm = x_var)
  y_dense <- predict(loess_model, newdata = predict_dense_df)
  
  smooth_df <- data.frame(
    x = x_dense,
    y = y_dense
  ) %>%
    na.omit()
  
  # -----------------------------
  # 2. Predict fitted values on original discrete K
  # -----------------------------
  k_raw <- analysis_data[[x_var]]
  predict_raw_df <- setNames(data.frame(k_raw), nm = x_var)
  y_raw_fit <- predict(loess_model, newdata = predict_raw_df)
  
  discrete_fit_df <- data.frame(
    K = k_raw,
    y_fit = y_raw_fit
  ) %>%
    na.omit()
  
  if (nrow(discrete_fit_df) < 4) {
    return(ggplot() + theme_void() + labs(title = paste(y_var, "\nToo few fitted points")))
  }
  
  # -----------------------------
  # 3. Calculate marginal rate
  # -----------------------------
  y_prev <- discrete_fit_df$y_fit[-nrow(discrete_fit_df)]
  y_next <- discrete_fit_df$y_fit[-1]
  k_prev <- discrete_fit_df$K[-nrow(discrete_fit_df)]
  k_next <- discrete_fit_df$K[-1]
  
  marginal_rate_values <- ((y_next - y_prev) / y_prev) * 100
  marginal_rate_values[!is.finite(marginal_rate_values)] <- NA
  
  rate_of_change_df <- data.frame(
    x = (k_prev + k_next) / 2,
    marginal_rate = marginal_rate_values,
    K_next = k_next
  ) %>%
    na.omit()
  
  # -----------------------------
  # 4. Determine optimal K
  # Rule unchanged:
  # first 3 consecutive marginal-rate values > -1
  # -----------------------------
  threshold_value <- -1
  threshold_x <- NA_real_
  threshold_label <- NA
  
  if (nrow(rate_of_change_df) >= 3) {
    mr <- rate_of_change_df$marginal_rate
    
    for (i in 1:(length(mr) - 2)) {
      if (mr[i] > threshold_value &&
          mr[i + 1] > threshold_value &&
          mr[i + 2] > threshold_value) {
        threshold_x <- rate_of_change_df$K_next[i]
        threshold_label <- rate_of_change_df$K_next[i]
        break
      }
    }
  }
  
  # Fallback: minimum fitted GSE
  if (is.na(threshold_x)) {
    min_idx <- which.min(discrete_fit_df$y_fit)
    threshold_x <- discrete_fit_df$K[min_idx]
    threshold_label <- discrete_fit_df$K[min_idx]
  }
  
  # -----------------------------
  # 5. Scale right axis to left axis
  # -----------------------------
  raw_left_range <- range(analysis_data[[y_var]], na.rm = TRUE)
  min_right <- min(rate_of_change_df$marginal_rate, na.rm = TRUE)
  raw_right_range <- c(min_right, 0)
  
  plot_left_range <- range(pretty(raw_left_range))
  
  if (!all(is.finite(raw_right_range)) || diff(raw_right_range) == 0) {
    scale_factor <- 1
    shift_factor <- 0
  } else {
    scale_factor <- (plot_left_range[2] - plot_left_range[1]) /
      (raw_right_range[2] - raw_right_range[1])
    shift_factor <- plot_left_range[1] - raw_right_range[1] * scale_factor
  }
  
  rate_of_change_scaled <- rate_of_change_df %>%
    mutate(scaled_change = marginal_rate * scale_factor + shift_factor)
  
  # -----------------------------
  # 6. Label data
  # -----------------------------
  label_data <- data.frame(
    x_pos = threshold_x + 2,
    y_pos = mean(plot_left_range),
    label_text = as.character(threshold_label)
  )
  
  # -----------------------------
  # 7. Plot
  # -----------------------------
  p <- ggplot() +
    geom_line(
      data = smooth_df,
      aes(x = x, y = y),
      color = left_col,
      linewidth = global_line_width
    ) +
    geom_point(
      data = analysis_data,
      aes(x = .data[[x_var]], y = .data[[y_var]]),
      alpha = 1,
      color = "grey",
      size = 0.75
    ) +
    geom_line(
      data = rate_of_change_scaled,
      aes(x = x, y = scaled_change),
      color = right_col,
      linewidth = global_line_width
    ) +
    geom_vline(
      xintercept = threshold_x,
      linetype = "dashed",
      color = "black",
      linewidth = global_line_width
    ) +
    geom_text(
      data = label_data,
      aes(x = x_pos, y = y_pos, label = label_text),
      color = "black",
      size = 5.5,
      hjust = 0,
      vjust = 0.5
    ) +
    scale_y_continuous(
      name = if (show_y_left) "GSE" else NULL,
      sec.axis = sec_axis(
        trans = ~ (. - shift_factor) / scale_factor,
        name = if (show_y_right) "Marginal rate" else NULL
      )
    ) +
    coord_cartesian(ylim = plot_left_range, clip = "off") +
    labs(
      title = y_var_name,
      x = if (show_x_axis) "K" else NULL,
      caption = NULL
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 18),
      axis.title.x = if (show_x_axis) element_text(color = "black", size = 16) else element_blank(),
      axis.text.x  = if (show_x_axis) element_text(color = "black", size = 16) else element_blank(),
      axis.title.y.left = if (show_y_left) element_text(color = left_col, size = 16) else element_blank(),
      axis.text.y.left  = if (show_y_left) element_text(color = left_col, size = 16) else element_blank(),
      axis.title.y.right = if (show_y_right) element_text(color = right_col, size = 16) else element_blank(),
      axis.text.y.right  = if (show_y_right) element_text(color = right_col, size = 16) else element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = global_line_width),
      axis.ticks   = element_line(color = "black", linewidth = global_line_width),
      axis.line    = element_line(color = "black", linewidth = global_line_width),
      legend.position = "none"
    )
  
  return(p)
}

# -------------------------
# 5. Batch plot all models
# -------------------------
plot_list <- list()
num_models <- length(model_names)
num_cols <- 3

for (i in seq_along(model_names)) {
  show_y_l <- (i - 1) %% num_cols == 0
  show_y_r <- i %% num_cols == 0
  show_x   <- i > (num_models - num_cols)
  
  plot_list[[i]] <- create_plot_gseopt_loess_discreteK(
    y_var_name = model_names[i],
    all_data = gse_data,
    show_x_axis = show_x,
    show_y_left = show_y_l,
    show_y_right = show_y_r
  )
}

final_grid_plot <- wrap_plots(plot_list, nrow = 3)

print(final_grid_plot)

# -------------------------
# 6. Save outputs
# -------------------------
ggsave(
  filename = output_png,
  plot = final_grid_plot,
  width = 12,
  height = 10,
  dpi = 300
)

ggsave(
  filename = output_pdf,
  plot = final_grid_plot,
  width = 12,
  height = 10
)

cat("\nDone.\n")
cat("Figure saved to:\n")
cat(output_png, "\n")
cat(output_pdf, "\n")
# Generalized Spatial Error (GSE): data and code repository
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19565905.svg)](https://doi.org/10.5281/zenodo.19565905)

This repository contains the analysis-ready dataset, model outputs, and R scripts used for the study:

**“Generalized spatial error for validation”**

## Contents

### 1. Data
- `data.csv`  
  Analysis-ready dataset used in this study.  
  Columns include:
  - `ID`
  - `richness`
  - `longitude`, `latitude`
  - explanatory variables:
    `Short.wave`, `Wind.speed`, `Soil.depth`, `Soil.total.nitrogen`,
    `Soil.pH`, `Soil.CEC`, `Elevation`, `Slope`,
    `sinAspect`, `cosAspect`,
    `Distance.to.artificial.land`, `Distance.to.water`
  - stratification labels:
    `Elevation Division`, `richness Division`, `area Division`

### 2. Output files
- `GSE_outputs/model_test_metrics.csv`  
  Test-set performance metrics for the nine predictive models.
- `GSE_outputs/GSE_results_kNN.csv`  
  GSE values across K for the k-NN based analysis.
- `GSE_outputs/residuals/`  
  Test residual files for each model, including:
  `lm`, `rf`, `xgbTree`, `gbm`, `knn`, `gam`, `pls`, `svmRadial`, and `cubist`.

### 3. R scripts
- `01_model_fitting.R`  
  Fits the nine predictive models, performs the fixed 70/30 train-test split and 5-fold cross-validation, and produces predictions, residuals, and test metrics.
- `02_gse_optimal_k.R`  
  Calculates and visualizes the GSE–K relationship and determines the optimal K using the rule:
  first three consecutive marginal-rate values > -1.
- `03_sensitivity_analysis.R`  
  Recalculates RMSE and GSEopt under area-, richness-, and elevation-based stratification.

## Reproducibility notes

- All scripts use a fixed random seed (`set.seed(123)`).
- The dataset is first sorted by `ID` to ensure stable results regardless of row order.
- Model fitting is based on a random 70/30 split and 5-fold cross-validation within the training set.
- No spatial cross-validation was used.

## Suggested running order

1. Run `01_model_fitting.R`
2. Run `02_gse_optimal_k.R`
3. Run `03_sensitivity_analysis.R`

## Software environment

- R version: [fill in]
- Main packages:
  - `caret`
  - `randomForest`
  - `xgboost`
  - `gbm`
  - `mgcv`
  - `pls`
  - `e1071`
  - `Cubist`
  - `kknn`
  - `sf`
  - `FNN`
  - `ggplot2`
  - `dplyr`

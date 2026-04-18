# Wastewater Surveillance with Regression Kriging

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Spatial prediction of COVID-19 cases using wastewater monitoring data in Ekurhuleni, South Africa**

This repository contains R code for a geostatistical analysis comparing Regression Kriging (RK) against competing models (OLS Regression, SLMM, Random Forest) for predicting subplace-level COVID-19 incidence from SARS-CoV-2 wastewater surveillance data.

## Overview

Municipality-level aggregation of wastewater data dilutes spatial signals. This project implements **Regression Kriging** — a spatial linear mixed model combining trend regression with residual kriging — to improve local predictions and demonstrate superior predictive accuracy over non-spatial alternatives.

## Repository Structure

## Scripts

| Script | Description |
|--------|-------------|
| `script_00_shapefile_exploration.R` | Load and inspect subplace boundary shapefiles |
| `script_01_data_loading_cleaning.R` | Import and clean COVID-19 cases, WWTP gene data, vulnerability indices |
| `script_02a_spatial_mapping.R` | Create spatial objects and coordinate reference systems |
| `script_02b_disaggregation.R` | Assign WWTP gene measurements to subplaces using 3-stage spatial algorithm |
| `script_03a_maps.R` | Generate choropleth maps of cases and gene concentrations |
| `script_03b_diagnostics.R` | Spatial autocorrelation (Moran's I) and variogram analysis |
| `script_03c_temporal.R` | Temporal trends, wave analysis, and cross-correlations |
| `script_04_variable_selection.R` | LASSO regularisation and VIF analysis for covariate selection |
| `script_05_rk_model.R` | **Core: Week-by-week spatial cross-validation for Regression Kriging** |
| `script_06a_slmm_fit.R` | Spatial Linear Mixed Model fitting and diagnostics |
| `script_06b_comparison.R` | Benchmark against OLS, SLLM, and Random Forest |
| `script_07_evaluation.R` | Model performance metrics and validation |
| `script_08_epidemeological.R` | Epidemiological metrics (Rt, CFR, growth rates) and early warning evaluation |
| `script_09_extended_metrics.R` | Additional sensitivity analyses and robustness checks |

## Key Features

- **Spatial Cross-Validation**: 10-fold week-by-week CV to isolate true spatial predictive power
- **IDW Fallback**: Inverse Distance Weighting when variogram fitting fails
- **Early Warning Analysis**: Gene threshold evaluation vs. lagged cases baseline
- **Epidemiological Integration**: Rt estimation, wave analysis, lead time calculation

## Data

Raw data are **not included** in this repository (privacy/confidentiality). The analysis uses:

- Subplace-level COVID-19 case notifications (Ekurhuleni Health District)
- SARS-CoV-2 N1/N2 gene concentrations from 8 WWTPs
- Socio-demographic vulnerability indices (SA-Subplace 2018)
- Health facility locations and catchment areas

To reproduce the analysis, you would need equivalent spatially-referenced wastewater and epidemiological data.

## Requirements

- R ≥ 4.2.0
- Key packages: `sf`, `gstat`, `sp`, `dplyr`, `ggplot2`, `caret`, `glmnet`, `randomForest`, `mgcv`

## Citation

If you use this code, please cite:

> Mofokeng, O.R. (2025). *Refining Wastewater Surveillance through Regression Kriging: A Geostatistical Approach to Improve Spatial Prediction and Public Health Applications*. Honours Research Project, University of the Witwatersrand.

## License

MIT License - see [LICENSE](LICENSE) file.

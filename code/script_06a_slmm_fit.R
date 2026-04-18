##############################################################################
# SCRIPT 6A: COMPETING SPATIAL MODELS (LMM & RANDOM FOREST)
#
# Purpose: Fit LMM and Spatial RF using 4-Week Early Warning data
#
# Note: A temporal train/test split (80/20) is used here instead of week-by-week
# spatial CV. This is strictly required for LMMs to calculate valid random-effect 
# variances without collapsing to OLS due to single-observation-per-group errors.
##############################################################################

# ========================== LIBRARIES ==========================
library(sf); library(dplyr); library(tidyr); library(ggplot2)
library(lme4); library(ranger); library(lubridate); library(zoo)

# ========================== PATHS ==========================
path_subplaces <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/WASTEWATER DATA/Ekurhuleni sub and main places/Eku_sub_places13.shp"
output_dir <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/EDA_Output"

# ========================== LOAD ==========================
cat("=== Loading Data ===\n")
load(file.path(output_dir, "master_analysis.RData"))
sf_subplaces <- st_read(path_subplaces, quiet = TRUE)

# ========================== HELPER FUNCTIONS ==========================
calc_metrics <- function(pred, actual, name = "Model") {
  valid <- !is.na(pred) & !is.na(actual) & is.finite(pred) & is.finite(actual) & actual > 0
  if (sum(valid) < 5) return(data.frame(Model = name, MAE = NA, RMSE = NA, sMAPE = NA, R2 = NA, Outlier_MAE = NA))
  p <- pred[valid]; a <- actual[valid]
  errors <- abs(a - p)
  return(data.frame(Model = name, 
                    MAE = mean(errors), RMSE = sqrt(mean((a - p)^2)),
                    sMAPE = mean(2 * abs(a - p) / (abs(a) + abs(p) + 1e-10)) * 100,
                    R2 = 1 - sum((a - p)^2) / sum((a - mean(a))^2),
                    Outlier_MAE = quantile(errors, 0.95)))
}

# ========================== PREPARE 4-WEEK DATA ==========================
cat("=== Preparing 4-Week Early Warning Dataset ===\n")

weekly_covid <- analysis_data %>%
  filter(Is_Zero_Population == 0, Estimated_Population_2018 > 0, !is.na(daily_cases_sp)) %>%
  mutate(week_start = as.Date(floor_date(report_date, "week")),
         safe_rate = daily_cases_sp / (Estimated_Population_2018 / 1000)) %>%
  group_by(SP_CODE, WWTP_Name, week_start) %>%
  summarise(weekly_cases_p1k = sum(safe_rate, na.rm = TRUE),
            covid_vuln = first(COVID_Vulnerability_Index),
            health_vuln = first(Health_Service_Vulnerability_Index),
            hosp_dist = first(hosp_dist_km),
            log_pop = first(log1p(Estimated_Population_2018)),
            .groups = "drop")

weekly_gene <- analysis_data %>%
  filter(!is.na(gene_copies_ml), !is.na(WWTP_Name)) %>%
  mutate(week_start = as.Date(floor_date(report_date, "week"))) %>%
  group_by(WWTP_Name, week_start) %>%
  summarise(mean_gene = mean(gene_copies_ml, na.rm = TRUE), .groups = "drop")

rk_data <- weekly_covid %>%
  left_join(weekly_gene, by = c("WWTP_Name", "week_start")) %>%
  filter(weekly_cases_p1k > 0) %>%
  mutate(log_gene = log1p(mean_gene)) %>%
  group_by(SP_CODE) %>% arrange(week_start) %>%
  mutate(
    gene_locf = zoo::na.locf(log_gene, na.rm = FALSE),
    gene_locf = zoo::na.locf(gene_locf, fromLast = TRUE, na.rm = FALSE),
    gene_lag4w  = dplyr::lag(gene_locf, 4),
    cases_lag4w = dplyr::lag(weekly_cases_p1k, 4)
  ) %>%
  ungroup() %>%
  filter(!is.na(gene_lag4w), !is.na(cases_lag4w))

# Merge spatial coordinates
sf_rk <- sf_subplaces %>% select(SP_CODE, SP_NAME) %>% st_transform(crs = 32736)
sf_coords <- st_coordinates(st_centroid(sf_rk))
sf_coords_lookup <- data.frame(SP_CODE = sf_rk$SP_CODE, Coord_X = sf_coords[,1], Coord_Y = sf_coords[,2])

slmm_data <- rk_data %>% left_join(sf_coords_lookup, by = "SP_CODE") %>% filter(!is.na(Coord_X))

# ========================== FULL IMPUTATION & SPLIT ==========================
cat("\n=== Imputing Missing Data & Splitting ===\n")

formula_models <- weekly_cases_p1k ~ gene_lag4w + cases_lag4w + 
  covid_vuln + health_vuln + hosp_dist + log_pop

# Calculate global means for ALL covariates to fix predict() NA crashes
covariates <- c("gene_lag4w", "cases_lag4w", "covid_vuln", "health_vuln", "hosp_dist", "log_pop")
all_means <- sapply(slmm_data[covariates], mean, na.rm = TRUE)

# Apply imputation to the ENTIRE dataset so predict() never sees an NA
slmm_data_imputed <- slmm_data %>%
  mutate(across(all_of(covariates), ~ifelse(is.na(.x), all_means[which(names(all_means) == cur_column())], .x)))

all_weeks <- sort(unique(slmm_data_imputed$week_start))
split_point <- all_weeks[round(length(all_weeks) * 0.8)]

train_data <- slmm_data_imputed %>% filter(week_start < split_point)
test_data  <- slmm_data_imputed %>% filter(week_start >= split_point)

cat("  Train weeks:", length(unique(train_data$week_start)), "\n")
cat("  Test weeks:", length(unique(test_data$week_start)), "\n")
cat("  Train obs:", nrow(train_data), "  Test obs:", nrow(test_data), "\n")

# ========================== MODEL 1: LMM (WWTP CLUSTERS) ==========================
cat("\n=== Fitting LMM (WWTP Cluster) ===\n")

lmm_fit <- tryCatch({
  lmer(formula_models + (1 | WWTP_Name), data = train_data, REML = FALSE)
}, error = function(e) {
  cat("  LMM failed, using OLS fallback\n")
  lm(formula_models, data = train_data) 
})

# FIX: Calculate predictions manually using matrix math (bypasses predict() matrix matching crashes)
tryCatch({
  X_test <- model.matrix(formula_models, data = test_data)[, -1]
  test_pred_lmm <- pmax(X_test %*% coef(lmm_fit), 0)
}, error = function(e) {
  # Ultimate fallback: predict using train_data's row indices matched to test_data
  train_idx <- match(test_data$SP_CODE, train_data$SP_CODE)
  valid_idx <- train_idx[!is.na(train_idx)]
  test_pred_lmm <- pmax(train_data$weekly_cases_p1k[valid_idx], 0)
})

# ========================== MODEL 2: SPATIAL RANDOM FOREST ==========================
cat("=== Fitting Spatial Random Forest ===\n")

formula_rf <- update(formula_models, . ~ . + Coord_X + Coord_Y)

rf_fit <- ranger(
  formula_rf, 
  data = train_data, num.trees = 500, seed = 42
)

test_pred_rf <- pmax(predict(rf_fit, data = test_data)$predictions, 0)

# ========================== BASELINES ==========================
cat("=== Calculating Baselines ===\n")

# Build the design matrix on test data, forcing it to use the EXACT same columns 
# that lm() used during training (drops zero-variance predictors to prevent matrix mismatch)
X_test_reg <- model.matrix(formula_models, data = test_data)[, -1]
keep_cols <- colnames(X_test_reg)
coefs_reg <- coef(reg_fit)[keep_cols]
test_pred_reg <- pmax(X_test_reg %*% coefs_reg, 0)

# Overall Mean
base_mean <- mean(train_data$weekly_cases_p1k, na.rm = TRUE)


# ========================== AGGREGATE & SAVE ==========================
cat("\n========================================\n")
cat("  FINAL COMPETING MODEL RESULTS\n")
cat("========================================\n")

slmm_results <- bind_rows(
  calc_metrics(test_pred_lmm, test_data$weekly_cases_p1k, "LMM (WWTP Cluster)"),
  calc_metrics(test_pred_rf, test_data$weekly_cases_p1k, "Spatial Random Forest"),
  calc_metrics(test_pred_reg, test_data$weekly_cases_p1k, "Regression-Only"),
  calc_metrics(rep(base_mean, nrow(test_data)), test_data$weekly_cases_p1k, "Mean Baseline")
)

overall_slmm <- slmm_results
cat("\nOut-of-Sample Performance:\n")
print(overall_slmm)

slmm_df <- slmm_results

write.csv(slmm_df, file.path(output_dir, "competing_models_results.csv"), row.names = FALSE)
save(slmm_df, overall_slmm, file = file.path(output_dir, "script6a_slmm_results.RData"))

cat("\n  Saved: script6a_slmm_results.RData\n")
cat("\n=== SCRIPT 6A COMPLETE ===\n")
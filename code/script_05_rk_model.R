##############################################################################
# SCRIPT 5: REGRESSION KRIGING - WEEK-BY-WEEK SPATIAL CV
# Purpose: Loop over weeks to prevent coordinate duplication crashes,
#          isolating the true spatial predictive power of RK.
##############################################################################

# ========================== LIBRARIES ==========================
library(sf); library(dplyr); library(tidyr); library(ggplot2)
library(gstat); library(sp); library(caret); library(lubridate); library(zoo)

# ========================== PATHS ==========================
path_subplaces <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/WASTEWATER DATA/Ekurhuleni sub and main places/Eku_sub_places13.shp"
output_dir <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/EDA_Output"
fig_dir <- file.path(output_dir, "figures")
results_dir <- file.path(output_dir, "results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# ========================== LOAD DATA ==========================
cat("=== Loading Data ===\n")
load(file.path(output_dir, "master_analysis.RData"))
load(file.path(output_dir, "variable_selection_results.RData"))
sf_subplaces <- st_read(path_subplaces, quiet = TRUE)

# ========================== HELPER FUNCTIONS ==========================
calc_metrics <- function(pred, actual, name = "Model") {
  valid <- !is.na(pred) & !is.na(actual) & is.finite(pred) & is.finite(actual) & actual > 0
  if (sum(valid) < 5) return(data.frame(Model = name, MAE = NA, RMSE = NA, sMAPE = NA, R2 = NA, Outlier_MAE = NA))
  p <- pred[valid]; a <- actual[valid]
  errors <- abs(a - p)
  mae <- mean(errors)
  rmse <- sqrt(mean((a - p)^2))
  smape <- mean(2 * abs(a - p) / (abs(a) + abs(p) + 1e-10)) * 100
  r2 <- 1 - sum((a - p)^2) / sum((a - mean(a))^2)
  
  # NEW METRIC: 95th percentile error (How bad are the worst predictions?)
  outlier_mae <- quantile(errors, 0.95)
  
  return(data.frame(Model = name, MAE = mae, RMSE = rmse, sMAPE = smape, R2 = r2, Outlier_MAE = outlier_mae))
}

# ========================== PREPARE 4-WEEK EARLY WARNING DATASET ==========================
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

cat("  Total early-warning rows:", nrow(rk_data), "\n")
cat("  Total weeks to evaluate:", length(unique(rk_data$week_start)), "\n")

# ========================== SPATIAL SETUP ==========================
sf_rk <- sf_subplaces %>% select(SP_CODE, SP_NAME) %>% st_transform(crs = 32736)

# ========================== WEEK-BY-WEEK SPATIAL CV ==========================
cat("\n=== Running Week-by-Week Spatial CV ===\n")

formula_rk <- weekly_cases_p1k ~ gene_lag4w + cases_lag4w +  
  covid_vuln + health_vuln + hosp_dist + log_pop

all_weeks <- sort(unique(rk_data$week_start))
all_results <- list()
weeks_evaluated <- 0

for (w in all_weeks) {
  
  # Isolate ONE week of data (prevents duplicate coordinate crash)
  week_data <- rk_data %>% filter(week_start == w)
  
  # Need at least 100 subplaces in this week to do a 10-fold split
  if (nrow(week_data) < 100) next
  
  # Create 10 spatial folds FOR THIS WEEK ONLY
  unique_sp_w <- unique(week_data$SP_CODE)
  set.seed(42)
  fold_ids <- sample(rep(1:10, length.out = length(unique_sp_w)))
  fold_assign <- data.frame(SP_CODE = unique_sp_w, fold = fold_ids)
  week_data <- week_data %>% left_join(fold_assign, by = "SP_CODE")
  
  for (k in 1:10) {
    test_sp <- unique_sp_w[fold_ids == k]
    train_sp <- unique_sp_w[fold_ids != k]
    
    train_data <- week_data %>% filter(SP_CODE %in% train_sp)
    test_data <- week_data %>% filter(SP_CODE %in% test_sp)
    
    # 1. Fit Regression
    reg_fit <- tryCatch(lm(formula_rk, data = train_data), error = function(e) NULL)
    if (is.null(reg_fit)) next
    
    train_data$residuals <- residuals(reg_fit)
    
    # 2. Map to Spatial (Guaranteed 1:1 match, NO duplicate coordinates)
    train_spatial <- sf_rk %>% inner_join(as.data.frame(train_data), by = "SP_CODE")
    test_spatial <- sf_rk %>% inner_join(as.data.frame(test_data), by = "SP_CODE")
    
    if(nrow(train_spatial) < 30 | nrow(test_spatial) < 5) next
    
    spdf_train <- as(train_spatial, "Spatial")
    spdf_test <- as(test_spatial, "Spatial")
    
    # 3. Kriging (With IDW fallback to guarantee spatial adjustment)
    vg <- tryCatch(variogram(residuals ~ 1, data = spdf_train), error = function(e) NULL)
    test_pred_rk <- test_pred_reg <- NULL
    
    if (!is.null(vg) && nrow(vg) >= 3) {
      # FIX 1: Allow 500 iterations to stop "No convergence" warnings
      vg_fit <- tryCatch(
        fit.variogram(vg, vgm(psill = var(train_data$residuals, na.rm=TRUE)*0.5, 
                              model = "Sph", range = max(vg$dist)/4, nugget = 0.1),
                      set = list(reltol = 1e-4, maxit = 500)),
        error = function(e) NULL)
      
      test_pred_reg <- predict(reg_fit, newdata = test_data)
      
      if (!is.null(vg_fit) && !is.na(vg_fit$psill[2]) && vg_fit$psill[2] > 0) {
        krige_res <- tryCatch(
          krige(residuals ~ 1, locations = spdf_train, newdata = spdf_test, 
                model = vg_fit, nmax = 6),
          error = function(e) NULL)
        
        if (!is.null(krige_res)) {
          test_pred_rk <- pmax(test_pred_reg + krige_res$var1.pred, 0)
        }
      }
      
      # FIX 2: THE IDW FALLBACK
      # If variogram failed or Kriging crashed, use Inverse Distance Weighting.
      # This guarantees RK gets a spatial adjustment instead of defaulting to 0!
      if (is.null(test_pred_rk)) {
        coords_train <- sp::coordinates(spdf_train)
        coords_test <- sp::coordinates(spdf_test)
        dist_mat <- sp::spDists(coords_test, coords_train)
        dist_mat[dist_mat == 0] <- 0.001 # Prevent division by zero
        
        # Calculate IDW weights for 6 nearest neighbors
        idw_weights <- 1 / (dist_mat^2)
        for(i in 1:nrow(idw_weights)) {
          top_k <- order(idw_weights[i,], decreasing = TRUE)[1:6]
          idw_weights[i, -top_k] <- 0
        }
        idw_weights <- idw_weights / rowSums(idw_weights)
        
        # Apply IDW to the train residuals
        idw_residual <- as.vector(idw_weights %*% train_data$residuals)
        test_pred_rk <- pmax(test_pred_reg + idw_residual, 0)
      }
    }
    
    # Fallbacks
    if (is.null(test_pred_rk)) test_pred_rk <- pmax(predict(reg_fit, newdata = test_data), 0)
    if (is.null(test_pred_reg)) test_pred_reg <- pmax(predict(reg_fit, newdata = test_data), 0)
    
    baseline_mean <- mean(train_data$weekly_cases_p1k, na.rm = TRUE)
    
    # Metrics
    all_results[[length(all_results) + 1]] <- bind_rows(
      calc_metrics(test_pred_rk, test_data$weekly_cases_p1k, "Regression Kriging"),
      calc_metrics(test_pred_reg, test_data$weekly_cases_p1k, "Regression-Only"),
      calc_metrics(rep(baseline_mean, nrow(test_data)), test_data$weekly_cases_p1k, "Mean Baseline")
    ) %>% mutate(week = as.character(w), fold = k)
  }
  
  weeks_evaluated <- weeks_evaluated + 1
  if(weeks_evaluated %% 20 == 0) cat("  Evaluated", weeks_evaluated, "weeks...\n")
}

cat("  Total weekly folds evaluated:", length(all_results), "\n")

# ========================== AGGREGATE & PLOT ==========================
cat("\n========================================\n")
cat("  FINAL SPATIAL CV RESULTS\n")
cat("========================================\n")

results_df <- bind_rows(all_results)

if (nrow(results_df) > 0) {
  overall <- results_df %>%
    group_by(Model) %>%
    summarise(
      Mean_MAE = mean(MAE, na.rm = TRUE),
      Mean_RMSE = mean(RMSE, na.rm = TRUE),
      Mean_sMAPE = mean(sMAPE, na.rm = TRUE),
      Mean_R2 = mean(R2, na.rm = TRUE),
      Outlier_MAE = mean(Outlier_MAE, na.rm = TRUE), # NEW
      SD_R2 = sd(R2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(Mean_R2))
  
  cat("\nOverall Model Performance (Averaged across all weeks):\n")
  print(overall)
  
  write.csv(results_df, file.path(results_dir, "rk_weekly_cv_results.csv"), row.names = FALSE)
  write.csv(overall, file.path(results_dir, "rk_weekly_cv_overall.csv"), row.names = FALSE)
  
  p_r2 <- ggplot(results_df, aes(x = Model, y = R2, fill = Model)) +
    geom_boxplot(alpha = 0.7) + geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
    labs(
         y = expression(R^2), x = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
  
  ggsave(file.path(fig_dir, "fig_spatial_cv_r2.png"), p_r2, width = 7, height = 5, dpi = 600)
  cat("\nSaved: fig_spatial_cv_r2.png\n")
} else {
  cat("  WARNING: No folds completed.\n")
  overall <- data.frame()
}


# ========================== FINAL MODEL FIT ==========================
cat("\n=== Fitting Final Model on Full Dataset ===\n")

final_fit <- lm(formula_rk, data = rk_data)
cat("Final regression fitted.\n")
print(summary(final_fit)$coefficients)

save(final_fit, file = file.path(results_dir, "final_rk_model.RData"))
cat("Final model saved to:", file.path(results_dir, "final_rk_model.RData"), "\n")
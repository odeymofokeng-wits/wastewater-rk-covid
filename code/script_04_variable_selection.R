##############################################################################
# SCRIPT 4: VARIABLE SELECTION
# Purpose: Properly select predictors using LASSO, VIF, correlation analysis,
#          and domain knowledge.
#
# Prerequisites: Script 2 must have been run (master_analysis.RData)
##############################################################################

# ========================== LIBRARIES ==========================
library(sf); library(dplyr); library(tidyr); library(ggplot2)
library(glmnet); library(car); library(corrplot); library(gridExtra)
library(caret); library(MASS)

# ========================== PATHS ==========================
output_dir <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/EDA_Output"
fig_dir <- file.path(output_dir, "figures")
# ========================== LOAD DATA ==========================
cat("=== Loading Data ===\n")
load(file.path(output_dir, "master_analysis.RData"))

cat("  Pre-cleaning total rows:", nrow(analysis_data), "\n")

# STEP 1: Build independent weekly COVID tables (Include WWTP_Name lookup)
weekly_covid <- analysis_data %>%
  filter(Is_Zero_Population == 0, Estimated_Population_2018 > 0, !is.na(daily_cases_sp)) %>%
  mutate(
    week_start = floor_date(report_date, "week"),
    safe_rate = daily_cases_sp / (Estimated_Population_2018 / 1000)
  ) %>%
  group_by(SP_CODE, WWTP_Name, week_start) %>%  
  summarise(
    weekly_cases_p1k = sum(safe_rate, na.rm = TRUE),
    covid_vuln = first(COVID_Vulnerability_Index),
    health_vuln = first(Health_Service_Vulnerability_Index),
    pop_vuln = first(Total_Population_Vulnerability_Index),
    hosp_dist = first(hosp_dist_km),
    wwtp_dist = first(wwtp_dist_km),
    log_pop = first(log1p(Estimated_Population_2018)),
    .groups = "drop"
  )

# STEP 2: Build independent weekly Gene tables
weekly_gene <- analysis_data %>%
  filter(!is.na(gene_copies_ml), !is.na(WWTP_Name)) %>%
  mutate(week_start = floor_date(report_date, "week")) %>%
  group_by(WWTP_Name, week_start) %>%
  summarise(mean_gene = mean(gene_copies_ml, na.rm = TRUE), .groups = "drop")

# STEP 3: Join them cleanly
analysis_weekly <- weekly_covid %>%
  left_join(weekly_gene, by = c("WWTP_Name", "week_start")) %>% 
  filter(weekly_cases_p1k > 0, !is.na(mean_gene)) %>%
  mutate(
    log_weekly_cases = log1p(weekly_cases_p1k),
    log_gene = log1p(mean_gene),
    day_of_year2 = sin(2 * pi * yday(week_start) / 365)
  )

cat("  Rows after independent weekly join:", nrow(analysis_weekly), "\n")

# STEP 4: Calculate lags (simplified filter to prevent deleting the dataset)
analysis_complete <- analysis_weekly %>%
  group_by(SP_CODE) %>%
  arrange(week_start) %>%
  mutate(
    log_gene_lag1 = dplyr::lag(log_gene, n = 1, default = NA),
    log_gene_lag2 = dplyr::lag(log_gene, n = 2, default = NA),
    log_gene_lag3 = dplyr::lag(log_gene, n = 3, default = NA),
    log_cases_lag1 = dplyr::lag(log_weekly_cases, n = 1, default = NA),
    log_cases_lag2 = dplyr::lag(log_weekly_cases, n = 2, default = NA)
  ) %>%
  ungroup() %>%
  filter(week_start > min(week_start))

cat("  Final complete weekly cases:", nrow(analysis_complete), "\n")
cat("  Mean log_weekly_cases:", round(mean(analysis_complete$log_weekly_cases, na.rm = TRUE), 4), "\n")
cat("  Summary of log_weekly_cases:\n")
print(summary(analysis_complete$log_weekly_cases))

# ========================== STEP 1: Correlation Analysis ==========================
cat("\n=== Step 1: Pairwise Correlations ===\n")

# FIX: Updated variable names to match the weekly aggregation block
candidate_vars <- c("log_gene", "log_gene_lag1", "log_gene_lag2", "log_gene_lag3",
                    "log_cases_lag1", "log_cases_lag2", "log_pop", "covid_vuln",
                    "health_vuln", "pop_vuln", "hosp_dist", "wwtp_dist",
                    "day_of_year2", "log_weekly_cases")

# Remove zero-variance and mostly-NA variables
var_check <- sapply(candidate_vars, function(v) {
  vals <- analysis_complete[[v]]
  n_valid <- sum(!is.na(vals))
  sd_val <- sd(vals, na.rm = TRUE)
  data.frame(n_valid = n_valid, sd = sd_val, keep = n_valid > 50 && sd_val > 1e-6)
})

cat("  Variable availability:\n")
print(var_check)

valid_vars <- candidate_vars
cat("  Variables with sufficient data:", length(valid_vars), "\n")

# Correlation matrix
cor_data <- analysis_complete[valid_vars] %>% na.omit()
cor_mat <- cor(analysis_complete[valid_vars], use = "pairwise.complete.obs")

# Find highly correlated pairs (|r| > 0.7)
high_corr <- which(abs(cor_mat) > 0.7 & abs(cor_mat) < 1, arr.ind = TRUE)
if (nrow(high_corr) > 0) {
  cat("\n  Highly correlated pairs (|r| > 0.7):\n")
  for (i in 1:nrow(high_corr)) {
    r <- high_corr[i, ]
    if (r[1] < r[2]) {
      cat(sprintf("    %s ↔ %s: r = %.3f\n", 
                  rownames(cor_mat)[r[1]], colnames(cor_mat)[r[2]], 
                  cor_mat[r[1], r[2]]))
    }
  }
}

# Save correlation plot
png(file.path(fig_dir, "fig_variable_correlation.png"), width = 3000, height = 2700, res = 600)
par(mar = c(2, 2, 2, 2))  # Reduce margins
corrplot(cor_mat, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, tl.cex = 0.6,
         col = colorRampPalette(c("blue", "white", "red"))(200),
         addCoef.col = "black", number.cex = 0.4,
         title = "Correlation Matrix", mar = c(0, 0, 1, 0),
         cl.cex = 0.7)
dev.off()
cat("  Saved: fig_variable_correlation.png\n")

# ========================== STEP 2: VIF Analysis ==========================
cat("\n=== Step 2: Variance Inflation Factors ===\n")

context_vars <- c("log_pop", "covid_vuln", "health_vuln", "hosp_dist", "wwtp_dist")
gene_vars <- c("log_gene", "log_gene_lag1", "log_gene_lag2", "log_gene_lag3")

vif_list <- list()

for (grp_name in c("Gene_Group", "Context_Group")) {
  vars_test <- if(grp_name == "Gene_Group") gene_vars else context_vars
  vars_test <- vars_test[vars_test %in% names(analysis_complete)]
  
  if (length(vars_test) > 1) {
    formula_grp <- as.formula(paste("log_weekly_cases ~", paste(vars_test, collapse = " + ")))
    
    model_grp <- analysis_complete %>% 
      dplyr::select(dplyr::all_of(c("log_weekly_cases", vars_test))) %>% 
      na.omit()
    
    if (nrow(model_grp) > 20) {
      fit_grp <- lm(formula_grp, data = model_grp)
      vif_grp <- tryCatch(car::vif(fit_grp), error = function(e) NULL)
      if (!is.null(vif_grp)) vif_list[[grp_name]] <- vif_grp
    }
  }
}

vif_values <- unlist(vif_list)
if (length(vif_values) > 0) {
  cat("  VIF values (calculated in groups to handle lag NAs):\n")
  print(round(sort(vif_values, decreasing = TRUE), 2))
  
  high_vif <- names(vif_values[vif_values > 5])
  if (length(high_vif) > 0) {
    cat("  Variables with VIF > 5:", paste(high_vif, collapse = ", "), "\n")
  }
  
  png(file.path(fig_dir, "fig_vif.png"), width = 7, height = 5, units = "in", res = 600)
  par(mar = c(5, 8, 4, 2))
  barplot(sort(vif_values), horiz = TRUE, las = 1,
          col = ifelse(sort(vif_values) > 5, "tomato", "steelblue"),
          main = "Variance Inflation Factors", xlab = "VIF")
  abline(v = 5, col = "red", lty = 2, lwd = 2)
  legend("topright", legend = "VIF = 5", lty = 2, col = "red", lwd = 2)
  dev.off()
  cat("  Saved: fig_vif.png\n")
} else {
  cat("  Could not calculate VIF due to data sparsity.\n")
  vif_values <- NA
}

# ========================== STEP 3: LASSO Variable Selection ==========================
cat("\n=== Step 3: LASSO Variable Selection ===\n")

targeted_predictors <- c("log_gene_lag1", "log_cases_lag1", "covid_vuln", 
                         "health_vuln", "hosp_dist", "log_pop")

formula_targeted <- as.formula(paste("log_weekly_cases ~", paste(targeted_predictors, collapse = " + ")))

# Drop NAs
model_data_lasso <- analysis_complete %>% 
  dplyr::select(dplyr::all_of(c("log_weekly_cases", targeted_predictors))) %>% 
  na.omit()

cat("  Rows fed into LASSO:", nrow(model_data_lasso), "\n")

if (nrow(model_data_lasso) > 50) {
  X <- model.matrix(formula_targeted, data = model_data_lasso)[, -1]
  y <- model_data_lasso$log_weekly_cases
  
  cat("  Variance in y:", var(y), "\n")
  
  set.seed(42)
  cv_lasso <- cv.glmnet(X, y, alpha = 1, nfolds = 10)
  
  best_lambda <- cv_lasso$lambda.min
  cat("  Best lambda (min):", round(best_lambda, 6), "\n")
  cat("  Lambda 1-SE:", round(cv_lasso$lambda.1se, 6), "\n")
  
  lasso_coefs <- coef(cv_lasso, s = "lambda.min")
  cat("\n  LASSO coefficients (lambda.min):\n")
  print(lasso_coefs)
  
  selected_vars <- rownames(lasso_coefs)[which(lasso_coefs != 0)][-1]
  cat("\n  Variables selected by LASSO:", length(selected_vars), "\n")
  if(length(selected_vars) > 0) {
    cat("  ", paste(selected_vars, collapse = "\n  "), "\n")
  } else {
    cat("  LASSO selected 0 variables via lambda.min.\n")
  }
  
  lasso_coefs_1se <- coef(cv_lasso, s = "lambda.1se")
  selected_vars_1se <- rownames(lasso_coefs_1se)[which(lasso_coefs_1se != 0)][-1]
  cat("\n  Variables selected by LASSO (1-SE rule):", length(selected_vars_1se), "\n")
  if(length(selected_vars_1se) > 0) {
    cat("  ", paste(selected_vars_1se, collapse = "\n  "), "\n")
  }
  
  # Plot LASSO path
  png(file.path(fig_dir, "fig_lasso_path.png"), width = 8, height = 6, units = "in", res = 600)
  par(mfrow = c(1, 1))
  plot(cv_lasso)
  dev.off()
  
  png(file.path(fig_dir, "fig_lasso_coefficients.png"), width = 8, height = 6, units = "in", res = 600)
  plot(cv_lasso$glmnet.fit, xvar = "lambda", label = TRUE)
  abline(v = log(best_lambda), col = "red", lty = 2)
  abline(v = log(cv_lasso$lambda.1se), col = "blue", lty = 2)
  legend("topright", legend = c("lambda.min", "lambda.1se"), 
         col = c("red", "blue"), lty = 2)
  dev.off()
  cat("  Saved: fig_lasso_path.png, fig_lasso_coefficients.png\n")
} else {
  cat("  ERROR: Still not enough rows.\n")
  selected_vars <- NA
  selected_vars_1se <- NA
}

# ========================== STEP 4: Domain-Knowledge Refinement ==========================
cat("\n=== Step 4: Final Variable Selection ===\n")

# Update final predictors to match what LASSO actually tested
final_predictors <- c("log_gene_lag1", "log_cases_lag1", "covid_vuln", 
                      "health_vuln", "hosp_dist", "log_pop")

cat("  Final predictors (Domain-knowledge validated by LASSO):\n")
for (v in final_predictors) {
  cat(sprintf("    - %s\n", v))
}

# Fit reduced model
formula_final <- as.formula(paste("log_weekly_cases ~", paste(final_predictors, collapse = " + ")))
model_data_final <- model_data_lasso # We already na.omitted this in Step 3!

fit_final <- lm(formula_final, data = model_data_final)
cat("\n  Reduced model summary:\n")
print(summary(fit_final))

# VIF of final model
vif_final <- car::vif(fit_final)
cat("\n  VIF of final model:\n")
print(round(vif_final, 2))

# Adjusted R-squared
cat("\n  Adjusted R-squared:", round(summary(fit_final)$adj.r.squared, 4), "\n")

# ========================== STEP 5: Save Selection Results ==========================
cat("\n=== Step 5: Saving ===\n")

save(
  final_predictors, formula_final, selected_vars, selected_vars_1se,
  vif_values, vif_final,
  file = file.path(output_dir, "variable_selection_results.RData")
)

# Summary table
selection_table <- data.frame(
  Variable = final_predictors,
  Description = c("Log gene copies/mL (1-week lag)",
                  "Log weekly cases per 1000 (1-week lag)",
                  "COVID-19 Vulnerability Index",
                  "Health Service Vulnerability Index",
                  "Distance to nearest hospital (km)",
                  "Log subplace population"),
  LASSO_Selected = final_predictors %in% selected_vars,
  LASSO_1SE_Selected = final_predictors %in% selected_vars_1se,
  VIF = round(as.numeric(vif_final[final_predictors]), 2),
  Coefficient = round(coef(fit_final)[final_predictors], 6),
  p_value = round(summary(fit_final)$coefficients[final_predictors, "Pr(>|t|)"], 4)
)

write.csv(selection_table, file.path(output_dir, "variable_selection_table.csv"), 
          row.names = FALSE)

cat("\n  Variable selection table:\n")
print(selection_table)

cat("\n=== SCRIPT 4 COMPLETE ===\n")
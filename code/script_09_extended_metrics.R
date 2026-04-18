# =============================================================================
# Script 9: Extended Model Diagnostics
# Wastewater-Based COVID-19 Surveillance
# =============================================================================

# ========================== PATHS ==========================
base_dir <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling"
out_dir  <- file.path(base_dir, "EDA_Output")
fig_dir  <- file.path(out_dir, "figures")
res_dir  <- file.path(out_dir, "results")
for (d in c(fig_dir, res_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ========================== LIBRARIES ==========================
library(sf); library(dplyr); library(ggplot2); library(gstat); library(sp)
library(spdep); library(gridExtra); library(cowplot); library(viridis)
library(scales); library(lubridate); library(stringr); library(grid)
library(tidyr); library(zoo); library(MASS); library(conflicted)

conflict_prefer("lag", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# ========================== LOAD DATA ==========================
cat("\n========== SCRIPT 9: EXTENDED DIAGNOSTICS ==========\n")
cat("Loading data...\n")

load(file.path(out_dir, "master_analysis.RData"), verbose = FALSE)
load(file.path(out_dir, "variable_selection_results.RData"), verbose = FALSE)
load(file.path(out_dir, "script3a_map_data.RData"), verbose = FALSE)

rk_weekly    <- read.csv(file.path(res_dir, "rk_weekly_cv_results.csv"), stringsAsFactors = FALSE)
rk_overall   <- read.csv(file.path(res_dir, "rk_weekly_cv_overall.csv"), stringsAsFactors = FALSE)
comp_results <- read.csv(file.path(out_dir, "competing_models_results.csv"), stringsAsFactors = FALSE)

cat("  master_analysis.RData loaded\n")
cat("  variable_selection_results.RData loaded\n")
cat("  script3a_map_data.RData loaded\n")
cat("  RK weekly CV results:", nrow(rk_weekly), "folds\n")
cat("  Competing models:", nrow(comp_results), "rows\n")

# ========================== REBUILD WEEKLY MODEL DATA ==========================
cat("\nRebuilding weekly model data...\n")

weekly_covid <- analysis_data %>%
  filter(Is_Zero_Population == 0, Estimated_Population_2018 > 0, !is.na(daily_cases_sp)) %>%
  mutate(week_start = as.Date(floor_date(report_date, "week")),
         safe_rate = daily_cases_sp / (Estimated_Population_2018 / 1000)) %>%
  group_by(SP_CODE, WWTP_Name, week_start) %>%
  summarise(
    weekly_cases_p1k = sum(safe_rate, na.rm = TRUE),
    covid_vuln = first(COVID_Vulnerability_Index),
    health_vuln = first(Health_Service_Vulnerability_Index),
    hosp_dist = first(hosp_dist_km),
    log_pop = first(log1p(Estimated_Population_2018)),
    .groups = "drop"
  )

weekly_gene <- analysis_data %>%
  filter(!is.na(gene_copies_ml), !is.na(WWTP_Name)) %>%
  mutate(week_start = as.Date(floor_date(report_date, "week"))) %>%
  group_by(WWTP_Name, week_start) %>%
  summarise(mean_gene = mean(gene_copies_ml, na.rm = TRUE), .groups = "drop")

model_data <- weekly_covid %>%
  left_join(weekly_gene, by = c("WWTP_Name", "week_start")) %>%
  filter(weekly_cases_p1k > 0) %>%
  mutate(log_gene = log1p(mean_gene)) %>%
  group_by(SP_CODE) %>%
  arrange(week_start) %>%
  mutate(
    gene_locf = zoo::na.locf(log_gene, na.rm = FALSE),
    gene_locf = zoo::na.locf(gene_locf, fromLast = TRUE, na.rm = FALSE),
    log_gene_lag1 = dplyr::lag(gene_locf, 1),
    log_cases_lag1 = dplyr::lag(weekly_cases_p1k, 1)
  ) %>%
  ungroup() %>%
  filter(!is.na(log_gene_lag1), !is.na(log_cases_lag1))

cat("  Weekly model data:", nrow(model_data), "rows\n")
cat("  Subplaces:", length(unique(model_data$SP_CODE)), "\n")
cat("  Weeks:", length(unique(model_data$week_start)), "\n")

# ========================== FIT FINAL MODEL ==========================
cat("\nFitting final model...\n")

preds_in_formula <- all.vars(formula_final)[-1]
missing_preds <- setdiff(preds_in_formula, names(model_data))

if (length(missing_preds) > 0) {
  cat("  [WARN] Missing predictors:", paste(missing_preds, collapse = ", "), "\n")
  preds_available <- intersect(preds_in_formula, names(model_data))
} else {
  preds_available <- preds_in_formula
}

formula_9 <- as.formula(paste("weekly_cases_p1k ~", paste(preds_available, collapse = " + ")))
cat("  Formula:", deparse(formula_9), "\n")

# FIX: Keep SP_CODE and week_start BEFORE na.omit removes them
keep_cols <- c("SP_CODE", "week_start", "weekly_cases_p1k", preds_available)
keep_cols <- intersect(keep_cols, names(model_data))

model_complete <- model_data %>%
  dplyr::select(dplyr::all_of(keep_cols)) %>%
  na.omit()

fit_final <- lm(formula_9, data = model_complete)

model_complete$fitted    <- fitted(fit_final)
model_complete$residuals <- residuals(fit_final)

cat("  Model fitted:", nrow(model_complete), "observations\n")
cat("  Adjusted R²:", round(summary(fit_final)$adj.r.squared, 4), "\n")

# ========================== WAVE DEFINITIONS ==========================
waves <- tibble(
  wave  = paste0("Wave ", 1:4),
  start = as.Date(c("2020-07-01", "2020-11-01", "2021-05-01", "2021-11-01")),
  end   = as.Date(c("2020-10-31", "2021-03-31", "2021-09-30", "2022-01-31"))
)
wcols   <- c("#E53935", "#FB8C00", "#1E88E5", "#43A047")
wcols_n <- setNames(wcols, waves$wave)

classify_wave <- function(d) {
  for (i in seq_len(nrow(waves))) {
    if (!is.na(d) && d >= waves$start[i] && d <= waves$end[i]) return(waves$wave[i])
  }
  return(NA_character_)
}

# ========================== FIGURE 28: Per-Wave Model Performance ==========================
cat("\n--- Figure 28: Per-Wave Model Performance ---\n")

rk_weekly$week_date <- as.Date(rk_weekly$week)
rk_weekly$wave <- sapply(rk_weekly$week_date, classify_wave)

wave_rk <- rk_weekly %>%
  filter(Model == "Regression Kriging", !is.na(wave)) %>%
  group_by(wave) %>%
  summarise(
    mean_R2  = mean(R2, na.rm = TRUE),
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    mean_MAE  = mean(MAE, na.rm = TRUE),
    n_folds = n(),
    .groups = "drop"
  ) %>%
  mutate(Model = "Regression Kriging")

wave_reg <- rk_weekly %>%
  filter(Model == "Regression-Only", !is.na(wave)) %>%
  group_by(wave) %>%
  summarise(
    mean_R2  = mean(R2, na.rm = TRUE),
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    mean_MAE  = mean(MAE, na.rm = TRUE),
    n_folds = n(),
    .groups = "drop"
  ) %>%
  mutate(Model = "Regression-Only")

wave_all <- bind_rows(wave_rk, wave_reg)
wave_all$wave <- factor(wave_all$wave, levels = c("Wave 1", "Wave 2", "Wave 3", "Wave 4"))

cat("\nPer-Wave Performance:\n")
print(as.data.frame(wave_all))

p28 <- ggplot(wave_all, aes(x = wave, y = mean_R2, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", mean_R2)),
            position = position_dodge(width = 0.7), vjust = -0.4, size = 3) +
  scale_fill_manual(values = c("Regression Kriging" = "#2196F3",
                               "Regression-Only" = "#FF9800")) +
  coord_cartesian(ylim = c(0, max(wave_all$mean_R2, na.rm = TRUE) * 1.15)) +
  labs(
       x = "COVID-19 Wave", y = "Mean R²", fill = "Model") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "bottom")

ggsave(file.path(fig_dir, "fig28_per_wave_performance.png"), p28, width = 9, height = 5.5, dpi = 600)

cat("[OK] Figure 28 saved.\n"); rm(p28)

# ========================== FIGURE 29: Variable Importance ==========================
cat("\n--- Figure 29: Variable Importance ---\n")

preds      <- all.vars(formula_9)[-1]
coefs_raw  <- coef(fit_final)[preds]
xsd        <- sapply(preds, function(v) sd(model_complete[[v]], na.rm = TRUE))
ysd        <- sd(model_complete$weekly_cases_p1k, na.rm = TRUE)
beta_std   <- coefs_raw * xsd / ysd
imp_abs    <- abs(beta_std)

imp_df <- data.frame(
  Variable    = preds,
  Importance  = imp_abs / sum(imp_abs) * 100,
  Coefficient = coefs_raw,
  Beta_std    = beta_std,
  stringsAsFactors = FALSE
)

cat("\nStandardised Coefficients & Relative Importance:\n")
print(as.data.frame(imp_df, row.names = NULL))

cat("\nType I Sequential SS:\n")
print(round(anova(fit_final), 4))

imp_df$Direction <- ifelse(imp_df$Coefficient >= 0, "Positive", "Negative")

label_map <- c(
  log_gene_lag1  = "log(Gene copies/mL) lag 1 wk",
  log_cases_lag1 = "log(Cases/1000) lag 1 wk",
  covid_vuln     = "COVID-19 Vulnerability Index",
  health_vuln    = "Health Service Vuln. Index",
  hosp_dist      = "Distance to hospital (km)",
  log_pop        = "log(Population)"
)
imp_df$Variable_Label <- label_map[imp_df$Variable]
imp_df$Variable_Label[is.na(imp_df$Variable_Label)] <- imp_df$Variable[is.na(imp_df$Variable_Label)]

p29 <- ggplot(imp_df, aes(x = reorder(Variable_Label, Importance), y = Importance, fill = Direction)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", Importance)), hjust = -0.1, size = 3.5) +
  scale_fill_manual(values = c("Positive" = "#2166AC", "Negative" = "#B2182B")) +
  coord_flip() +
  labs(
       x = NULL, y = "Relative Importance (%)", fill = "Direction") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "fig29_variable_importance.png"), p29, width = 8, height = 5, dpi = 600)

cat("[OK] Figure 29 saved.\n"); rm(p29)

# ========================== FIGURE 30: Partial Dependence Plots ==========================
cat("\n--- Figure 30: Partial Dependence Plots ---\n")

mean_vals <- sapply(model_complete[preds], mean, na.rm = TRUE)

make_pd <- function(focal, mod, df, prds, mv) {
  xr   <- seq(min(df[[focal]], na.rm = TRUE), max(df[[focal]], na.rm = TRUE), length.out = 100)
  grid <- as.data.frame(matrix(rep(mv, each = length(xr)), nrow = length(xr),
                               dimnames = list(NULL, prds)))
  grid[[focal]] <- xr
  y_pred <- predict(mod, newdata = grid)
  data.frame(x = xr, y = y_pred)
}

pd_list <- setNames(lapply(preds, make_pd, mod = fit_final, df = model_complete,
                           prds = preds, mv = mean_vals), preds)

pd_all <- bind_rows(Map(function(d, v) cbind(d, Variable = v), pd_list, preds))
pd_all$Variable_Label <- label_map[pd_all$Variable]
pd_all$Variable_Label[is.na(pd_all$Variable_Label)] <- pd_all$Variable[is.na(pd_all$Variable_Label)]

p30 <- ggplot(pd_all, aes(x = x, y = y)) +
  geom_line(colour = "#2196F3", linewidth = 0.8) +
  facet_wrap(~ Variable_Label, scales = "free_x", nrow = 2) +
  labs(
       x = "Predictor value", y = "Predicted cases per 1,000") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 9))

ggsave(file.path(fig_dir, "fig30_partial_dependence.png"), p30, width = 11, height = 7, dpi = 600)

cat("[OK] Figure 30 saved.\n"); rm(p30, pd_all, pd_list)

# ========================== FIGURE 31: Temporal Stability ==========================
cat("\n--- Figure 31: Temporal Stability ---\n")

weekly_r2 <- model_complete %>%
  group_by(week_start) %>%
  filter(n() > 5) %>%
  summarise(
    R2    = cor(fitted, weekly_cases_p1k, use = "complete.obs")^2,
    n_obs = n(),
    .groups = "drop"
  )

weekly_r2$wave <- sapply(weekly_r2$week_start, classify_wave)

wave_rects <- waves %>% 
  mutate(ymin = -Inf, ymax = Inf) %>%
  rename(wave_label = wave)

p31 <- ggplot(weekly_r2, aes(x = week_start, y = R2)) +
  geom_rect(data = wave_rects,
            aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = wave_label),
            alpha = 0.15, inherit.aes = FALSE) +
  geom_point(size = 1.5, colour = "#333333", na.rm = TRUE) +
  geom_smooth(method = "loess", span = 0.3, se = TRUE,
              colour = "#D32F2F", linewidth = 0.8, na.rm = TRUE) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 0,   linetype = "solid",   colour = "grey30") +
  scale_fill_manual(values = wcols_n,
                    guide = guide_legend(title = "Wave", override.aes = list(alpha = 0.3))) +
  scale_y_continuous(limits = c(-0.2, 1.05)) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "2 months") +
  labs(
       x = "Date", y = "Weekly R²", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "fig31_temporal_stability.png"), p31, width = 10, height = 5, dpi = 600)

cat("[OK] Figure 31 saved.\n"); rm(p31, wave_rects)

# ========================== FIGURE 32: Residual Spatial Autocorrelation ==========================
cat("\n--- Figure 32: Residual Spatial Autocorrelation ---\n")

sp_resid <- model_complete %>%
  group_by(SP_CODE) %>%
  summarise(mean_resid = mean(residuals, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(mean_resid))

cat("  Subplaces with residuals:", nrow(sp_resid), "\n")

resid_map <- sf_subplaces %>%
  left_join(sp_resid, by = "SP_CODE") %>%
  filter(!is.na(mean_resid))

sf_proj <- st_transform(resid_map, crs = 32736)

# FIX: Compute centroids for the knn fallback
sp_centroids <- st_centroid(sf_proj)

nb <- tryCatch(
  poly2nb(sf_proj, queen = TRUE, zero.policy = TRUE),
  error = function(e) {
    cat("  [WARN] poly2nb failed, using knn fallback on centroids\n")
    knn2nb(knearneigh(sp_centroids, k = 4))
  }
)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

moran_res <- tryCatch(
  moran.test(resid_map$mean_resid, lw, zero.policy = TRUE, na.action = na.exclude),
  error = function(e) {
    cat("[WARN] Moran test failed:", e$message, "\n")
    NULL
  }
)

if (!is.null(moran_res)) {
  mi <- moran_res$estimate["Moran I statistic"]
  mp <- moran_res$p.value
  cat(sprintf("\nMoran's I = %.4f (p = %.4f)\n", mi, mp))
  if (mp < 0.05) {
    cat("  -> Significant residual spatial autocorrelation detected.\n")
    cat("     RK captures most but not all spatial structure.\n")
  } else {
    cat("  -> No significant residual spatial autocorrelation.\n")
    cat("     RK adequately captures spatial dependence.\n")
  }
  mi_val <- round(mi, 3)
  mp_val <- round(mp, 3)
} else {
  mi_val <- NA
  mp_val <- NA
}

brk <- c(-max(abs(sp_resid$mean_resid), na.rm = TRUE), 
         max(abs(sp_resid$mean_resid), na.rm = TRUE))

p32 <- ggplot(resid_map) +
  geom_sf(aes(fill = mean_resid), colour = "grey50", linewidth = 0.05) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, limits = brk,
                       name = "Mean residual") +
  theme_void(base_size = 12) +
  theme(legend.position = "right", plot.title = element_text(face = "bold"))

ggsave(file.path(fig_dir, "fig32_residual_spatial.png"), p32, width = 9, height = 7, dpi = 600)
ggsave(file.path(fig_dir, "fig32_residual_spatial.pdf"), p32, width = 9, height = 7)
cat("[OK] Figure 32 saved.\n"); rm(p32, resid_map, sf_proj, sp_centroids)

# ========================== Tables ==========================

cat("\n% Table: Per-Wave Model Performance\n")
cat("\\begin{table}[htbp]\n\\centering\n")
cat("\\caption{Mean spatial CV performance by COVID-19 wave}\n\\label{tab:wave_perf}\n")
cat("\\begin{tabular}{llrrr}\n\\toprule\nWave & Model & R$^2$ & RMSE & N folds \\\\\n\\midrule\n")
for (i in seq_len(nrow(wave_all))) {
  cat(sprintf("  %s & %s & %.3f & %.4f & %d \\\\\n",
              wave_all$wave[i], wave_all$Model[i],
              wave_all$mean_R2[i], wave_all$mean_RMSE[i], wave_all$n_folds[i]))
}
cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n")

cat("\n% Table: Variable Importance\n")
cat("\\begin{table}[htbp]\n\\centering\n")
cat("\\caption{Predictor Relative Importance (Standardised Coefficients)}\n\\label{tab:var_imp}\n")
cat("\\begin{tabular}{lcr}\n\\toprule\nPredictor & Importance (\\%) & Coefficient \\\\\n\\midrule\n")
for (i in seq_len(nrow(imp_df))) {
  lab <- ifelse(is.na(imp_df$Variable_Label[i]), imp_df$Variable[i], imp_df$Variable_Label[i])
  cat(sprintf("  %s & %.1f & %.4f \\\\\n", lab, imp_df$Importance[i], imp_df$Coefficient[i]))
}
cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n")

cat("\n% Table: Moran's I Test\n")
cat("\\begin{table}[htbp]\n\\centering\n")
cat("\\caption{Moran's I test on model residuals}\n\\label{tab:moran}\n")
cat("\\begin{tabular}{lr}\n\\toprule\nStatistic & Value \\\\\n\\midrule\n")
if (!is.null(moran_res)) {
  cat(sprintf("  Moran's I & %.4f \\\\\n", mi))
  cat(sprintf("  Expected I & %.4f \\\\\n", moran_res$estimate["Expectation"]))
  cat(sprintf("  p-value & %.4f \\\\\n", mp))
  cat(sprintf("  Interpretation & %s \\\\\n",
              ifelse(mp < 0.05, "Significant residual autocorrelation",
                     "No significant autocorrelation")))
} else {
  cat("  Moran's I & N/A \\\\\n")
}
cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n")

# ========================== Save ==========================
save(imp_df, wave_all, wave_rk, weekly_r2, moran_res, fit_final,
     model_complete,
     file = file.path(out_dir, "extended_diagnostics_results.RData"))

cat("\n========== SCRIPT 9 COMPLETE ==========\n")
cat("Figures 28-32 saved to:", fig_dir, "\n")
cat("Results saved to:", out_dir, "\n")
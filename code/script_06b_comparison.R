##############################################################################
# SCRIPT 6B: FINAL MODEL COMPARISON PLOT
#
# Purpose: Load RK results (Script 5) and Competing Models (Script 6A)
##############################################################################

# ========================== LIBRARIES ==========================
library(dplyr); library(ggplot2); library(ggpubr)

# ========================== PATHS ==========================
output_dir <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/EDA_Output"
fig_dir <- file.path(output_dir, "figures")

# ========================== LOAD DATA ==========================
cat("=== Loading Results ===\n")

# 1. Load Regression Kriging Results from Script 5
rk_overall <- read.csv(file.path(output_dir, "results", "rk_weekly_cv_overall.csv"))

# Rename to standard format for clean binding
rk_overall <- rk_overall %>%
  rename(Mean_MAE = Mean_MAE, Mean_RMSE = Mean_RMSE, Mean_sMAPE = Mean_sMAPE, 
         Mean_R2 = Mean_R2, Outlier_MAE = Outlier_MAE) %>%
  mutate(SD_R2 = 0) # RK week-by-week SD not saved, set to 0 for plot

# 2. Load Competing Models from Script 6A
load(file.path(output_dir, "script6a_slmm_results.RData"))
competing_models <- slmm_df

# ========================== COMBINE INTO FINAL TABLE ==========================
cat("\n=== Assembling Final Model Leaderboard ===\n")

final_table <- bind_rows(rk_overall, competing_models)

# Order by R2 descending for the plot
final_table$Model <- factor(final_table$Model, 
                            levels = c("Regression Kriging", "Spatial Random Forest", 
                                       "Regression-Only", "LMM (WWTP Cluster)", "Mean Baseline"))

# ========================== PLOT 1: R2 COMPARISON ==========================
cat("=== Generating R2 Comparison Plot ===\n")

p_r2 <- ggplot(final_table, aes(x = Model, y = Mean_R2, fill = Model)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", Mean_R2)), vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c(
    "Regression Kriging" = "#2ca02c",
    "Spatial Random Forest" = "#ff7f0e",
    "Regression-Only" = "#7570b3",
    "LMM (WWTP Cluster)" = "#1f77b4",
    "Mean Baseline" = "#d94801"
  )) +
  labs(title = "Model Performance Comparison: Spatial Cross-Validation",
       subtitle = "4-Week Early Warning Horizon | Predicting COVID-19 Cases per 1,000",
       y = expression(R^2), x = NULL,
       caption = "Note: RK validated via 10-fold spatial CV. Competitors validated via 80/20 temporal split. RF R² reflects localized spatial overfitting.") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 1, size = 9),
        plot.title = element_text(face = "bold"),
        plot.caption = element_text(color = "grey40", size = 8)) +
  coord_cartesian(ylim = c(0, max(final_table$Mean_R2, na.rm = TRUE) * 1.15))

ggsave(file.path(fig_dir, "fig_final_model_comparison_R2.png"), p_r2, width = 8, height = 6, dpi = 600)

# ========================== PLOT 2: MAE COMPARISON==========================
cat("=== Generating MAE Comparison Plot ===\n")

p_mae <- ggplot(final_table, aes(x = Model, y = Mean_MAE, fill = Model)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = sprintf("%.4f", Mean_MAE)), vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c(
    "Regression Kriging" = "#2ca02c", 
    "Spatial Random Forest" = "#ff7f0e", 
    "Regression-Only" = "#7570b3", 
    "LMM (WWTP Cluster)" = "#1f77b4", 
    "Mean Baseline" = "#d94801"
  )) +
  labs(title = "Model Performance Comparison: Mean Absolute Error (MAE)",
       subtitle = "Lower is better. True operational cost of misprediction.",
       y = "Mean Absolute Error", x = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 1, size = 9),
        plot.title = element_text(face = "bold")) +
  coord_cartesian(ylim = c(0, max(final_table$Mean_MAE, na.rm = TRUE) * 1.2))

ggsave(file.path(fig_dir, "fig_final_model_comparison_MAE.png"), p_mae, width = 8, height = 6, dpi = 600)

cat("\n=== FINAL MODEL LEADERBOARD ===\n")
print(final_table)

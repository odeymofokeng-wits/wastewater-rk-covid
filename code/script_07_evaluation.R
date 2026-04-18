# SCRIPT 7: FINAL EVALUATION, DIAGNOSTICS
##############################################################################
            
# ========================== LIBRARIES ==========================
            library(sf); library(dplyr); library(tidyr); library(ggplot2)
            library(gstat); library(sp); library(spdep); library(MASS)
            library(gridExtra); library(viridis); library(scales)
            library(lubridate); library(knitr); library(kableExtra)
            library(cowplot)
            
# ========================== PATHS ==========================
            path_subplaces <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/WASTEWATER DATA/Ekurhuleni sub and main places/Eku_sub_places13.shp"
            output_dir <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/EDA_Output"
            fig_dir <- file.path(output_dir, "figures")
            results_dir <- file.path(output_dir, "results")
            
            # ========================== LOAD DATA ==========================
            cat("=== Loading Data ===\n")
            
            # Load the base spatial data
            sf_subplaces <- st_read(path_subplaces, quiet = TRUE)
            
            rk_csv <- file.path(results_dir, "rk_weekly_cv_overall.csv")
            cat("Loading RK results from:", basename(rk_csv), "\n")
            rk_df <- read.csv(rk_csv)
            cat("  Loaded", nrow(rk_df), "rows\n")
            
            comp_csv <- file.path(output_dir, "competing_models_results.csv")
            cat("Loading Competing Models from:", basename(comp_csv), "\n")
            competing_df <- read.csv(comp_csv)
            cat("  Loaded", nrow(competing_df), "rows\n")
            
            # Add Mean Baseline if missing
            if (!"Mean Baseline" %in% competing_df$Model) {
              mean_baseline <- data.frame(
                Model = "Mean Baseline",
                MAE = 0.8688773,  # From your slmm_df
                RMSE = 1.370184,
                sMAPE = 101.67972,
                R2 = -0.122899886,
                Outlier_MAE = 3.203014
              )
              competing_df <- rbind(competing_df, mean_baseline)
            }
            
            # Load variable selection
            var_table_path <- file.path(output_dir, "variable_selection_table.csv")
            var_table <- read.csv(var_table_path)
            cat("  Loaded variable selection table:", nrow(var_table), "rows\n")
            
            # Load master analysis for model fitting
            load(file.path(output_dir, "master_analysis.RData"))
            cat("  Loaded master_analysis.RData\n")
            
            # ========================== PREPARE MODEL DATA ==========================
            cat("\n=== Preparing Model Data ===\n")
            
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
            
            model_data <- weekly_covid %>%
              left_join(weekly_gene, by = c("WWTP_Name", "week_start")) %>%
              filter(weekly_cases_p1k > 0) %>%
              mutate(log_gene = log1p(mean_gene)) %>%
              group_by(SP_CODE) %>% arrange(week_start) %>%
              mutate(
                gene_locf = zoo::na.locf(log_gene, na.rm = FALSE),
                gene_locf = zoo::na.locf(gene_locf, fromLast = TRUE, na.rm = FALSE),
                gene_lag1 = dplyr::lag(gene_locf, 1),
                cases_lag1 = dplyr::lag(weekly_cases_p1k, 1)
              ) %>%
              ungroup() %>%
              filter(!is.na(gene_lag1), !is.na(cases_lag1))
            
            cat("  Prepared model data:", nrow(model_data), "rows\n")
            
            # Fit final model for diagnostics
            fit_final <- lm(weekly_cases_p1k ~ gene_lag1 + cases_lag1 + covid_vuln + health_vuln + hosp_dist + log_pop,
                            data = model_data)
            cat("  Model fitted with", nrow(model_data), "observations\n")
            cat("  Model R-squared:", round(summary(fit_final)$r.squared, 4), "\n")
            
            # ========================== FIGURE 1: RK Model Diagnostics ==========================
            cat("\n=== Figure 1: RK Residual Diagnostics ===\n")
            
            final_fit_residuals <- residuals(fit_final)
            final_data_pred <- fitted(fit_final)
            
            diag_data <- data.frame(
              residuals = final_fit_residuals,
              fitted = final_data_pred,
              observed = model_data$weekly_cases_p1k,
              gene_lag1 = model_data$gene_lag1,
              hosp_dist = model_data$hosp_dist
            )
            
            p1 <- ggplot(diag_data, aes(sample = residuals)) +
              stat_qq(alpha = 0.5, color = "steelblue") +
              stat_qq_line(color = "red", linewidth = 0.8) +
              labs(
                   x = "Theoretical Quantiles", y = "Sample Quantiles") +
              theme_minimal(base_size = 10)
            
            p2 <- ggplot(diag_data, aes(x = fitted, y = residuals)) +
              geom_point(alpha = 0.3, size = 0.8, color = "steelblue") +
              geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
              geom_smooth(method = "loess", color = "darkgreen", se = FALSE, linewidth = 0.8) +
              labs(x = "Fitted Values", y = "Residuals") +
              theme_minimal(base_size = 10)
            
            p3 <- ggplot(diag_data, aes(x = residuals)) +
              geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "steelblue", 
                             color = "white", alpha = 0.7) +
              geom_density(color = "red", linewidth = 0.8) +
              labs( x = "Residuals", y = "Density") +
              theme_minimal(base_size = 10)
            
            p4 <- ggplot(diag_data, aes(x = fitted, y = observed)) +
              geom_point(alpha = 0.3, size = 0.8, color = "steelblue") +
              geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
              geom_smooth(method = "lm", color = "darkgreen", se = TRUE, alpha = 0.2) +
              labs( x = "Predicted Cases/1000", y = "Observed Cases/1000") +
              coord_equal() +
              theme_minimal(base_size = 10)
            
            p5 <- ggplot(diag_data, aes(x = gene_lag1, y = residuals)) +
              geom_point(alpha = 0.3, size = 0.8, color = "steelblue") +
              geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
              geom_smooth(method = "loess", color = "darkgreen", se = FALSE, linewidth = 0.8) +
              labs(x = "Log Gene Copies/mL", y = "Residuals") +
              theme_minimal(base_size = 10)
            
            p6 <- ggplot(diag_data, aes(x = hosp_dist, y = residuals)) +
              geom_point(alpha = 0.3, size = 0.8, color = "steelblue") +
              geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
              geom_smooth(method = "loess", color = "darkgreen", se = FALSE, linewidth = 0.8) +
              labs( x = "Distance (km)", y = "Residuals") +
              theme_minimal(base_size = 10)
            
            diag_grid <- plot_grid(p1, p2, p3, p4, p5, p6, ncol = 3, 
                                   labels = paste0("(", letters[1:6], ")"))
            
            ggsave(file.path(fig_dir, "fig_rk_diagnostics.png"), diag_grid, 
                   width = 14, height = 10, dpi = 600)
           
            
            cat("  Saved: fig_rk_diagnostics\n")
            
            # ========================== FIGURE 2: Model Performance Comparison ==========================
            cat("\n=== Figure 2: Model Performance Comparison ===\n")
            
        
            rk_summary <- rk_df %>%
              select(Model, Mean_MAE, Mean_RMSE, Mean_sMAPE, Mean_R2) %>%
              mutate(Analysis = "Spatial CV (RK)")
            
        
            comp_summary <- competing_df %>%
              select(Model, MAE, RMSE, sMAPE, R2) %>%
              rename(Mean_MAE = MAE, Mean_RMSE = RMSE, Mean_sMAPE = sMAPE, Mean_R2 = R2) %>%
              mutate(Analysis = "Temporal Test (RF/LMM/OLS)")
            
            comparison_data <- bind_rows(rk_summary, comp_summary)
            
            # R2 Comparison
            p_r2 <- ggplot(comparison_data, aes(x = Model, y = Mean_R2, fill = Model)) +
              geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
              geom_text(aes(label = sprintf("%.3f", Mean_R2)), vjust = -0.5, size = 3.5, fontface = "bold") +
              facet_wrap(~Analysis, scales = "free_x") +
              scale_fill_brewer(palette = "Set1") +
              labs(
                   y = expression(R^2), x = NULL) +
              theme_minimal(base_size = 11) +
              theme(axis.text.x = element_text(angle = 30, hjust = 1),
                    legend.position = "none")
            
            ggsave(file.path(fig_dir, "fig_model_comparison_R2.png"), p_r2, width = 10, height = 6, dpi = 600)
            
            cat("  Saved: fig_model_comparison_R2\n")
            
            # MAE Comparison
            p_mae <- ggplot(comparison_data, aes(x = Model, y = Mean_MAE, fill = Model)) +
              geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
              geom_text(aes(label = sprintf("%.4f", Mean_MAE)), vjust = -0.5, size = 3.5, fontface = "bold") +
              facet_wrap(~Analysis, scales = "free_x") +
              scale_fill_brewer(palette = "Set1") +
              labs(
                   y = "MAE", x = NULL) +
              theme_minimal(base_size = 11) +
              theme(axis.text.x = element_text(angle = 30, hjust = 1),
                    legend.position = "none")
            
            ggsave(file.path(fig_dir, "fig_model_comparison_MAE.png"), p_mae, width = 10, height = 6, dpi = 600)
            cat("  Saved: fig_model_comparison_MAE\n")
            
            # ========================== FIGURE 3: Spatial Residual Variogram ==========================
            cat("\n=== Figure 3: Spatial Residual Variogram ===\n")
            
            sf_coords <- sf_subplaces %>%
              st_transform(crs = 32736) %>%
              st_centroid() %>%
              st_coordinates()
            
            sp_residuals <- data.frame(
              SP_CODE = model_data$SP_CODE,
              residuals = final_fit_residuals
            ) %>%
              group_by(SP_CODE) %>%
              summarise(resid_mean = mean(residuals, na.rm = TRUE), .groups = "drop")
            
            coord_idx <- match(sp_residuals$SP_CODE, sf_subplaces$SP_CODE)
            valid_coords <- !is.na(coord_idx)
            
            cat("  Valid coordinates:", sum(valid_coords), "of", length(coord_idx), "\n")
            
            if (sum(valid_coords) > 10) {
              # Create SpatialPointsDataFrame
              spdf_resid <- SpatialPointsDataFrame(
                coords = sf_coords[coord_idx[valid_coords], ],
                data = sp_residuals[valid_coords, ]
              )
              
              # Compute variogram
              vg_resid <- tryCatch(variogram(resid_mean ~ 1, data = spdf_resid), error = function(e) NULL)
              
              if (!is.null(vg_resid) && nrow(vg_resid) >= 3) {
                cat("  Variogram computed with", nrow(vg_resid), "bins\n")
                
                # Convert to data frame for ggplot
                vg_df <- as.data.frame(vg_resid)
                
                # Start with empirical variogram plot
                p_vg <- ggplot(vg_df, aes(x = dist, y = gamma)) +
                  geom_point(color = "steelblue", size = 3, alpha = 0.7) +
                  geom_segment(aes(x = dist - dist/10, xend = dist + dist/10, y = gamma, yend = gamma), 
                               color = "steelblue", alpha = 0.5) +
                  labs(
                       x = "Distance (m)", 
                       y = "Semivariance") +
                  theme_minimal(base_size = 12) +
                  theme(plot.title = element_text(face = "bold"))
                
                # Try to fit variogram model
                vg_fit_resid <- tryCatch(
                  fit.variogram(vg_resid, vgm(psill = var(sp_residuals$resid_mean[valid_coords], na.rm = TRUE) * 0.5,
                                              "Exp", range = max(vg_resid$dist) / 3, nugget = 0.1)),
                  error = function(e) {
                    cat("  Variogram fit warning:", conditionMessage(e), "\n")
                    NULL
                  }
                )
                
                # Add fitted line if successful
                if (!is.null(vg_fit_resid) && !is.na(vg_fit_resid$psill[2]) && vg_fit_resid$psill[2] > 0) {
                  # Generate smooth fitted line using gstat's variogramLine
                  fit_df <- variogramLine(vg_fit_resid, maxdist = max(vg_df$dist, na.rm = TRUE), n = 200)
                  
                  p_vg <- p_vg + 
                    geom_line(data = fit_df, aes(x = dist, y = gamma), 
                              color = "red", linewidth = 1.2, linetype = "dashed") +
                    annotate("text", 
                             x = max(vg_df$dist, na.rm = TRUE) * 0.65, 
                             y = max(vg_df$gamma, na.rm = TRUE) * 0.9,
                             label = sprintf("Nugget: %.3f\nPartial Sill: %.3f\nRange: %.0f m", 
                                             vg_fit_resid$psill[1],
                                             vg_fit_resid$psill[2],
                                             vg_fit_resid$range[2]),
                             hjust = 0, size = 3, color = "darkred")
                  
                  cat("  Fitted exponential model added\n")
                } else {
                  p_vg <- p_vg + annotate("text", 
                                          x = max(vg_df$dist, na.rm = TRUE) * 0.5,
                                          y = max(vg_df$gamma, na.rm = TRUE) * 0.9,
                                          label = "Model fit failed",
                                          color = "red", size = 4)
                }
                
                # Save 
                save_path <- file.path(fig_dir, "fig_residual_variogram.png")
                cat("  Saving to:", save_path, "\n")
                
                
              } else {
                cat("  Could not compute variogram (insufficient data)\n")
              }
            } else {
              cat("  Insufficient spatial matches for variogram:", sum(valid_coords), "\n")
            }
            
            # ========================== FIGURE 4: Variable Importance ==========================
            cat("\n=== Figure 4: Variable Importance ===\n")
            
            coef_summary <- summary(fit_final)$coefficients
            coef_df <- data.frame(
              Variable = rownames(coef_summary)[-1],
              Coefficient = coef_summary[-1, 1],
              StdError = coef_summary[-1, 2],
              PValue = coef_summary[-1, 4]
            ) %>%
              mutate(
                Significant = PValue < 0.05,
                AbsCoef = abs(Coefficient)
              )
            
            p_coef <- ggplot(coef_df, aes(x = reorder(Variable, AbsCoef), y = Coefficient, fill = Significant)) +
              geom_bar(stat = "identity", alpha = 0.8) +
              geom_errorbar(aes(ymin = Coefficient - 1.96*StdError, 
                                ymax = Coefficient + 1.96*StdError), width = 0.2) +
              coord_flip() +
              scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "grey70"),
                                name = "p < 0.05") +
              labs(
                   x = NULL, y = "Coefficient Value") +
              theme_minimal(base_size = 11)
            
            ggsave(file.path(fig_dir, "fig_variable_importance.png"), p_coef, width = 8, height = 6, dpi = 600)
            ggsave(file.path(fig_dir, "fig_variable_importance.pdf"), p_coef, width = 8, height = 6)
            cat("  Saved: fig_variable_importance\n")
            
            # ========================== FIGURE 5: Temporal Trends ==========================
            cat("\n=== Figure 5: Temporal Trends ===\n")
            
            temporal_summary <- analysis_data %>%
              group_by(report_date) %>%
              summarise(
                total_cases = sum(daily_cases_sp, na.rm = TRUE),
                mean_gene = mean(gene_copies_ml, na.rm = TRUE),
                n_obs = n(),
                .groups = "drop"
              ) %>%
              mutate(
                cases_ma7 = zoo::rollmeanr(total_cases, k = 7, fill = NA),
                gene_ma7 = zoo::rollmeanr(mean_gene, k = 7, fill = NA)
              )
            
            p_temp <- ggplot(temporal_summary) +
              geom_bar(aes(x = report_date, y = total_cases / 1000),
                       stat = "identity", fill = "steelblue", alpha = 0.3, width = 1) +
              geom_line(aes(x = report_date, y = cases_ma7 / 1000), color = "blue", linewidth = 0.8) +
              geom_line(aes(x = report_date, y = gene_ma7 * 100), color = "red", linewidth = 0.8) +
              scale_x_date(date_labels = "%b %Y", date_breaks = "2 months") +
              scale_y_continuous(
                name = "COVID-19 Cases (thousands)",
                sec.axis = sec_axis(~ . / 100, name = "Gene Copies/mL (scaled)")
              ) +
              
              theme_minimal(base_size = 11) +
              theme(axis.text.x = element_text(angle = 45, hjust = 1))
            
            ggsave(file.path(fig_dir, "fig_temporal_trends.png"), p_temp, width = 14, height = 6, dpi = 600)
            ggsave(file.path(fig_dir, "fig_temporal_trends.pdf"), p_temp, width = 14, height = 6)
            cat("  Saved: fig_temporal_trends\n")
            
            # ========================== SUMMARY TABLES ==========================
            cat("\n=== Generating Summary Tables ===\n")
            
            cat("\n--- Table 1: RK Performance Summary ---\n")
            rk_table <- rk_df %>%
              select(Model, Mean_MAE, Mean_RMSE, Mean_sMAPE, Mean_R2, Outlier_MAE)
            print(rk_table)
            
            cat("\n--- Table 2: Competing Models Summary ---\n")
            comp_table <- competing_df %>%
              select(Model, MAE, RMSE, sMAPE, R2, Outlier_MAE)
            print(comp_table)
            
            cat("\n--- Table 3: Final Model Comparison ---\n")
            print(comparison_data %>% select(Analysis, Model, Mean_R2, Mean_MAE, Mean_RMSE))
            
            cat("\n--- Table 4: Variable Selection Results ---\n")
            print(var_table)
            
            # ========================== SAVE FINAL RESULTS ==========================
            cat("\n=== Saving Final Results ===\n")
            
            save(rk_df, competing_df, comparison_data, var_table, fit_final,
                 file = file.path(output_dir, "script7_final_results.RData"))
            cat("  Saved: script7_final_results.RData\n")
            
            write.csv(comparison_data, file.path(results_dir, "final_model_comparison.csv"), row.names = FALSE)
            cat("  Saved: final_model_comparison.csv\n")
            
          
            cat("  SCRIPT 7 COMPLETE\n")
           
            

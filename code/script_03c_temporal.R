##############################################################################
# SCRIPT 3C: EDA TEMPORAL & DISTRIBUTIONS (Part 3 of 3) — Figures 13 to 18
#
# Figures:
#   13. COVID wave analysis
#   14. Sub-district comparison
#   15. Population distribution
#   16. Cases per 1000 distribution
#   17. Cross-correlation (gene vs cases)
#   18. Vulnerability vs Cases scatter
#
# REQUIRES: master_analysis.RData + script3b_temporal.RData
##############################################################################

# ========================== LIBRARIES ==========================
library(sf); library(dplyr); library(ggplot2)
library(viridis); library(scales); library(lubridate); library(zoo)

# ========================== PATHS ==========================
output_dir <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/EDA_Output"
fig_dir <- file.path(output_dir, "figures")

# ========================== LOAD ==========================
cat("=== Loading Data ===\n")
load(file.path(output_dir, "master_analysis.RData"))
load(file.path(output_dir, "script3b_temporal.RData"))

# ========================== FIGURE 13: Wave Analysis ==========================
cat("=== Figure 13: COVID Wave Analysis ===\n")

max_cases <- max(temporal_data$total_daily_cases, na.rm = TRUE) / 1000

p13 <- ggplot(temporal_data %>% filter(!is.na(cases_ma7))) +
  # Wave shading
  annotate("rect", xmin = as.Date("2020-06-01"), xmax = as.Date("2020-10-01"),
           ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.1) +
  annotate("rect", xmin = as.Date("2020-11-01"), xmax = as.Date("2021-03-01"),
           ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.1) +
  annotate("rect", xmin = as.Date("2021-05-01"), xmax = as.Date("2021-08-01"),
           ymin = -Inf, ymax = Inf, fill = "purple", alpha = 0.1) +
  annotate("rect", xmin = as.Date("2021-11-01"), xmax = as.Date("2022-01-31"),
           ymin = -Inf, ymax = Inf, fill = "darkred", alpha = 0.1) +
  # Data
  geom_bar(aes(x = report_date, y = total_daily_cases / 1000),
           stat = "identity", fill = "steelblue", alpha = 0.4, width = 1) +
  geom_line(aes(x = report_date, y = cases_ma7 / 1000), color = "blue", linewidth = 0.8) +
  geom_line(aes(x = report_date, y = gene_ma7 * 100), color = "red", linewidth = 0.8) +
  # Wave labels
  annotate("text", x = as.Date("2020-08-01"), y = max_cases * 0.9,
           label = "Wave 1", size = 3, fontface = "bold") +
  annotate("text", x = as.Date("2021-01-15"), y = max_cases * 0.9,
           label = "Wave 2", size = 3, fontface = "bold") +
  annotate("text", x = as.Date("2021-06-15"), y = max_cases * 0.9,
           label = "Wave 3\n(Delta)", size = 3, fontface = "bold") +
  annotate("text", x = as.Date("2021-12-15"), y = max_cases * 0.9,
           label = "Wave 4\n(Omicron)", size = 3, fontface = "bold") +
  scale_x_date(date_labels = "%b %Y", date_breaks = "2 months") +
  labs(
       x = NULL, y = "Cases (thousands) / Gene signal (scaled)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "fig13_wave_analysis.png"), p13, width = 14, height = 6, dpi = 600)
cat("  Saved: fig13_wave_analysis\n")

# ========================== FIGURE 14: Sub-District Comparison ==========================
cat("=== Figure 14: Sub-District Comparison ===\n")

sd_temporal <- analysis_data %>%
  group_by(health_subdistrict, report_date) %>%
  summarise(daily_cases = sum(daily_cases_sp, na.rm = TRUE), .groups = "drop")

p14 <- ggplot(sd_temporal, aes(x = report_date, y = daily_cases, color = health_subdistrict)) +
  geom_line(linewidth = 0.5, alpha = 0.7) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
  facet_wrap(~health_subdistrict, scales = "free_y", ncol = 2) +
  labs(
       x = NULL, y = "Daily Cases") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "none", strip.text = element_text(face = "bold"))

ggsave(file.path(fig_dir, "fig14_subdistrict_comparison.png"), p14, width = 12, height = 10, dpi = 600)

cat("  Saved: fig14_subdistrict_comparison\n")

# ========================== FIGURE 15: Population Distribution ==========================
cat("=== Figure 15: Population Distribution ===\n")

med_pop <- median(sp_totals$pop, na.rm = TRUE)

p15 <- ggplot(sp_totals, aes(x = pop / 1000)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_vline(xintercept = med_pop / 1000, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
       x = "Estimated Population (thousands)", y = "Number of Subplaces") +
  theme_minimal(base_size = 11)

ggsave(file.path(fig_dir, "fig15_population_distribution.png"), p15, width = 8, height = 5, dpi = 600)
cat("  Saved: fig15_population_distribution\n")

# ========================== FIGURE 16: Cases per 1000 Distribution ==========================
cat("=== Figure 16: Cases per 1000 Distribution ===\n")

p16 <- ggplot(sp_totals %>% filter(mean_cases_per_1000 > 0),
              aes(x = mean_cases_per_1000)) +
  geom_histogram(bins = 40, fill = "darkorange", color = "white", alpha = 0.8) +
  scale_x_continuous(trans = "log10") +
  geom_vline(xintercept = median(sp_totals$mean_cases_per_1000, na.rm = TRUE),
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(
       x = "Mean Daily Cases per 1,000 (log)", y = "Number of Subplaces") +
  theme_minimal(base_size = 11)

ggsave(file.path(fig_dir, "fig16_cases_per_1000_distribution.png"), p16, width = 8, height = 5, dpi = 600)
cat("  Saved: fig16_cases_per_1000_distribution\n")

# ========================== FIGURE 17: Cross-Correlation ==========================
cat("=== Figure 17: Cross-Correlation ===\n")

weekly_data <- analysis_data %>%
  group_by(weekly_number) %>%
  summarise(
    total_weekly_cases = sum(daily_cases_sp, na.rm = TRUE),
    mean_weekly_gene = mean(gene_copies_ml, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(mean_weekly_gene), mean_weekly_gene > 0, total_weekly_cases > 0)

cat("  Weekly data with gene:", nrow(weekly_data), "\n")

if (nrow(weekly_data) > 20) {
  ccf_result <- ccf(weekly_data$mean_weekly_gene, weekly_data$total_weekly_cases,
                    lag.max = 8, plot = FALSE)
  ccf_df <- data.frame(lag = ccf_result$lag, correlation = ccf_result$acf)

  best_lag <- ccf_df$lag[which.max(ccf_df$correlation)]
  best_cor <- max(ccf_df$correlation)
  cat("  Best lag:", best_lag, "weeks, r =", round(best_cor, 3), "\n")

  p17 <- ggplot(ccf_df, aes(x = lag, y = correlation)) +
    geom_col(aes(fill = correlation > 0), width = 0.6, alpha = 0.7) +
    geom_hline(yintercept = c(-0.2, 0.2), color = "grey50", linetype = "dashed") +
    geom_vline(xintercept = best_lag, color = "red", linetype = "dashed", linewidth = 0.8) +
    annotate("text", x = best_lag, y = max(ccf_df$correlation) * 0.95,
             label = sprintf("Peak: lag=%d, r=%.3f", best_lag, best_cor),
             size = 3, hjust = 0) +
    scale_fill_manual(values = c("steelblue", "darkorange"), guide = "none") +
    labs(
         x = "Lag (weeks)", y = "Cross-Correlation") +
    theme_minimal(base_size = 11)

  ggsave(file.path(fig_dir, "fig17_cross_correlation.png"), p17, width = 8, height = 5, dpi = 600)
  ggsave(file.path(fig_dir, "fig17_cross_correlation.pdf"), p17, width = 8, height = 5)
  cat("  Saved: fig17_cross_correlation\n")
}

# ========================== FIGURE 18: Vulnerability vs Cases ==========================
cat("=== Figure 18: Vulnerability vs Cases ===\n")

p18 <- ggplot(sp_totals %>% filter(mean_cases_per_1000 > 0)) +
  geom_point(aes(x = vulnerability, y = mean_cases_per_1000,
                 color = covid_vuln, size = pop), alpha = 0.6) +
  geom_smooth(aes(x = vulnerability, y = mean_cases_per_1000),
              method = "lm", color = "red", se = TRUE, linewidth = 0.8) +
  scale_y_continuous(trans = "log10") +
  scale_size_continuous(name = "Population", range = c(1, 8), labels = comma) +
  scale_color_viridis_c(name = "COVID Vuln. Index", option = "plasma") +
  labs(
       x = "Composite Vulnerability Score",
       y = "Mean Cases per 1 000 Population") +
  theme_minimal(base_size = 11)

ggsave(file.path(fig_dir, "fig18_vulnerability_vs_cases.png"), p18, width = 8, height = 6, dpi = 600)
cat("  Saved: fig18_vulnerability_vs_cases\n")

# ========================== SUMMARY ==========================
cat("\n=== Summary Statistics ===\n")

cat("\n  Descriptive stats:\n")
cat("  Total subplaces:", nrow(sp_totals), "\n")
cat("  Date range:", min(analysis_data$report_date), "to", max(analysis_data$report_date), "\n")
cat("  Total obs:", format(nrow(analysis_data), big.mark = ","), "\n")
cat("  Mean daily cases/sp:", round(mean(analysis_data$daily_cases_sp, na.rm = TRUE), 4), "\n")
cat("  Mean gene/mL:", round(mean(analysis_data$gene_copies_ml, na.rm = TRUE), 2), "\n")

cat("\n  Moran's I:\n")
cat("  Total cases:  I =", round(moran_cases$estimate[1], 4), "\n")
cat("  Mean gene:    I =", round(moran_gene$estimate[1], 4), "\n")
cat("  Cases/1000:   I =", round(moran_cp1k$estimate[1], 4), "\n")

cat("\n=== SCRIPT 3C COMPLETE ===\n")
cat("All 18 EDA figures saved to:", fig_dir, "\n")


cat("\n=== Running Script 3C Extras ===\n")
suppressPackageStartupMessages(library(patchwork))

# ---- Extra Figure 1: Per-WWTP Cross-Correlation Heatmap ----
cat("Computing per-WWTP CCFs...\n")
weekly_wwtp <- analysis_data %>%
  group_by(WWTP_Name, weekly_number) %>%
  summarise(weekly_gene = mean(gene_copies_ml, na.rm = TRUE),
            weekly_cases = sum(daily_cases_sp, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(weekly_gene), !is.na(weekly_cases))

wwtps <- sort(unique(weekly_wwtp$WWTP_Name))
lag_max <- 6
ccf_mat <- matrix(NA, length(wwtps), 2 * lag_max + 1,
                  dimnames = list(wwtps, -lag_max:lag_max))

for (w in wwtps) {
  sub <- weekly_wwtp %>% filter(WWTP_Name == w) %>% arrange(weekly_number)
  if (nrow(sub) > lag_max + 2) {
    cr <- ccf(sub$weekly_gene, sub$weekly_cases, lag.max = lag_max, plot = FALSE)
    ccf_mat[w, as.character(cr$lag)] <- cr$acf
  }
}

ccf_df <- as.data.frame(as.table(ccf_mat), stringsAsFactors = FALSE)
colnames(ccf_df) <- c("WWTP", "lag", "correlation")
ccf_df$lag <- as.numeric(ccf_df$lag)
ccf_df$WWTP <- factor(ccf_df$WWTP, levels = wwtps)

best_lag <- ccf_df %>% group_by(WWTP) %>% slice_max(abs(correlation), n = 1)
ccf_df$is_best <- interaction(ccf_df$WWTP, ccf_df$lag) %in%
  interaction(best_lag$WWTP, best_lag$lag)

p_ccf <- ggplot(ccf_df, aes(x = lag, y = WWTP, fill = correlation)) +
  geom_tile(color = "grey70") +
  geom_tile(data = filter(ccf_df, is_best), fill = NA, color = "black", linewidth = 1) +
  geom_text(aes(label = ifelse(is_best, sprintf("%.2f*", correlation),
                               ifelse(abs(correlation) > 0.3, sprintf("%.2f", correlation), ""))),
            size = 2.3) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, limits = c(-1, 1), name = "Correlation") +
  labs(
       x = "Lag (weeks)", y = NULL) +
  theme_minimal(base_size = 10) + theme(axis.text.y = element_text(size = 7))

ggsave(file.path(fig_dir, "fig19_ccf_by_wwtp.png"), p_ccf,
       width = 11, height = max(5, length(wwtps) * 0.3), dpi = 600)
cat("  Saved fig19_ccf_by_wwtp\n")

# ---- Extra Figure 2: Lagged Regression R-squared vs Lead Time ----
cat("Computing lagged regressions...\n")
weekly_all <- analysis_data %>%
  group_by(weekly_number) %>%
  summarise(gene = mean(gene_copies_ml, na.rm = TRUE),
            cases = sum(daily_cases_sp, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(gene), !is.na(cases))

lag_res <- data.frame(lag = 0:lag_max, r_sq = NA_real_, slope = NA_real_)
for (l in 0:lag_max) {
  df <- if (l == 0) weekly_all else
    data.frame(gene = head(weekly_all$gene, nrow(weekly_all) - l),
               cases = tail(weekly_all$cases, nrow(weekly_all) - l))
  fit <- lm(cases ~ gene, data = df)
  lag_res$r_sq[l + 1]  <- summary(fit)$r.squared
  lag_res$slope[l + 1] <- coef(fit)[2]
}

opt <- which.max(lag_res$r_sq)
cat(sprintf("  Optimal lag: %d weeks (R² = %.4f, slope = %.4f)\n",
            lag_res$lag[opt], lag_res$r_sq[opt], lag_res$slope[opt]))

p_r2 <- ggplot(lag_res, aes(x = lag, y = r_sq)) +
  geom_line(linewidth = 1, color = "#2166ac") + 
  geom_point(size = 3) +
  geom_point(data = lag_res[opt, ], color = "red", size = 5, shape = 18) +
  geom_text(data = lag_res[opt, ], 
            label = paste0("R2=", round(lag_res$r_sq[opt], 3)),
            vjust = -1.2, color = "red", fontface = "bold") +
  labs(x = "Lag (weeks)", y = "R-squared") +
  ylim(0, max(lag_res$r_sq, na.rm = TRUE) * 1.15)

p_sl <- ggplot(lag_res, aes(x = lag, y = slope)) +
  geom_line(linewidth = 1, color = "#b2182b") + 
  geom_point(size = 3) +
  labs(x = "Lag (weeks)", y = "Slope")

p_lag <- p_r2 + p_sl + 
  plot_layout(ncol = 2) +
  plot_annotation(title = NULL)


ggsave(file.path(fig_dir, "fig20_lagged_r2.png"), p_lag, width = 12, height = 5, dpi = 600)

cat("  Saved fig20_lagged_r2\n")



# ========================== FIGURE 9, 19, 11, 12 COMBINED ==========================

# ---- Figure (a): Correlation Heatmap (ggplot version matching corrplot style) ----
cat("=== Figure (a): Correlation Heatmap ===\n")

corr_vars <- analysis_data %>%
  select(daily_cases_sp, gene_copies_ml, cases_per_1000,
         wwtp_dist_km, hosp_dist_km,
         COVID_Vulnerability_Index, Health_Service_Vulnerability_Index,
         Total_Population_Vulnerability_Index, Estimated_Population_2018) %>%
  na.omit()

cor_mat <- cor(corr_vars, use = "pairwise.complete.obs")
short_names <- c("Daily Cases", "Gene/mL", "Cases/1000", "WWTP Dist",
                 "Hosp Dist", "COVID Vuln.", "Health Vuln.", "Pop Vuln.", "Pop")
colnames(cor_mat) <- short_names
rownames(cor_mat) <- short_names

# Melt to long format for ggplot
library(reshape2)
cor_df <- melt(cor_mat)
cor_df$Var1 <- factor(cor_df$Var1, levels = rev(short_names))
cor_df$Var2 <- factor(cor_df$Var2, levels = short_names)

# Keep only upper triangle (including diagonal for full matrix look)
cor_df$keep <- as.numeric(cor_df$Var1) <= as.numeric(cor_df$Var2)
cor_df$value[cor_df$keep == FALSE] <- NA

p_a <- ggplot(cor_df, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(is.na(value), "", sprintf("%.2f", value))),
            size = 2.8, color = "black") +
  scale_fill_gradientn(colors = colorRampPalette(brewer.pal(11, "RdBu"))(200),
                       limits = c(-1, 1), name = "r") +
  coord_fixed() +
  labs(x = NULL, y = NULL, tag = "(a)") +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    plot.tag = element_text(size = 12, margin = margin(t = 10)),
    plot.tag.position = "bottom",
    panel.grid = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    plot.margin = margin(5, 5, 10, 5)
  )

# ---- Figure (b): CCF by WWTP ----
cat("=== Figure (b): CCF by WWTP ===\n")

weekly_wwtp <- analysis_data %>%
  group_by(WWTP_Name, weekly_number) %>%
  summarise(weekly_gene = mean(gene_copies_ml, na.rm = TRUE),
            weekly_cases = sum(daily_cases_sp, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(weekly_gene), !is.na(weekly_cases))

wwtps <- sort(unique(weekly_wwtp$WWTP_Name))
lag_max <- 6
ccf_mat <- matrix(NA, length(wwtps), 2 * lag_max + 1,
                  dimnames = list(wwtps, -lag_max:lag_max))

for (w in wwtps) {
  sub <- weekly_wwtp %>% filter(WWTP_Name == w) %>% arrange(weekly_number)
  if (nrow(sub) > lag_max + 2) {
    cr <- ccf(sub$weekly_gene, sub$weekly_cases, lag.max = lag_max, plot = FALSE)
    ccf_mat[w, as.character(cr$lag)] <- cr$acf
  }
}

ccf_df <- as.data.frame(as.table(ccf_mat), stringsAsFactors = FALSE)
colnames(ccf_df) <- c("WWTP", "lag", "correlation")
ccf_df$lag <- as.numeric(ccf_df$lag)
ccf_df$WWTP <- factor(ccf_df$WWTP, levels = wwtps)

best_lag <- ccf_df %>% group_by(WWTP) %>% slice_max(abs(correlation), n = 1)
ccf_df$is_best <- interaction(ccf_df$WWTP, ccf_df$lag) %in%
  interaction(best_lag$WWTP, best_lag$lag)

p_b <- ggplot(ccf_df, aes(x = lag, y = WWTP, fill = correlation)) +
  geom_tile(color = "grey70") +
  geom_tile(data = filter(ccf_df, is_best), fill = NA, color = "black", linewidth = 0.8) +
  geom_text(aes(label = ifelse(is_best, sprintf("%.2f*", correlation),
                               ifelse(abs(correlation) > 0.3, sprintf("%.2f", correlation), ""))),
            size = 2.3) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, limits = c(-1, 1), name = "Correlation") +
  labs(x = "Lag (weeks)", y = NULL, tag = "(b)") +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.y = element_text(size = 7),
    plot.tag = element_text(size = 12, margin = margin(t = 10)),
    plot.tag.position = "bottom",
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 7)
  )

# ---- Figure (c): Temporal Comparison ----
cat("=== Figure (c): Temporal Comparison ===\n")

p_c <- ggplot(temporal_data) +
  geom_bar(aes(x = report_date, y = total_daily_cases / 1000),
           stat = "identity", fill = "steelblue", alpha = 0.3, width = 1) +
  geom_line(aes(x = report_date, y = cases_ma7 / 1000), color = "blue", linewidth = 0.6) +
  geom_line(aes(x = report_date, y = gene_ma7 * 100), color = "red", linewidth = 0.6) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "2 months") +
  scale_y_continuous(name = "Cases (thousands)",
                     sec.axis = sec_axis(~ . / 100, name = "Gene/mL (×100)")) +
  labs(x = NULL, tag = "(c)") +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    plot.tag = element_text(size = 12, margin = margin(t = 10)),
    plot.tag.position = "bottom",
    axis.title.y.left = element_text(size = 8),
    axis.title.y.right = element_text(size = 8)
  )

# ---- Figure (d): Gene Signal by WWTP ----
cat("=== Figure (d): Gene Signal by WWTP ===\n")

gene_ts <- analysis_data %>%
  group_by(report_date, WWTP_Name) %>%
  summarise(gene_mean = mean(gene_copies_ml, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(gene_mean))

p_d <- ggplot(gene_ts) +
  geom_line(aes(x = report_date, y = gene_mean, color = WWTP_Name),
            linewidth = 0.5, alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~WWTP_Name, scales = "free_y", ncol = 2) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
  labs(x = "Date", y = "Gene Copies/mL", tag = "(d)") +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    strip.text = element_text(face = "bold", size = 7),
    plot.tag = element_text(size = 12, margin = margin(t = 10)),
    plot.tag.position = "bottom",
    panel.spacing = unit(0.3, "lines")
  )

# ---- Combine all 4 figures ----
library(patchwork)

combined_fig <- (p_a + p_b) / (p_c + p_d) +
  plot_layout(guides = "collect") &
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank()
  )

ggsave(
  filename = file.path(fig_dir, "fig9_ccf_combined.png"),
  plot = combined_fig,
  width = 16,
  height = 14,
  dpi = 600,
  bg = "white"
)

cat("  Saved: fig9_ccf_combined.png\n")

cat("\n=== SCRIPT 3C EXTRAS COMPLETE ===\n")





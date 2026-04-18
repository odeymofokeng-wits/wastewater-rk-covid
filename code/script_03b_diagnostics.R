##############################################################################
# SCRIPT 3B: EDA DIAGNOSTICS (Part 2 of 3) — Figures 7 to 12
#
# Figures:
#   7a-c. Moran's I scatter plots
#   8.  Empirical variograms (3-panel)
#   9.  Correlation heatmap
#   10. Gene vs Cases scatter by WWTP
#   11. Temporal comparison (cases + gene)
#   12. Gene signal by WWTP
#
# REQUIRES: master_analysis.RData + script3a_map_data.RData
# SAVES:    script3b_temporal.RData (loaded by Script 3C)
##############################################################################

# ========================== LIBRARIES ==========================
library(sf); library(dplyr); library(ggplot2)
library(spdep); library(gstat); library(sp)
library(corrplot); library(RColorBrewer); library(viridis)
library(scales); library(lubridate); library(zoo)

library(dplyr)
library(tidyr)
library(conflicted)

# Force R to always use dplyr versions of common functions
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("arrange", "dplyr")
conflict_prefer("rename", "dplyr")

# ========================== PATHS ==========================
output_dir <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/EDA_Output"
fig_dir <- file.path(output_dir, "figures")

# ========================== LOAD ==========================
cat("=== Loading Data ===\n")
load(file.path(output_dir, "master_analysis.RData"))
load(file.path(output_dir, "script3a_map_data.RData"))

# ========================== FIGURE 7: Moran's I ==========================
cat("=== Figure 7: Spatial Autocorrelation ===\n")

sf_map_proj <- sf_map %>% st_transform(crs = 32736)
nb <- poly2nb(sf_map_proj, queen = TRUE)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
cat("  Neighbours:", length(nb), "subplaces\n")

moran_cases <- moran.test(sf_map$total_cases, lw, 
                          zero.policy = TRUE, 
                          na.action = na.omit)
moran_gene  <- moran.test(sf_map$mean_gene, lw, zero.policy = TRUE, na.action = na.exclude)
moran_cp1k  <- moran.test(sf_map$mean_cases_per_1000, lw, zero.policy = TRUE, na.action = na.exclude)

cat("  Moran I (total cases):", round(moran_cases$estimate[1], 4), "p:", format(moran_cases$p.value, digits = 4), "\n")
cat("  Moran I (mean gene):",   round(moran_gene$estimate[1], 4),  "p:", format(moran_gene$p.value, digits = 4), "\n")
cat("  Moran I (cases/1000):",  round(moran_cp1k$estimate[1], 4),  "p:", format(moran_cp1k$p.value, digits = 4), "\n")

# Create valid index and subset for plotting
valid_idx <- !is.na(sf_map$total_cases)
cases_valid <- sf_map$total_cases[valid_idx]
nb_valid <- subset(nb, valid_idx)
lw_valid <- nb2listw(nb_valid, style = "W", zero.policy = TRUE)

valid_gene <- !is.na(sf_map$mean_gene)
gene_valid <- sf_map$mean_gene[valid_gene]
nb_gene <- subset(nb, valid_gene)
lw_gene <- nb2listw(nb_gene, style = "W", zero.policy = TRUE)

valid_cp1k <- !is.na(sf_map$mean_cases_per_1000)
cp1k_valid <- sf_map$mean_cases_per_1000[valid_cp1k]
nb_cp1k <- subset(nb, valid_cp1k)
lw_cp1k <- nb2listw(nb_cp1k, style = "W", zero.policy = TRUE)

# ================= MORAN SCATTER DATA =================

moran_df_cases <- data.frame(
  value = cases_valid,
  lag   = lag.listw(lw_valid, cases_valid, zero.policy = TRUE)
)

moran_df_gene <- data.frame(
  value = gene_valid,
  lag   = lag.listw(lw_gene, gene_valid, zero.policy = TRUE)
)

moran_df_cp1k <- data.frame(
  value = cp1k_valid,
  lag   = lag.listw(lw_cp1k, cp1k_valid, zero.policy = TRUE)
)

# Custom Moran plot function with quadrant lines
create_moran_plot <- function(data, x_label, y_label, point_color, tag_label) {
  
  # Calculate means for quadrant lines
  x_mean <- mean(data$value, na.rm = TRUE)
  y_mean <- mean(data$lag, na.rm = TRUE)
  
  # Calculate Moran's I for slope line
  moran_i <- cor(data$value, data$lag, use = "complete.obs")
  
  ggplot(data, aes(x = value, y = lag)) +
    # Add quadrant lines (dotted cross)
    geom_vline(xintercept = x_mean, linetype = "dotted", color = "gray40", linewidth = 0.8) +
    geom_hline(yintercept = y_mean, linetype = "dotted", color = "gray40", linewidth = 0.8) +
    # Points
    geom_point(color = point_color, size = 2, alpha = 0.7, stroke = 0.5) +
    # Regression line (slope = Moran's I)
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
    # Labels
    labs(
      x = x_label,
      y = y_label,
      tag = tag_label
    ) +
    # Clean theme
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.tag = element_text(size = 14, face = "plain", margin = margin(t = 10)),
      plot.tag.position = "bottom",
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 20, 10)
    )
}

# Create individual plots
p_a <- create_moran_plot(
  moran_df_cases, 
  "Total COVID-19 Cases", 
  "Spatial Lag of Cases", 
  "steelblue", 
  "(a)"
)

p_b <- create_moran_plot(
  moran_df_gene, 
  "Mean Gene Copies per mL", 
  "Spatial Lag of Gene Copies", 
  "darkgreen", 
  "(b)"
)

p_c <- create_moran_plot(
  moran_df_cp1k, 
  "Cases per 1 000 Population", 
  "Spatial Lag of Cases per 1 000", 
  "firebrick", 
  "(c)"
)

# Combine plots using patchwork or cowplot
library(patchwork)

combined_plot <- p_a + p_b + p_c + 
  plot_layout(ncol = 3) &
  theme(plot.margin = margin(10, 10, 20, 10))


# Save combined figure
ggsave(
  filename = file.path(fig_dir, "fig7_moran_combined.png"),
  plot = combined_plot,
  width = 18,
  height = 6,
  dpi = 600,
  bg = "white"
)

cat("  Saved: fig7_moran_combined.png\n")

# ========================== FIGURE 8: Variograms ==========================
cat("=== Figure 8: Empirical Variograms ===\n")

coords <- st_coordinates(st_centroid(sf_map_proj))
warning("st_centroid assumes attributes are constant over geometries")

# FIX: Remove NAs before creating SpatialPointsDataFrame
valid_idx <- !is.na(sf_map$total_cases) & !is.na(sf_map$mean_cases_per_1000)
coords_valid <- coords[valid_idx, ]
data_valid <- data.frame(
  total_cases = sf_map$total_cases[valid_idx],
  mean_gene = sf_map$mean_gene[valid_idx],
  mean_cp1k = sf_map$mean_cases_per_1000[valid_idx]
)

cat("  Valid observations for variogram:", sum(valid_idx), "of", nrow(sf_map), "\n")

spdf <- SpatialPointsDataFrame(coords = coords_valid, data = data_valid)

# 8a: Variogram for total cases
vg_cases <- variogram(total_cases ~ 1, data = spdf)

# 8b: Variogram for gene (use same valid_idx to ensure alignment)
valid_gene <- !is.na(sf_map$mean_gene)
cat("  Valid gene observations:", sum(valid_gene), "of", nrow(sf_map), "\n")

if(sum(valid_gene) > 10) {
  coords_gene <- coords[valid_gene, ]
  spdf_gene <- SpatialPointsDataFrame(
    coords = coords_gene,
    data = data.frame(mean_gene = sf_map$mean_gene[valid_gene])
  )
  vg_gene <- variogram(mean_gene ~ 1, data = spdf_gene)
  cat("  Gene variogram computed with", nrow(vg_gene), "bins\n")
} else {
  cat("  WARNING: Insufficient gene data for variogram\n")
  vg_gene <- NULL
}

# 8c: Variogram for cases per 1000
vg_cp1k <- variogram(mean_cp1k ~ 1, data = spdf)

# ================= GGPLOT VARIOGRAMS WITH PATCHWORK =================

library(ggplot2)

# Convert variogram objects to data frames for ggplot
vg_cases_df <- data.frame(
  dist = vg_cases$dist,
  gamma = vg_cases$gamma,
  np = vg_cases$np
)

vg_cp1k_df <- data.frame(
  dist = vg_cp1k$dist,
  gamma = vg_cp1k$gamma,
  np = vg_cp1k$np
)

# Create base plot function
create_variogram_plot <- function(vg_df, x_label, y_label, point_color, tag_label) {
  ggplot(vg_df, aes(x = dist, y = gamma)) +
    geom_point(color = point_color, size = 2.5, alpha = 0.8) +
    geom_line(color = point_color, linewidth = 0.5, alpha = 0.5) +
    labs(
      x = x_label,
      y = y_label,
      tag = tag_label
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.tag = element_text(size = 14, margin = margin(t = 10)),
      plot.tag.position = "bottom",
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 20, 10)
    )
}

p_8a <- create_variogram_plot(
  vg_cases_df,
  "Distance (meters)",
  "Semivariance (Total Cases)",
  "steelblue",
  "(a)"
)

p_8c <- create_variogram_plot(
  vg_cp1k_df,
  "Distance (meters)",
  "Semivariance (Cases per 1,000)",
  "firebrick",
  "(c)"
)

# Handle gene plot
if(!is.null(vg_gene)) {
  vg_gene_df <- data.frame(
    dist = vg_gene$dist,
    gamma = vg_gene$gamma,
    np = vg_gene$np
  )
  
  p_8b <- create_variogram_plot(
    vg_gene_df,
    "Distance (meters)",
    "Semivariance (Gene Copies per mL)",
    "darkgreen",
    "(b)"
  )
  
  combined_vg <- p_8a + p_8b + p_8c + plot_layout(ncol = 3)
} else {
  # If no gene data, just show cases and cp1k with empty placeholder or skip
  combined_vg <- p_8a + p_8c + plot_layout(ncol = 2)
  cat("  Combined plot: 2 panels only (gene data insufficient)\n")
}

ggsave(
  filename = file.path(fig_dir, "fig8_variograms_combined.png"),
  plot = combined_vg,
  width = 18,
  height = 6,
  dpi = 600,
  bg = "white"
)

cat("  Saved: fig8_variograms_combined.png\n")

# ========================== FIGURE 9: Correlation Heatmap ==========================
cat("=== Figure 9: Correlation Heatmap ===\n")

corr_vars <- analysis_data %>%
  select(daily_cases_sp, gene_copies_ml, cases_per_1000,
         wwtp_dist_km, hosp_dist_km,
         COVID_Vulnerability_Index, Health_Service_Vulnerability_Index,
         Total_Population_Vulnerability_Index, Estimated_Population_2018) %>%
  na.omit()

cat("  Correlation data:", nrow(corr_vars), "rows\n")

if (nrow(corr_vars) > 50) {
  cor_mat <- cor(corr_vars, use = "pairwise.complete.obs")
  short_names <- c("Daily Cases", "Gene/mL", "Cases/1000", "WWTP Dist",
                   "Hosp Dist", "COVID Vuln.", "Health Vuln.", "Pop Vuln.", "Pop")
  colnames(cor_mat) <- short_names; rownames(cor_mat) <- short_names
  
  cat("  Correlations with daily cases:\n")
  for (v in short_names[-1]) cat(sprintf("    %s: r = %.3f\n", v, cor_mat["Daily Cases", v]))
  
  # FIX: Use pixels for dimensions, fix color palette, reduce margins
  png(file.path(fig_dir, "fig9_correlation_heatmap.png"), width = 2400, height = 2100, res = 600)
  
  # FIX: Use colorRampPalette to expand RdBu to 200 colors, or use fewer colors
  col_palette <- colorRampPalette(brewer.pal(11, "RdBu"))(200)
  
  # FIX: Reduce margins and title size
  par(mar = c(2, 2, 2, 2))
  corrplot(cor_mat, method = "color", type = "upper", tl.col = "black", tl.srt = 45,
           tl.cex = 0.7, col = col_palette,
           addCoef.col = "black", number.cex = 0.5,
           cl.cex = 0.7)
  dev.off()
  cat("  Saved: fig9_correlation_heatmap\n")
}

# ========================== FIGURE 10: Gene vs Cases ==========================
cat("=== Figure 10: Gene vs Cases Scatter ===\n")

gc_data <- analysis_data %>%
  filter(!is.na(gene_copies_ml), gene_copies_ml > 0, daily_cases_sp >= 0)

cat("  Gene-cases points:", nrow(gc_data), "\n")

if (nrow(gc_data) > 10) {
  p10 <- ggplot(gc_data) +
    geom_point(aes(x = log10(gene_copies_ml + 1), y = log10(daily_cases_sp + 1),
                   color = WWTP_Name), alpha = 0.3, size = 1) +
    geom_smooth(aes(x = log10(gene_copies_ml + 1), y = log10(daily_cases_sp + 1)),
                method = "lm", color = "red", se = TRUE, linewidth = 0.8) +
    facet_wrap(~WWTP_Name, scales = "free") +
    labs(x = expression(log[10](Gene~Copies/mL + 1)),
         y = expression(log[10](Daily~Cases + 1))) +
    theme_minimal(base_size = 10) + theme(legend.position = "none")

  ggsave(file.path(fig_dir, "fig10_gene_vs_cases.png"), p10, width = 12, height = 8, dpi = 600)
  
  cat("  Saved: fig10_gene_vs_cases\n")
}

# ========================== FIGURE 11: Temporal Comparison ==========================
cat("=== Figure 11: Temporal Comparison ===\n")

temporal_data <- analysis_data %>%
  group_by(report_date) %>%
  summarise(
    total_daily_cases = sum(daily_cases_sp, na.rm = TRUE),
    mean_gene = mean(gene_copies_ml, na.rm = TRUE),
    n_gene_obs = sum(!is.na(gene_copies_ml)),
    .groups = "drop"
  ) %>%
  mutate(
    cases_ma7 = rollmeanr(total_daily_cases, k = 7, fill = NA),
    gene_ma7  = rollmeanr(mean_gene, k = 7, fill = NA)
  )

cat("  Temporal data:", nrow(temporal_data), "days\n")

p11 <- ggplot(temporal_data) +
  geom_bar(aes(x = report_date, y = total_daily_cases / 1000),
           stat = "identity", fill = "steelblue", alpha = 0.3, width = 1) +
  geom_line(aes(x = report_date, y = cases_ma7 / 1000), color = "blue", linewidth = 0.8) +
  geom_line(aes(x = report_date, y = gene_ma7 * 100), color = "red", linewidth = 0.8) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "2 months") +
  scale_y_continuous(name = "COVID-19 Cases (thousands)",
                    sec.axis = sec_axis(~ . / 100, name = "Gene Copies/mL (scaled)")) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "fig11_temporal_comparison.png"), p11, width = 14, height = 6, dpi = 600)

cat("  Saved: fig11_temporal_comparison\n")

# ========================== FIGURE 12: Gene by WWTP ==========================
cat("=== Figure 12: Gene Signal by WWTP ===\n")

gene_ts <- analysis_data %>%
  filter(!is.na(gene_copies_ml)) %>%
  group_by(WWTP_Name, report_date) %>%
  summarise(gene_mean = mean(gene_copies_ml, na.rm = TRUE), .groups = "drop")

p12 <- ggplot(gene_ts) +
  geom_line(aes(x = report_date, y = gene_mean, color = WWTP_Name),
            linewidth = 0.6, alpha = 0.8) +
  facet_wrap(~WWTP_Name, scales = "free_y", ncol = 2) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
  labs(
       x = "Date", y = "Gene Copies/mL") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "none", strip.text = element_text(face = "bold"))

ggsave(file.path(fig_dir, "fig12_gene_by_wwtp.png"), p12, width = 12, height = 10, dpi = 600)

cat("  Saved: fig12_gene_by_wwtp\n")

# ========================== SAVE ==========================
save(temporal_data, sp_totals,
     file = file.path(output_dir, "script3b_temporal.RData"))

cat("\n=== SCRIPT 3B COMPLETE ===\n")


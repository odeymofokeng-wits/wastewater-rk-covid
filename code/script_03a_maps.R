##############################################################################
# SCRIPT 3A: EDA MAPS (Part 1 of 3) — Figures 1 to 6
#
# Figures:
#   1. Study area map (WWTP service areas)
#   2. COVID-19 cases choropleth
#   3. Vulnerability index map
#   4. Hospital accessibility map
#   5. WWTP accessibility map
#   6. Gene data coverage map
#
# REQUIRES: master_analysis.RData (from Script 2B)
# SAVES:    script3a_map_data.RData (loaded by Script 3B)
##############################################################################

# ========================== LIBRARIES ==========================
library(sf); library(dplyr); library(ggplot2)
library(ggspatial); library(viridis); library(scales); library(lubridate)


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
path_subplaces    <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/WASTEWATER DATA/Ekurhuleni sub and main places/Eku_sub_places13.shp"
path_wwtp         <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/WASTEWATER DATA/Ekurhuleni WWTP/Ekurhuleni_WWTP.shp"
path_wwtp_service <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/WASTEWATER DATA/WWTP Serviced Areas-20260328T133206Z-1-001/WWTP Serviced Areas/Ekurhuleni_WWTP_Service_Areas.shp"
output_dir <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/EDA_Output"
fig_dir <- file.path(output_dir, "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# ========================== LOAD ==========================
cat("=== Loading Data ===\n")
load(file.path(output_dir, "master_analysis.RData"))
sf_subplaces    <- st_read(path_subplaces, quiet = TRUE)
sf_wwtp         <- st_read(path_wwtp, quiet = TRUE)

# Per-subplace aggregates
sp_totals <- analysis_data %>%
  group_by(SP_CODE) %>%
  summarise(
    total_cases = sum(daily_cases_sp, na.rm = TRUE),
    mean_cases_per_1000 = mean(cases_per_1000, na.rm = TRUE),
    mean_gene = mean(gene_copies_ml, na.rm = TRUE),
    n_gene_obs = sum(!is.na(gene_copies_ml)),
    pop = first(Estimated_Population_2018),
    vulnerability = first(Composite_Vulnerability_Score),
    covid_vuln = first(COVID_Vulnerability_Index),
    hosp_dist = first(hosp_dist_km),
    wwtp_dist = first(wwtp_dist_km),
    .groups = "drop"
  )

sf_map <- sf_subplaces %>%
  left_join(sp_totals, by = "SP_CODE") %>%
  left_join(wwtp_sp_mapping %>% select(SP_CODE, WWTP_Name), by = "SP_CODE")

cat("  Map data:", nrow(sf_map), "subplaces\n")

# ========================== FIGURE 1: Study Area ==========================
cat("=== Figure 1: Study Area ===\n")

sf_map_proj  <- st_transform(sf_map, 32735)
sf_wwtp_proj <- st_transform(sf_wwtp, 32735)


library(ggplot2)
library(sf)
library(ggspatial)
library(cowplot)  # for ggdraw
library(grid)

# Read the provincial boundaries
sa_provinces <- st_read("C:/Users/odeyr/OneDrive - University of Witwatersrand/2026 Mathematical Statistics Masters Thesis/KZN Data Sources/Socio-Economic Data/zaf_admin_boundaries_shapefiles/zaf_admin1.shp")

# Ensure same CRS
sa_provinces_proj <- st_transform(sa_provinces, st_crs(sf_map_proj))
names(sa_provinces)


# Create the inset map - showing SA with Gauteng highlighted
inset_map <- ggplot(sa_provinces_proj) +
  geom_sf(aes(fill = adm1_name == "Gauteng"), color = "grey40", linewidth = 0.2) +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "grey90"), 
                    guide = "none") +
  # Add a point or box showing where the main map is
  geom_sf(data = st_union(sf_map_proj), fill = NA, color = "black", linewidth = 0.8) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    plot.margin = margin(2, 2, 2, 2)
  )





p1 <- ggplot(sf_map_proj) +
  geom_sf(aes(fill = WWTP_Name), color = "grey40", linewidth = 0.1) +
  geom_sf(data = sf_wwtp_proj, shape = 17, size = 3, color = "red", fill = "red") +
  geom_sf_text(data = sf_wwtp_proj, aes(label = Water_Care), size = 2.5,
               nudge_x = 1000, nudge_y = 1000, fontface = "bold") +
  scale_fill_brewer(palette = "Set2", name = "WWTP Service Area", na.value = "grey80") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),
    plot.title = element_text(face = "bold"),
    
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    
    panel.background = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         height = unit(1, "cm"), width = unit(0.5, "cm"))

# Combine using cowplot
combined_plot <- ggdraw(p1) +
  draw_plot(inset_map, 
            x = -0.01, y = 0.61,   # <-- higher up (y=0.72), more to left (x=0.02)
            width = 0.345, height = 0.35)  # <-- smaller size so it fits


ggsave(file.path(fig_dir, "fig1_study_area.png"), combined_plot, width = 12, height = 8, dpi = 600)
cat("  Saved: fig1_study_area\n")

# ========================== FIGURE 2: COVID Cases ==========================
cat("=== Figure 2: COVID Cases Choropleth ===\n")

library(RColorBrewer)

p2 <- ggplot(sf_map) +
  geom_sf(aes(fill = cut_number(mean_cases_per_1000, 5,
                                labels = paste("Range", 1:5))),
          color = "grey50", linewidth = 0.05) +
  geom_sf(data = sf_wwtp, shape = 17, size = 2, color = "black") +
  scale_fill_brewer(palette = "YlOrRd", name = "Mean Cases\nper 1000",
                    na.value = "grey80") +
  coord_sf(expand = FALSE) +   
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    
    
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    
    
    panel.background = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )


ggsave(file.path(fig_dir, "fig2_covid_cases_map.png"), p2, width = 10, height = 7, dpi = 600)
cat("  Saved: fig2_covid_cases_map\n")

# ========================== FIGURE 3: Vulnerability ==========================
cat("=== Figure 3: Vulnerability Index ===\n")

p3 <- ggplot(sf_map) +
  geom_sf(aes(fill = covid_vuln), color = "grey50", linewidth = 0.05) +
  scale_fill_viridis_c(
    name = "COVID-19\nVulnerability\nIndex",
    option = "magma",
    na.value = "grey80"
  ) +
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(fig_dir, "fig3_vulnerability_map.png"), p3, width = 10, height = 7, dpi = 600)

cat("  Saved: fig3_vulnerability_map\n")

# ========================== FIGURE 4: Hospital Distance ==========================
cat("=== Figure 4: Hospital Distance ===\n")

p4 <- ggplot(sf_map) +
  geom_sf(aes(fill = hosp_dist), color = "grey50", linewidth = 0.05) +
  geom_sf(data = sf_wwtp, shape = 17, size = 2, color = "blue") +
  scale_fill_viridis_c(
    name = "Distance to\nNearest Hospital\n(km)",
    option = "plasma"
  ) +
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(fig_dir, "fig4_hospital_distance.png"), p4, width = 10, height = 7, dpi = 600)

cat("  Saved: fig4_hospital_distance\n")

# ========================== FIGURE 5: WWTP Distance ==========================
cat("=== Figure 5: WWTP Distance ===\n")

p5 <- ggplot(sf_map) +
  geom_sf(aes(fill = wwtp_dist), color = "grey50", linewidth = 0.05) +
  geom_sf(data = sf_wwtp, shape = 17, size = 2, color = "red") +
  scale_fill_viridis_c(
    name = "Distance to\nNearest WWTP\n(km)",
    option = "inferno"
  ) +
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(fig_dir, "fig5_wwtp_distance.png"), p5, width = 10, height = 7, dpi = 600)

cat("  Saved: fig5_wwtp_distance\n")

# ========================== FIGURE 6: Gene Coverage ==========================
cat("=== Figure 6: Gene Data Coverage ===\n")

p6 <- ggplot(sf_map) +
  geom_sf(aes(fill = n_gene_obs), color = "grey50", linewidth = 0.05) +
  geom_sf(data = sf_wwtp, shape = 17, size = 2, color = "red") +
  scale_fill_viridis_c(
    name = "Number of\nGene Observations",
    option = "cividis",
    na.value = "grey80"
  ) +
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(fig_dir, "fig6_gene_coverage.png"), p6, width = 10, height = 7, dpi = 600)

cat("  Saved: fig6_gene_coverage\n")

combined <- (p4 + labs(tag = "(a)") |
               p5 + labs(tag = "(b)")) +
  plot_layout(ncol = 2) &
  theme(
    plot.tag = element_text(size = 11),
    plot.tag.position = "bottom"
  )

ggsave(
  file.path(fig_dir, "fig4_5_combined.png"),
  combined,
  width = 14,
  height = 7,
  dpi = 600
)

# ========================== SAVE FOR NEXT SCRIPT ==========================
save(sf_map, sf_subplaces, sf_wwtp, sp_totals,
     file = file.path(output_dir, "script3a_map_data.RData"))

cat("\n=== SCRIPT 3A COMPLETE ===\n")





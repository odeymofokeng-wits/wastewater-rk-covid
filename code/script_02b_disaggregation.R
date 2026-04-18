##############################################################################
# SCRIPT 2B: SPATIAL DISAGGREGATION (Part 2 of 2)
#
# This script:
#   Step 4: Compute distances (WWTP, hospital)
#   Step 5: Disaggregate COVID cases to subplaces (population weights)
#   Step 6: Disaggregate gene counts to subplaces (WWTP assignment)
#   Step 7: Create lagged variables
#   Step 8: Final cleaning and validation
#   Step 9: Save master_analysis.RData
#
# REQUIRES: Script 1 (cleaned_data.RData) + Script 2A (script2a_mappings.RData)
##############################################################################

# ========================== LIBRARIES ==========================
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(zoo)

# ========================== PATHS ==========================
path_subplaces <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/WASTEWATER DATA/Ekurhuleni sub and main places/Eku_sub_places13.shp"
path_wwtp      <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/WASTEWATER DATA/Ekurhuleni WWTP/Ekurhuleni_WWTP.shp"
data_dir   <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/WASTEWATER DATA/Cleaned datasets"
output_dir <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/EDA_Output"

# ========================== LOAD ==========================
cat("=== Loading Data ===\n")
load(file.path(output_dir, "cleaned_data.RData"))
load(file.path(output_dir, "script2a_mappings.RData"))

sf_subplaces <- st_read(path_subplaces, quiet = TRUE)
sf_wwtp      <- st_read(path_wwtp, quiet = TRUE)

cat("  Cleaned data loaded\n")
cat("  Mappings loaded: WWTP=", nrow(wwtp_sp_mapping),
    ", Ward-SD=", nrow(ward_sd_full),
    ", SP-SD=", nrow(sp_sd_assignment), "\n")

# ========================== STEP 4: Compute Distances ==========================
cat("\n=== Step 4: Computing Distances ===\n")

sf_sub_proj <- st_transform(sf_subplaces, crs = 32736)
sf_wwtp_proj <- st_transform(sf_wwtp, crs = 32736)

hospitals_sf_dist <- st_transform(
  st_as_sf(hospitals_enriched, coords = c("Longitude", "Latitude"), crs = 4326),
  crs = 32736
)

sp_centroids <- st_centroid(sf_sub_proj)

dmat_wwtp <- st_distance(sp_centroids, sf_wwtp_proj)
wwtp_dist <- apply(dmat_wwtp, 1, min, na.rm = TRUE) / 1000

dmat_hosp <- st_distance(sp_centroids, hospitals_sf_dist)
hosp_dist <- apply(dmat_hosp, 1, min, na.rm = TRUE) / 1000

cat("  WWTP dist: mean =", round(mean(wwtp_dist, na.rm = TRUE), 2), "km\n")
cat("  Hosp dist: mean =", round(mean(hosp_dist, na.rm = TRUE), 2), "km\n")

# ========================== STEP 5: Disaggregate COVID Cases ==========================
cat("\n=== Step 5: Disaggregating COVID Cases ===\n")

sp_master <- sf_subplaces %>%
  st_drop_geometry() %>%
  left_join(
    wwtp_sp_mapping %>% dplyr::select(SP_CODE, WWTP_Name),
    by = "SP_CODE"
  ) %>%
  left_join(
    sp_sd_assignment %>% dplyr::select(SP_CODE, health_subdistrict, n_wards),
    by = "SP_CODE"
  ) %>%
  left_join(
    vulnerability %>% dplyr::select(
      Subplace_Code,
      Estimated_Population_2018,
      COVID_Vulnerability_Index,
      Health_Service_Vulnerability_Index,
      Total_Population_Vulnerability_Index,
      Composite_Vulnerability_Score,
      Is_Zero_Population
    ),
    by = c("SP_CODE" = "Subplace_Code")
  ) %>%
  dplyr::mutate(
    wwtp_dist_km = wwtp_dist,
    hosp_dist_km = hosp_dist
  )

cat("  Master subplace table:", nrow(sp_master), "rows\n")
cat("  WWTP coverage:", sum(!is.na(sp_master$WWTP_Name)), "of", nrow(sp_master), "\n")
cat("  Health SD coverage:", sum(!is.na(sp_master$health_subdistrict)), "of", nrow(sp_master), "\n")

# Population share within each health sub-district
sp_master <- sp_master %>%
  group_by(health_subdistrict) %>%
  mutate(
    sd_total_pop = sum(Estimated_Population_2018, na.rm = TRUE),
    pop_share_sd = Estimated_Population_2018 / sd_total_pop
  ) %>%
  ungroup() %>%
  mutate(pop_share_sd = ifelse(Is_Zero_Population == 1, 0, pop_share_sd))

cat("  Population share check (should sum to ~1 per SD):\n")
print(sp_master %>% group_by(health_subdistrict) %>%
        summarise(total_share = sum(pop_share_sd, na.rm = TRUE), .groups = "drop"))

# Create subplace x date dataset
sd_cols <- c("daily_east1", "daily_east2", "daily_north1", "daily_north2",
             "daily_south1", "daily_south2")

covid_sd <- covid_daily %>%
  dplyr::filter(!is.na(daily_east1)) %>%
  dplyr::select(report_date, dplyr::all_of(sd_cols))
cat("  COVID dates:", nrow(covid_sd), "\n")

covid_long <- covid_sd %>%
  pivot_longer(cols = all_of(sd_cols), names_to = "health_subdistrict",
               values_to = "daily_cases_sd") %>%
  mutate(health_subdistrict = gsub("daily_", "", health_subdistrict))

analysis_data <- sp_master %>%
  dplyr::select(SP_CODE, SP_NAME, health_subdistrict, pop_share_sd,
                Estimated_Population_2018, Is_Zero_Population,
                WWTP_Name, wwtp_dist_km, hosp_dist_km,
                COVID_Vulnerability_Index, Health_Service_Vulnerability_Index,
                Total_Population_Vulnerability_Index, Composite_Vulnerability_Score) %>%
  dplyr::right_join(covid_long, by = "health_subdistrict") %>%
  dplyr::mutate(
    daily_cases_sp = daily_cases_sd * pop_share_sd,
    weekly_number = floor_date(report_date, unit = "week"),
    month = month(report_date),
    year = year(report_date),
    day_of_year = yday(report_date)
  )

# ========================== STEP 6: Disaggregate Gene Counts ==========================
cat("\n=== Step 6: Disaggregating Gene Counts ===\n")

analysis_data$weekly_number <- as.Date(analysis_data$weekly_number)

wwtp_name_map <- data.frame(
  gene_name = c("Daveyton WWTW", "Olifantsfontein WWTW", "Vlakplaats WWTW",
                "Carl Grundlingh WWTW", "Herbert Bickley WWTW", "Jan Smuts WWTW",
                "JP Marais WWTW", "H. Rynfield WWTW"),
  erwat_name = c("Daveyton", "Olifantsfontein", "Vlakplaats",
                 "Carl Grundlingh", "Herbert Bickley", "Jan Smuts",
                 "JP Marais", "Rynfield"),
  stringsAsFactors = FALSE
)

cat("  Unique site_name_standardized in gene_weekly:\n")
print(unique(gene_weekly$site_name_standardized))

gene_wwtp <- gene_weekly %>%
  dplyr::left_join(wwtp_name_map, by = c("site_name_standardized" = "gene_name")) %>%
  dplyr::rename(WWTP_Name = erwat_name) %>%
  dplyr::select(WWTP_Name, week_start, gene_copies_ml, n_samples, pct_positive)

cat("  Gene by WWTP-week:", nrow(gene_wwtp), "| Matched:",
    sum(!is.na(gene_wwtp$WWTP_Name)), "\n")
cat("  Gene per WWTP:\n")
print(gene_wwtp %>%
        dplyr::filter(!is.na(gene_copies_ml)) %>%
        dplyr::group_by(WWTP_Name) %>%
        dplyr::summarise(n = n(), .groups = "drop"))

analysis_data <- analysis_data %>%
  dplyr::left_join(
    gene_wwtp %>% dplyr::rename(weekly_number = week_start),
    by = c("WWTP_Name", "weekly_number")
  )

cat("  After gene merge:", nrow(analysis_data), "rows | Gene coverage:",
    sum(!is.na(analysis_data$gene_copies_ml)),
    sprintf("(%.1f%%)\n", 100 * sum(!is.na(analysis_data$gene_copies_ml)) / nrow(analysis_data)))

# ========================== STEP 7: Lagged Variables ==========================
cat("\n=== Step 7: Creating Lagged Variables ===\n")

analysis_data <- analysis_data %>%
  dplyr::group_by(SP_CODE) %>%
  dplyr::arrange(report_date) %>%
  dplyr::mutate(
    gene_lag1 = dplyr::lag(gene_copies_ml, 1),
    gene_lag2 = dplyr::lag(gene_copies_ml, 7),
    gene_lag3 = dplyr::lag(gene_copies_ml, 14),
    cases_lag1 = dplyr::lag(daily_cases_sp, 1),
    cases_lag7 = dplyr::lag(daily_cases_sp, 7),
    cases_ma7 = zoo::rollmeanr(daily_cases_sp, k = 7, fill = NA),
    cases_ma14 = zoo::rollmeanr(daily_cases_sp, k = 14, fill = NA)
  ) %>%
  dplyr::ungroup()

cat("  gene_lag1 non-NA:", sum(!is.na(analysis_data$gene_lag1)), "\n")
cat("  gene_lag2 non-NA:", sum(!is.na(analysis_data$gene_lag2)), "\n")
cat("  cases_lag1 non-NA:", sum(!is.na(analysis_data$cases_lag1)), "\n")

# ========================== STEP 8: Final Cleaning ==========================
cat("\n=== Step 8: Final Cleaning ===\n")

n_before <- nrow(analysis_data)
analysis_data <- analysis_data %>%
  dplyr::filter(Is_Zero_Population == 0 | is.na(Is_Zero_Population))
cat("  Removed", n_before - nrow(analysis_data), "zero-pop rows\n")

analysis_data <- analysis_data %>%
  dplyr::mutate(
    cases_per_1000 = daily_cases_sp / (Estimated_Population_2018 / 1000),
    cases_per_1000 = dplyr::if_else(is.infinite(cases_per_1000) | is.nan(cases_per_1000),
                                    NA_real_, cases_per_1000)
  )

cat("\n========================================\n")
cat("  FINAL DATASET SUMMARY\n")
cat("========================================\n")
cat("  Rows:", nrow(analysis_data), "\n")
cat("  Subplaces:", length(unique(analysis_data$SP_CODE)), "\n")
cat("  Dates:", length(unique(analysis_data$report_date)), "\n")
cat("  Range:", min(analysis_data$report_date), "to", max(analysis_data$report_date), "\n")
cat("  Columns:", ncol(analysis_data), "\n")
print(names(analysis_data))

cat("\n  Variable summaries:\n")
for (v in c("daily_cases_sp", "gene_copies_ml", "cases_per_1000",
            "wwtp_dist_km", "hosp_dist_km", "COVID_Vulnerability_Index")) {
  if (v %in% names(analysis_data)) {
    cat(sprintf("  %-25s n=%d mean=%.4f med=%.4f\n", v,
                sum(!is.na(analysis_data[[v]])),
                mean(analysis_data[[v]], na.rm = TRUE),
                median(analysis_data[[v]], na.rm = TRUE)))
  }
}

cat("\n  Subplaces per WWTP:\n")
print(table(analysis_data$WWTP_Name, useNA = "ifany"))

# ========================== STEP 9: Save ==========================
cat("\n=== Saving ===\n")

save(
  analysis_data, sp_master, wwtp_sp_mapping, ward_sd_full,
  sp_sd_assignment, wwtp_address_map,
  file = file.path(output_dir, "master_analysis.RData")
)

write.csv(
  analysis_data %>%
    dplyr::filter(report_date == max(report_date)) %>%
    dplyr::select(SP_CODE, SP_NAME, health_subdistrict, WWTP_Name,
                  Estimated_Population_2018, wwtp_dist_km, hosp_dist_km,
                  daily_cases_sp, gene_copies_ml, cases_per_1000),
  file.path(output_dir, "master_analysis_sample.csv"),
  row.names = FALSE
)

cat("  Saved: master_analysis.RData\n")
cat("  Saved: master_analysis_sample.csv\n")
cat("\n=== SCRIPT 2B COMPLETE ===\n")
cat(">>> Run Script 3A next <<<\n")

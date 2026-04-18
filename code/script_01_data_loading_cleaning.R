##############################################################################
# SCRIPT 1: DATA LOADING, CLEANING & MERGING
# Purpose: Load all CSV/Excel data, clean known issues, create lookup tables,
#          and produce a cleaned master RData file.
#
# Prerequisites:
#   - Run Script 0 first 
#   - Install required packages 
#
##############################################################################

# ========================== PACKAGE INSTALLATION ==========================
# Install packages if not already installed
required_packages <- c("sf", "dplyr", "tidyr", "ggplot2", "readxl", 
                       "lubridate", "stringr", "knitr", "kableExtra")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org")
  }
}

# ========================== LIBRARIES ==========================
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(lubridate)
library(stringr)

# ========================== SETUP ==========================
# === PATHS ===
# Update the data_dir to where your CSV/Excel files are stored locally.

data_dir <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/WASTEWATER DATA/Cleaned datasets"

# If using the default download location, change to:
# data_dir <- "C:/Users/odeyr/Downloads"
##############################################################################
# SCRIPT 1: DATA LOADING, CLEANING & MERGING
##############################################################################

# ========================== OUTPUT DIRECTORY ==========================
output_dir <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/EDA_Output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ========================== 1. LOAD VULNERABILITY DATA ==========================
cat("=== 1. Loading Vulnerability Data ===\n")
vulnerability <- read.csv(file.path(data_dir, "Ekurhuleni_Vulnerability_Subplaces_Cleaned.csv"),
                          stringsAsFactors = FALSE)

cat("  Rows:", nrow(vulnerability), "  Columns:", ncol(vulnerability), "\n")
cat("  SP range:", min(vulnerability$Subplace_Code), "-", max(vulnerability$Subplace_Code), "\n")
cat("  Zero-population subplaces:", sum(vulnerability$Is_Zero_Population == 1), "\n")
cat("  NAs in Composite_Vulnerability_Score:", sum(is.na(vulnerability$Composite_Vulnerability_Score)), "\n")

# ========================== 2. LOAD AND CLEAN GENE COUNT DATA ==========================
cat("\n=== 2. Loading & Cleaning Gene Count Data ===\n")
gene_raw <- read_excel(file.path(data_dir, "Gene_Count_data.xlsx"))

cat("  Raw rows:", nrow(gene_raw), "  Columns:", ncol(gene_raw), "\n")
cat("  Sites:", paste(unique(gene_raw$site_name_standardized), collapse = ", "), "\n")
# FIX: Explicitly convert and preserve Date class
gene_data <- gene_raw %>%
  mutate(
    date_sampled = as.Date(date_sampled, format = "%Y-%m-%d"),
    date_received = as.Date(date_received, format = "%Y-%m-%d"),
    date_tested = as.Date(date_tested, format = "%Y-%m-%d")
  ) %>%
  filter(!is.na(date_sampled))

# FIX: Format dates properly for display
cat("  Date range:", format(min(gene_data$date_sampled), "%Y-%m-%d"), 
    "to", format(max(gene_data$date_sampled), "%Y-%m-%d"), "\n")
cat("  After cleaning:", nrow(gene_data), "rows\n")


cat("  Date range:", min(gene_raw$date_sampled), "to", max(gene_raw$date_sampled), "\n")

gene_data <- gene_data %>%
  mutate(
    n1_gene_ct = ifelse(n1_gene_detected == "N", NA_real_, n1_gene_ct),
    n2_gene_ct = ifelse(grepl("^[Nn]$", as.character(`N2 gene`)), NA_real_, n2_gene_ct)
  )

gene_data <- gene_data %>%
  mutate(
    n2_gene_ct = ifelse(n2_gene_ct > 45, NA_real_, n2_gene_ct)
  )

gene_data <- gene_data %>%
  mutate(
    gene_copies_ml = n1_copies_ml,
    gene_copies_avg = rowMeans(cbind(n1_copies_ml, n2_copies_ml), na.rm = TRUE),
    gene_both_available = (!is.na(n1_copies_ml) & !is.na(n2_copies_ml))
  )

cat("  After cleaning:", nrow(gene_data), "rows\n")
cat("  N1 copies/mL: non-NA =", sum(!is.na(gene_data$gene_copies_ml)), 
    "  range:", range(gene_data$gene_copies_ml, na.rm = TRUE), "\n")
cat("  Average copies/mL: non-NA =", sum(!is.na(gene_data$gene_copies_avg)), "\n")
cat("  Both N1 & N2 available:", sum(gene_data$gene_both_available), "\n")

# Site-level summary
site_summary <- gene_data %>%
  group_by(site_name_standardized) %>%
  summarise(
    n_samples = n(),
    date_range = paste(min(date_sampled), "to", max(date_sampled)),
    median_gene_copies = median(gene_copies_ml, na.rm = TRUE),
    pct_positive = mean(is_positive == 1) * 100,
    .groups = "drop"
  )
cat("\n  Site summary:\n")
print(site_summary)

# Weekly aggregation
gene_weekly <- gene_data %>%
  mutate(week_start = floor_date(date_sampled, unit = "week")) %>%
  group_by(site_name_standardized, week_start) %>%
  summarise(
    gene_copies_ml = mean(gene_copies_ml, na.rm = TRUE),
    n_samples = n(),
    pct_positive = mean(is_positive == 1) * 100,
    .groups = "drop"
  )

# ========================== 3. LOAD AND CLEAN COVID DATA ==========================
cat("\n=== 3. Loading & Cleaning COVID Data ===\n")
covid_raw <- read.csv(file.path(data_dir, "ekurhuleni_covid_cleaned.csv"),
                      stringsAsFactors = FALSE)

cat("  Raw rows:", nrow(covid_raw), "  Columns:", ncol(covid_raw), "\n")

covid_data <- covid_raw %>%
  mutate(report_date = as.Date(report_date)) %>%
  filter(!is.na(report_date)) %>%
  arrange(report_date)

covid_data <- covid_data %>%
  distinct(report_date, .keep_all = TRUE)

cat("  After dedup:", nrow(covid_data), "rows\n")
cat("  Date range:", min(covid_data$report_date), "to", max(covid_data$report_date), "\n")

covid_daily <- covid_data %>%
  mutate(
    daily_total = c(NA, diff(eku_total_cases)),
    daily_east1 = c(NA, diff(east1_cases)),
    daily_east2 = c(NA, diff(east2_cases)),
    daily_north1 = c(NA, diff(north1_cases)),
    daily_north2 = c(NA, diff(north2_cases)),
    daily_south1 = c(NA, diff(south1_cases)),
    daily_south2 = c(NA, diff(south2_cases)),
    daily_deaths = c(NA, diff(eku_total_deaths))
  )

daily_cols <- c("daily_total", "daily_east1", "daily_east2", "daily_north1",
                "daily_north2", "daily_south1", "daily_south2", "daily_deaths")

covid_daily[daily_cols] <- lapply(covid_daily[daily_cols], function(x) {
  x[is.na(x)] <- 0
  pmax(x, 0)
})

covid_daily <- covid_daily %>%
  mutate(
    ma7_total = zoo::rollmeanr(daily_total, k = 7, fill = NA),
    ma7_east1 = zoo::rollmeanr(daily_east1, k = 7, fill = NA),
    ma7_east2 = zoo::rollmeanr(daily_east2, k = 7, fill = NA),
    ma7_north1 = zoo::rollmeanr(daily_north1, k = 7, fill = NA),
    ma7_north2 = zoo::rollmeanr(daily_north2, k = 7, fill = NA),
    ma7_south1 = zoo::rollmeanr(daily_south1, k = 7, fill = NA),
    ma7_south2 = zoo::rollmeanr(daily_south2, k = 7, fill = NA)
  )

cat("  Daily total cases summary:\n")
print(summary(covid_daily$daily_total))
cat("  Max daily cases:", max(covid_daily$daily_total), "\n")
cat("  Mean daily cases:", mean(covid_daily$daily_total), "\n")

# ========================== 4. LOAD ERWAT DATA ==========================
cat("\n=== 4. Loading ERWAT Data ===\n")
erwat <- read.csv(file.path(data_dir, "erwat_water_care_works.csv"),
                  stringsAsFactors = FALSE)

cat("  Rows:", nrow(erwat), "  Unique WWTWs:", length(unique(erwat$Water_Care_Works)), "\n")

wwtp_lookup <- erwat %>%
  group_by(Water_Care_Works) %>%
  summarise(
    design_capacity_ML_day = mean(Design_Hydraulic_Capacity_ML_per_day, na.rm = TRUE),
    domestic_pct = mean(Domestic_Influent_Percent, na.rm = TRUE),
    industrial_pct = mean(Industrial_Influent_Percent, na.rm = TRUE),
    n_areas_served = n(),
    total_pop_served = sum(Estimated_Population_2011_Census, na.rm = TRUE),
    total_catchment_km2 = sum(Catchment_Area_km2, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n  WWTW Lookup Table:\n")
print(wwtp_lookup)

# ========================== 11. SAVE EVERYTHING ==========================
cat("\n=== 11. Saving Cleaned Data ===\n")

save(
  vulnerability, gene_data, gene_weekly, covid_daily, erwat,
  wwtp_lookup,
  output_dir,
  file = file.path(output_dir, "cleaned_data.RData")
)

cat("  Saved to:", file.path(output_dir, "cleaned_data.RData"), "\n")

write.csv(gene_data, file.path(output_dir, "gene_data_cleaned.csv"), row.names = FALSE)
write.csv(covid_daily, file.path(output_dir, "covid_daily_cleaned.csv"), row.names = FALSE)
write.csv(wwtp_lookup, file.path(output_dir, "wwtp_lookup.csv"), row.names = FALSE)

cat("\n=== SCRIPT 1 COMPLETE ===\n")
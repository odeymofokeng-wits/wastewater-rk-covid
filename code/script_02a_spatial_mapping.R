##############################################################################
# SCRIPT 2A: SPATIAL MAPPING (Part 1 of 2)
#
# This script:
#   Step 1: Map subplaces to WWTPs 
#   Step 2: Map wards to health sub-districts
#   Step 3: Assign health sub-districts to subplaces
#
# SAVES: script2a_mappings.RData (loaded by Script 2B)
#
# Prerequisites: Script 1 must have been run (cleaned_data.RData exists)
##############################################################################

# ========================== LIBRARIES ==========================
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(readxl)

# ========================== PATHS ==========================
path_subplaces    <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/WASTEWATER DATA/Ekurhuleni sub and main places/Eku_sub_places13.shp"
path_wwtp         <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/WASTEWATER DATA/Ekurhuleni WWTP/Ekurhuleni_WWTP.shp"
path_wwtp_service <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/WASTEWATER DATA/WWTP Serviced Areas-20260328T133206Z-1-001/WWTP Serviced Areas/Ekurhuleni_WWTP_Service_Areas.shp"
path_ward_sp      <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/WASTEWATER DATA/Wards of Served Sub Places/Ekurhuleni_Ward_Sub_Places.shp"
data_dir   <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/WASTEWATER DATA/Cleaned datasets"
output_dir <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/EDA_Output"

# ========================== LOAD DATA ==========================
cat("=== Loading Data ===\n")
load(file.path(output_dir, "cleaned_data.RData"))
sf_subplaces    <- st_read(path_subplaces, quiet = TRUE)
sf_wwtp         <- st_read(path_wwtp, quiet = TRUE)
sf_wwtp_service <- st_read(path_wwtp_service, quiet = TRUE)
sf_ward_sp      <- st_read(path_ward_sp, quiet = TRUE)
cat("  Subplaces:", nrow(sf_subplaces), "| WWTPs:", nrow(sf_wwtp),
    "| Service areas:", nrow(sf_wwtp_service), "| Ward-SP:", nrow(sf_ward_sp), "\n")

# ========================== STEP 1: WWTP-to-Subplace Mapping ==========================
cat("\n=== Step 1: WWTP-to-Subplace Mapping ===\n")

# 1A: Address-based mapping (7 of 8 WWTPs have addresses; JP Marais is NA)
wwtp_address_map <- data.frame(
  Address = c("Vorsterkroon, Nigel", "Plot 4, Holfontein Road, Etwawa Ext 1, Benoni",
              "Heidelberg Road, Plot 14, Maraisdrift, Nigel",
              "Corner Escombe and Wanderers Street, Brakpan",
              "Keramiek Road, Olifantsfontein",
              "69 Sarel Cilliers Street, Rynfield, Benoni",
              "Corner Brickfield/Bierman Streets, Vosloorus"),
  WWTP_Name = c("Carl Grundlingh", "Daveyton", "Herbert Bickley",
                "Jan Smuts", "Olifantsfontein", "Rynfield", "Vlakplaats"),
  stringsAsFactors = FALSE
)

sf_wwtp_service <- sf_wwtp_service %>%
  left_join(wwtp_address_map, by = "Address")

cat("  Mapped by Address:", sum(!is.na(sf_wwtp_service$WWTP_Name)),
    "of", nrow(sf_wwtp_service), "\n")
cat("  Unmapped (NA Address):", sum(is.na(sf_wwtp_service$WWTP_Name)), "\n")

# 1B: Use Drainage column from Ward-Subplace shapefile for unmapped subplaces
cat("\n  Drainage values in Ward-Subplace shapefile:\n")
print(table(sf_ward_sp$Drainage, useNA = "ifany"))

drainage_wwtp_map <- data.frame(
  Drainage = c("Carl Grundlingh & Herbert Bickley", "Daveyton", "Jan Smuts",
               "JP Marais", "Olifantsfontein", "Rynfield", "Vlakplaats"),
  WWTP_Drainage = c("Carl Grundlingh", "Daveyton", "Jan Smuts",
                    "JP Marais", "Olifantsfontein", "Rynfield", "Vlakplaats"),
  stringsAsFactors = FALSE
)

# Get Drainage per subplace (largest area if multiple wards)
sp_drainage <- sf_ward_sp %>%
  as.data.frame() %>%
  filter(!is.na(Drainage)) %>%
  group_by(SP_CODE) %>%
  arrange(desc(Shape_Area)) %>%
  summarise(Drainage = first(Drainage), n_ward_sp = n(), .groups = "drop")

cat("  Subplaces with Drainage:", nrow(sp_drainage), "of", nrow(sf_wwtp_service), "\n")

sp_drainage <- sp_drainage %>%
  left_join(drainage_wwtp_map, by = "Drainage")

cat("  Drainage mapped:", sum(!is.na(sp_drainage$WWTP_Drainage)),
    "| Unmapped:", sum(is.na(sp_drainage$WWTP_Drainage)), "\n")

# Resolve combined Drainage values (e.g. "Carl Grundlingh & Herbert Bickley")
sf_wwtp_proj <- st_transform(sf_wwtp, crs = 32736)
wwtp_coords <- st_coordinates(sf_wwtp_proj)
wwtp_names_in_shp <- sf_wwtp$Water_Care

sf_sub_proj <- st_transform(sf_subplaces, crs = 32736)
sp_centroids_all <- st_centroid(sf_sub_proj)

combined_mask <- grepl("&", sp_drainage$Drainage) & !is.na(sp_drainage$Drainage)
if (any(combined_mask)) {
  cat("  Resolving combined Drainage values...\n")
  for (sp_code in sp_drainage$SP_CODE[combined_mask]) {
    row_idx <- which(sp_drainage$SP_CODE == sp_code)
    parts <- trimws(unlist(strsplit(sp_drainage$Drainage[row_idx], "&")))
    match_idx <- sapply(parts, function(p) {
      w <- which(grepl(p, wwtp_names_in_shp, ignore.case = TRUE))
      if (length(w) > 0) w[1] else NA
    })
    match_idx <- na.omit(match_idx)
    if (length(match_idx) > 0) {
      if (length(match_idx) == 1) {
        sp_drainage$WWTP_Drainage[row_idx] <- wwtp_names_in_shp[match_idx[1]]
      } else {
        sp_cx <- st_coordinates(sp_centroids_all[sp_centroids_all$SP_CODE == sp_code, ])
        if (length(sp_cx) > 0) {
          d <- sqrt(rowSums((sweep(wwtp_coords, 2, sp_cx))^2))
          sp_drainage$WWTP_Drainage[row_idx] <- wwtp_names_in_shp[match_idx[which.min(d[match_idx])]]
        }
      }
    }
  }
}

cat("  After resolving combined Drainage:\n")
print(table(sp_drainage$WWTP_Drainage, useNA = "ifany"))

# 1C: Merge Drainage mapping into service areas
sf_wwtp_service <- sf_wwtp_service %>%
  left_join(sp_drainage %>% select(SP_CODE, WWTP_Drainage), by = "SP_CODE")

# 1D: Final assignment logic
sf_wwtp_service <- sf_wwtp_service %>%
  mutate(WWTP_Final = case_when(
    !is.na(WWTP_Name) ~ WWTP_Name,
    !is.na(WWTP_Drainage) ~ WWTP_Drainage,
    TRUE ~ NA_character_
  ))

# Last resort: nearest WWTP by centroid
n_unmapped <- sum(is.na(sf_wwtp_service$WWTP_Final))
if (n_unmapped > 0) {
  cat("  Still unmapped:", n_unmapped, "-> assigning to nearest WWTP\n")
  
  # Transform sf_subplaces to same CRS as sf_wwtp_proj
  sf_sub_proj <- st_transform(sf_subplaces, crs = st_crs(sf_wwtp_proj))
  
  for (sp_code in sf_wwtp_service$SP_CODE[is.na(sf_wwtp_service$WWTP_Final)]) {
    d <- st_distance(sf_sub_proj[sf_sub_proj$SP_CODE == sp_code, ], sf_wwtp_proj)
    sf_wwtp_service$WWTP_Final[sf_wwtp_service$SP_CODE == sp_code] <-
      wwtp_names_in_shp[which.min(as.numeric(d))]
  }
}

sf_wwtp_service <- sf_wwtp_service %>%
  mutate(WWTP_Name = WWTP_Final) %>%
  select(-WWTP_Final, -WWTP_Drainage)

cat("\n  *** FINAL WWTP ASSIGNMENT ***\n")
print(table(sf_wwtp_service$WWTP_Name, useNA = "ifany"))

jp_count <- sum(sf_wwtp_service$WWTP_Name == "JP Marais", na.rm = TRUE)
cat("\n  SANITY CHECK: JP Marais =", jp_count,
    ifelse(jp_count > 100, " *** WARNING: TOO MANY ***", " *** OK (expected ~29) ***"), "\n")

wwtp_sp_mapping <- sf_wwtp_service %>%
  st_drop_geometry() %>%
  select(SP_CODE, SP_NAME, WWTP_Name, WWTP_CODE, Address)
write.csv(wwtp_sp_mapping, file.path(output_dir, "wwtp_subplace_mapping.csv"), row.names = FALSE)

# ========================== STEP 2: Health Sub-District Mapping ==========================
cat("\n=== Step 2: Health Sub-District Mapping ===\n")

hospitals_enriched <- read_excel(file.path(data_dir, "Merged_Hospital_Dataset.xlsx"))
cat("  Hospitals:", nrow(hospitals_enriched), "| With SD:",
    sum(!is.na(hospitals_enriched$Subdistrict)), "\n")

hospital_sd <- hospitals_enriched %>%
  filter(!is.na(Subdistrict)) %>%
  # Filter to only Ekurhuleni subdistricts
  filter(grepl("Ekurhuleni", Subdistrict, ignore.case = TRUE)) %>%
  mutate(
    sd_code = trimws(gsub(" Health Sub-District", "", Subdistrict, ignore.case = TRUE)),
    sd_code = trimws(gsub("Ekurhuleni ", "", sd_code, ignore.case = TRUE)),
    # Handle both full names and abbreviations
    covid_sd = case_when(
      grepl("[Nn]orth.*1|^[Nn]1$", sd_code) ~ "north1",
      grepl("[Nn]orth.*2|^[Nn]2$", sd_code) ~ "north2",
      grepl("[Ee]ast.*1|^[Ee]1$", sd_code) ~ "east1",
      grepl("[Ee]ast.*2|^[Ee]2$", sd_code) ~ "east2",
      grepl("[Ss]outh.*1|^[Ss]1$", sd_code) ~ "south1",
      grepl("[Ss]outh.*2|^[Ss]2$", sd_code) ~ "south2",
      TRUE ~ NA_character_
    ),
    ward_no = as.integer(gsub("^ZA7970*", "", Ward_Code)),
    ward_no = ifelse(grepl("^ZA797\\d{4}$", Ward_Code), ward_no, NA_integer_)
  ) %>%
  select(Hospital_Name, sd_code, covid_sd, ward_no, Longitude, Latitude)

ward_sd_lookup <- hospital_sd %>%
  filter(!is.na(ward_no), !is.na(covid_sd)) %>%
  distinct(ward_no, .keep_all = TRUE) %>%
  select(ward_no, covid_sd)

cat("  Directly mapped wards:", nrow(ward_sd_lookup), "\n")

# Nearest-hospital fallback for unmapped wards
hospitals_sf_all <- st_as_sf(hospitals_enriched, coords = c("Longitude", "Latitude"),
                              crs = 4326, remove = FALSE) %>%
  filter(!is.na(Subdistrict)) %>%
  left_join(hospital_sd %>% select(Hospital_Name, covid_sd), by = "Hospital_Name")

sf_ward_sp <- sf_ward_sp %>%
  mutate(ward_no = as.integer(gsub("[^0-9]", "", as.character(WardID))))

ward_centroids_sf <- sf_ward_sp %>%
  filter(!is.na(ward_no)) %>%
  group_by(ward_no) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_transform(crs = 4326)

hospitals_sf_proj <- st_transform(hospitals_sf_all, crs = 32736)
ward_centroids_proj <- st_transform(ward_centroids_sf, crs = 32736)

nearest_hosp <- st_nearest_feature(ward_centroids_proj, hospitals_sf_proj)
ward_centroids_sf$covid_sd_nearest <- hospitals_sf_all$covid_sd[nearest_hosp]

ward_sd_full <- ward_centroids_sf %>%
  left_join(ward_sd_lookup %>% rename(covid_sd_direct = covid_sd), by = "ward_no") %>%
  mutate(
    covid_sd = ifelse(!is.na(covid_sd_direct), covid_sd_direct, covid_sd_nearest),
    source = ifelse(!is.na(covid_sd_direct), "hospital_direct", "nearest_neighbor")
  ) %>%
  st_drop_geometry() %>%
  select(ward_no, covid_sd, source)

cat("  Ward sub-district mapping:\n")
print(table(ward_sd_full$covid_sd, useNA = "ifany"))
cat("  Sources:\n")
print(table(ward_sd_full$source, useNA = "ifany"))

write.csv(ward_sd_full, file.path(output_dir, "ward_subdistrict_full_mapping.csv"), row.names = FALSE)

# ========================== STEP 3: Assign Sub-districts to Subplaces ==========================
cat("\n=== Step 3: Assign Health Sub-districts to Subplaces ===\n")

sf_ward_sp_sd <- sf_ward_sp %>%
  left_join(ward_sd_full %>% select(ward_no, covid_sd), by = "ward_no")

sp_sd_assignment <- sf_ward_sp_sd %>%
  as.data.frame() %>%
  group_by(SP_CODE) %>%
  arrange(desc(Shape_Area)) %>%
  summarise(
    SP_NAME = first(SP_NAME),
    health_subdistrict = covid_sd[1],
    n_wards = n(),
    n_unique_sd = n_distinct(covid_sd, na.rm = TRUE),
    .groups = "drop"
  )

cat("  Mapped:", sum(!is.na(sp_sd_assignment$health_subdistrict)), "of",
    nrow(sp_sd_assignment), "\n")
print(table(sp_sd_assignment$health_subdistrict, useNA = "ifany"))

write.csv(sp_sd_assignment, file.path(output_dir, "sp_subdistrict_assignment.csv"), row.names = FALSE)

# ========================== SAVE INTERMEDIATE ==========================
cat("\n=== Saving intermediate results ===\n")

save(sf_subplaces, sf_wwtp, sf_wwtp_service, sf_ward_sp,
     wwtp_sp_mapping, ward_sd_full, sp_sd_assignment,
     file = file.path(output_dir, "script2a_mappings.RData"))

cat("  Saved: script2a_mappings.RData\n")
cat("\n=== SCRIPT 2A COMPLETE ===\n")

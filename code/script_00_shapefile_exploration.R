##############################################################################
# SCRIPT 0: SHAPEFILE EXPLORATION



# Shapefile paths 
path_subplaces   <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/WASTEWATER DATA/Ekurhuleni sub and main places/Eku_sub_places13.shp"

path_main_places <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/WASTEWATER DATA/Ekurhuleni sub and main places/Eku_main_places13.shp"

path_wwtp <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/WASTEWATER DATA/Ekurhuleni WWTP/Ekurhuleni_WWTP.shp"

path_treatment   <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/WASTEWATER DATA/Treatment Plants/Treatment Plants.shp"

path_ward_sp     <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/WASTEWATER DATA/Wards of Served Sub Places/Ekurhuleni_Ward_Sub_Places.shp"

path_wwtp_service <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/WASTEWATER DATA/WWTP Serviced Areas-20260328T133206Z-1-001/WWTP Serviced Areas/Ekurhuleni_WWTP_Service_Areas.shp"

# Output directory for reports
output_dir <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/EDA_Output"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ========================== LIBRARIES ==========================
library(sf)
library(dplyr)
library(readr)
library(knitr)

# ========================== HELPER ==========================
report_shapefile <- function(path, name) {
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("SHAPEFILE:", name, "\n")
  cat("Path:", path, "\n")
  cat(strrep("=", 70), "\n")
  
  sf_obj <- tryCatch(
    st_read(path, quiet = TRUE),
    error = function(e) {
      cat("ERROR reading file:", e$message, "\n")
      return(NULL)
    }
  )
  
  if (is.null(sf_obj)) return(NULL)
  
  cat("Geometry type:", st_geometry_type(sf_obj)[1], "\n")
  cat("CRS:", st_crs(sf_obj)$proj4string, "\n")
  cat("Dimensions:", nrow(sf_obj), "features x", ncol(sf_obj), "columns\n")
  cat("\nColumn names and types:\n")
  cat(paste(names(sf_obj), sapply(sf_obj, class), sep = ": "), sep = "\n")
  
  cat("\nFirst 5 rows of attribute table:\n")
  print(head(as.data.frame(sf_obj), 5))
  
  cat("\nUnique values in key columns (max 20):\n")
  for (col in names(sf_obj)) {
    if (is.character(sf_obj[[col]]) || is.factor(sf_obj[[col]])) {
      n_unique <- length(unique(sf_obj[[col]]))
      cat(sprintf("  %s: %d unique values", col, n_unique), "\n")
      if (n_unique <= 30) {
        cat("    Values:", paste(head(sort(unique(sf_obj[[col]])), 20), collapse = ", "), "\n")
      }
    }
  }
  
  # Check for SP_CODE 
  sp_cols <- grep("SP|sp_code|SP_CODE|SP_NAME|SUBPLACE|subplace", 
                  names(sf_obj), value = TRUE, ignore.case = TRUE)
  if (length(sp_cols) > 0) {
    cat("\nPotential SP_CODE columns found:", paste(sp_cols, collapse = ", "), "\n")
    for (col in sp_cols) {
      cat(sprintf("  %s: %d unique, sample:", col, length(unique(sf_obj[[col]]))),
          head(sort(unique(sf_obj[[col]])), 5), "\n")
    }
  }
  
  # Check for WWTP 
  wwtp_cols <- grep("WWTP|wwtp|Water|WATER|Plant|PLANT|TREAT|TREAT", 
                    names(sf_obj), value = TRUE, ignore.case = TRUE)
  if (length(wwtp_cols) > 0) {
    cat("\nPotential WWTP columns found:", paste(wwtp_cols, collapse = ", "), "\n")
    for (col in wwtp_cols) {
      cat(sprintf("  %s: %d unique, sample:", col, length(unique(sf_obj[[col]]))),
          head(sort(unique(sf_obj[[col]])), 10), "\n")
    }
  }
  
  # Check for ward columns
  ward_cols <- grep("WARD|ward|Ward", names(sf_obj), value = TRUE, ignore.case = TRUE)
  if (length(ward_cols) > 0) {
    cat("\nPotential Ward columns found:", paste(ward_cols, collapse = ", "), "\n")
    for (col in ward_cols) {
      cat(sprintf("  %s: %d unique", col, length(unique(sf_obj[[col]]))), "\n")
    }
  }
  
  return(sf_obj)
}

# ========================== RUN EXPLORATION ==========================
cat("SHAPEFILE EXPLORATION REPORT")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# 1. Subplaces
sf_subplaces <- report_shapefile(path_subplaces, "Subplaces")

# 2. Main Places  
sf_main <- report_shapefile(path_main_places, "Main Places")

# 3. WWTP locations
sf_wwtp <- report_shapefile(path_wwtp, "WWTP Locations")

# 4. Treatment Plants
sf_treatment <- report_shapefile(path_treatment, "Treatment Plants")

# 5. Ward-Subplace intersections
sf_ward_sp <- report_shapefile(path_ward_sp, "Ward-Subplace")

# 6. WWTP Service Areas
sf_wwtp_service <- report_shapefile(path_wwtp_service, "WWTP Service Areas")

# ========================== CROSS-CHECK ==========================
cat("\n\n")
cat(strrep("=", 70), "\n")
cat("CROSS-CHECKING IDENTIFIERS\n")
cat(strrep("=", 70), "\n")

# Check SP_CODE match between vulnerability CSV and shapefile
vulnerability <- read.csv("C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/WASTEWATER DATA/Cleaned datasets/Ekurhuleni_Vulnerability_Subplaces_Cleaned.csv")
cat("\nVulnerability file has", nrow(vulnerability), "subplaces")
cat("Subplace_Code range:", min(vulnerability$Subplace_Code), "to", max(vulnerability$Subplace_Code))
cat("Shapefile has", nrow(sf_subplaces), "features")

# If shapefile loaded, try to match
if (!is.null(sf_subplaces)) {
  # Try various possible column names for SP code
  possible_cols <- c("SP_CODE", "SP_CODE_", "Sp_Code", "SPCODE", 
                     "SUBPLACE_C", "CODE", "SP_COD", "OBJECTID", "FID",
                     names(sf_subplaces)[grepl("CODE|code|SP_|SP", names(sf_subplaces))])
  
  for (col in possible_cols) {
    if (col %in% names(sf_subplaces)) {
      cat(sprintf("\nTrying column '%s':", col))
      cat("  Unique values:", length(unique(sf_subplaces[[col]])))
      cat("  Sample:", head(sort(unique(sf_subplaces[[col]])), 5))
      
      # Try matching with vulnerability Subplace_Code
      if (is.numeric(sf_subplaces[[col]])) {
        match_count <- sum(sf_subplaces[[col]] %in% vulnerability$Subplace_Code)
        cat("  Match with vulnerability Subplace_Code:", match_count)
      } else if (is.character(sf_subplaces[[col]])) {
        match_count <- sum(as.character(sf_subplaces[[col]]) %in% as.character(vulnerability$Subplace_Code))
        cat("  Match with vulnerability Subplace_Code:", match_count)
      }
    }
  }
}

colnames(vulnerability)

# Check WWTP name match between ERWAT CSV and shapefile
erwat <- read.csv("C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/WASTEWATER DATA/Cleaned datasets/erwat_water_care_works.csv")
cat("\n\nERWAT has", length(unique(erwat$Water_Care_Works)), "unique WWTWs:")
cat(paste(unique(erwat$Water_Care_Works), collapse = "\n  "), "\n")

if (!is.null(sf_wwtp)) {
  cat("\nWWTP shapefile has", nrow(sf_wwtp), "features")
  # Print all column names and all values for text columns
  for (col in names(sf_wwtp)) {
    if (is.character(sf_wwtp[[col]]) || is.factor(sf_wwtp[[col]])) {
      cat(sprintf("\n  Column '%s' values:", col))
      print(unique(sf_wwtp[[col]]))
    }
  }
}

if (!is.null(sf_wwtp_service)) {
  cat("\nWWTP Service Areas shapefile has", nrow(sf_wwtp_service), "features")
  for (col in names(sf_wwtp_service)) {
    if (is.character(sf_wwtp_service[[col]]) || is.factor(sf_wwtp_service[[col]])) {
      cat(sprintf("\n  Column '%s': %d unique values", col, length(unique(sf_wwtp_service[[col]]))))
      if (length(unique(sf_wwtp_service[[col]])) <= 30) {
        print(unique(sf_wwtp_service[[col]]))
      }
    }
  }
}

cat("\n\nDONE.\n")

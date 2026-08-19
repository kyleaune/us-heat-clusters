#---- Script Metadata #----
# Title: Data Wrangling
# Author: Kyle T. Aune, PhD, MPH
# Date: 03/03/2026
#--------------------------

# Creating spatial weighted average values for CONUS counties of:
# 1)  Green space (% coverage) -
# 2)  Tree canopy (% coverage) - dl, process, save
# 3)  Water bodies (% coverage) - dl, process, save
# 4)  Land use - impervious
# 5)  Land use - manmade impervious
# 6)  Albedo
# 7)  Albedo - impervious
# 8)  Albedo - manmade imperviousm
# 9)  Air conditioning ownership - dl
# 10) Residential energy cost : HH income
# 11) Combine all measures


#---- Setup #----

pkgs <- c("tidyverse", "sf", "stars", "raster", "tidycensus", "tigris", "exactextractr",
          "parallel", "doParallel", "foreach", "readxl", "httr", "jsonlite")

install.packages(setdiff(pkgs, rownames(installed.packages())))

lapply(pkgs, library, character.only = TRUE, quietly = TRUE, verbose = FALSE)

setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Research/260204 - US Heat Clusters/Data")

# Downloading US county shapefile
co <- counties(state = state.abb, cb = TRUE) %>%
  # Filtering for contiguous states (to match extent of heat data)
  filter(!STUSPS %in% c("AK", "HI"))


#---- Green Space #----

# Workflow
# 1) Read in monthly rasters
# 2) Mask by quality assurance layer
# 3) Calculate median NDVI value for each cell (7/1/2020-8/31/2020)
# 4) Reclassify to 0, 1 based on green space vs. not according to
#   - 0.2 (less conservative, https://www.mdpi.com/2073-445X/11/3/351#B1-land-11-00351)
#   - 0.3 (more conservative, unclear ref but apparently common)
# 5) Calculate % areal coverage by green space of each CONUS county


#---- Tree Canopy Coverage #----

# Workflow
# 1) Read in 2020 USFS raster
# 2) Calculate % areal coverage by tree canopy coverage of each CONUS county

# Read in TCC data
tcc <- raster("TCC/science_tcc_conus_wgs84_v2023-5_20200101_20201231.tif")

# Calculate mean % areal coverage
tcc.comean <- co %>%
  mutate(tcc = exact_extract(tcc, co %>% st_transform(st_crs(tcc)), "mean"))

# Saving county TCC shapefile
st_write(tcc.comean, "TCC/tcc_county_conus_2020.shp")


#---- Water Bodies #----

# Calculating water body percentages from census TIGER/Lines
co.water <- co %>%
  mutate(water.pct = AWATER / (ALAND + AWATER))

# Saving county water area shapefile
st_write(co.water, "Water/water_area_county_conus_2024.shp")


#---- Land Use #----

# Read in annual land use raster (downloaded from https://www.mrlc.gov/data)


#---- Residential Energy Cost #----

# Read in Form 861 (2020) - annual sales ($ and mWh) by utility
eng.20 <- read_xlsx("EIA Form 861 - 2020/Sales_Ult_cust_2020.xlsx",
                 sheet = "States",
                 skip = 3,
                 n_max = 2656,
                 col_names = c("year", "utility_no", "utility", "part", "service",
                               "data_type", "state", "own", "ba",
                               paste(rep(c("rev", "sales", "cust"), times = 5),
                                     rep(c("res", "com", "ind", "trans", "total"), each = 3), sep = "."))) %>%
  mutate(across(.cols = c(rev.res:cust.trans),
                .fns = ~ as.numeric(na_if(., "." ))))

# Read in Form 861 (2024) - annual sales ($ and mWh) by utility
eng.24 <- read_xlsx("EIA Form 861 - 2024/Sales_Ult_cust_2024.xlsx",
                    sheet = "States",
                    skip = 3,
                    n_max = 2656,
                    col_names = c("year", "utility_no", "utility", "part", "service",
                                  "data_type", "state", "own", "ba",
                                  paste(rep(c("rev", "sales", "cust"), times = 5),
                                        rep(c("res", "com", "ind", "trans", "total"), each = 3), sep = "."))) %>%
  mutate(across(.cols = c(rev.res:cust.trans),
                .fns = ~ as.numeric(na_if(., "." )))) %>%
  # Dropping non-Part A (A = bundled utilities (producer/supplier), B = energy production
  # only, C = energy delivery only, D = largely unused code)
  filter(part == "A")

# Cleaning and calculating 2020 energy data
eng.20 <- eng.20 %>%
  # Dropping Alaska and Hawaii
  filter(state %in% c("AK", "HI") == FALSE) %>%
    ## Dropped 15 utilities
  # Dropping non-Part A (A = bundled utilities (producer/supplier), B = energy production
  # only, C = energy delivery only, D = largely unused code)
  filter(part == "A") %>%
    ## 1,039 utilities dropped
  # Dropping utilities with no residential sales
  filter(rev.res > 0) %>%
  filter(is.na(rev.res) == FALSE) %>%
    ## 98 utilities dropped
  # Dropping 'Behind the Meter' companies
  filter(own != "Behind the Meter") %>%
    ## 192 utilities dropped, 1,312 utilities in 2020
  # Calculate cost per kWh
  mutate(cost.res = rev.res / sales.res)

# 3) Read in utility territory shapefile
util.url <- "https://services5.arcgis.com/HDRa0B57OVrv2E1q/ArcGIS/rest/services/Electric_Retail_Service_Territories/FeatureServer/0/query"

# Function to download most up to date shapefile of utility boundaries
dlfx <- function(url, batch_size = 2000) {

  all_chunks <- list()
  offset <- 0

  repeat {

    query <- list(
      where = "1=1",
      outFields = "*",
      outSR = "4326",
      f = "json",
      resultOffset = offset,
      resultRecordCount = batch_size
    )

    res <- GET(url, query = query)
    stop_for_status(res)

    txt <- content(res, as = "text", encoding = "UTF-8")
    json <- fromJSON(txt, simplifyVector = FALSE)

    # Stop if no features present
    if (length(json$features) == 0) break

    # Convert json to sf
    tmp_file <- tempfile(fileext = ".json")
    writeLines(txt, tmp_file)

    sf_chunk <- st_read(tmp_file, quiet = TRUE)

    all_chunks[[length(all_chunks) + 1]] <- sf_chunk

    # Download progress message
    message(paste("Downloaded", offset + nrow(sf_chunk), "features"))

    if (nrow(sf_chunk) < batch_size) break

    offset <- offset + batch_size
  }

  bind_rows(all_chunks)
}

util.sf <- dlfx(util.url)

st_write(util.sf, ("Electricity Utility Providers/ERST.shp"))


territories_sf <- st_read(
  "https://services5.arcgis.com/HDRa0B57OVrv2E1q/ArcGIS/rest/services/Electric_Retail_Service_Territories/FeatureServer/0",
  quiet = FALSE
)






temp <- tempfile(fileext = ".geojson")
httr::GET(url, httr::write_disk(temp, overwrite = TRUE))
st_read(temp)
st_read("/Users/kta5166/Downloads/erst.json")

util.sf <- st_read("https://services3.arcgis.com/OYP7N6mAJJCyH6hd/ArcGIS/rest/services/Electric_Retail_Service_Territories_HIFLD/FeatureServer/0?f=json")

util.sf <- st_read("Electricity Utility Providers/Electric-Retail-Service-Territories.shp")

table(eng.20$utility_no %in% util.sf$ID)
table(str_to_upper(eng.20$utility) %in% util.sf$NAME)

eng.20[eng.20$utility_no %in% util.sf$ID == FALSE, ]
  ## Gulf Power in Florida is missing
eng.20[str_to_upper(eng.20$utility) %in% util.sf$NAME == FALSE, ]$utility


url <- "https://services3.arcgis.com/OYP7N6mAJJCyH6hd/ArcGIS/rest/services/Electric_Retail_Service_Territories_HIFLD/FeatureServer/0?f=json"

  jsonlite::validate(readLines(url, n = 1))

# 4) Spatial join to counties to determine customer-weighted spatial average cost












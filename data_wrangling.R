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
          "parallel", "doParallel", "foreach")

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

# Workflow:
# 1) Read in Form 861 (2020) - annual sales ($ and mWh) by utility
# 2) Calculate cost per kWh (Sales_Ult_Cust_2020.xslx)
# 3) Read in utility territory shapefile
# 4) Spatial join to counties to determine customer-weighted spatial average cost
#---- Script Metadata #----
# Title: heat cluster
# Author: Hang Li
# Date: 07/08/2026


#---- Setup #----

pkgs <- c("tidyverse", "sf", "stars", "raster", "tidycensus", "tigris", "exactextractr",
          "parallel", "doParallel", "foreach", "readxl", "httr", "jsonlite")
lapply(pkgs, library, character.only = TRUE, quietly = TRUE, verbose = FALSE)

setwd("C:/Users/velvet/The Pennsylvania State University/Aune, Kyle T - 260204 - US Heat Clusters/Data/")

#---- Heat clustering#----
heat <- readRDS("Heatvars_County_2000-2020_v1.2.Rds")

heat2 <- heat %>%
  mutate(
    Date = as.Date(Date),
    doy = yday(Date),
    month = month(Date)
  )
wbgt_clim <- heat %>%
  mutate(
    Date = as.Date(Date),
    month = month(Date),
    day = day(Date)
  ) %>%
  filter(!(month == 2 & day == 29)) %>%      # remove leap day
  group_by(StCoFIPS, month, day) %>%
  summarise(
    WBGT = mean(WBGTmax_C, na.rm = TRUE),
    .groups = "drop"
  )
head(wbgt_clim)
county_clusters <- wbgt_clim %>%
  group_by(StCoFIPS) %>%
  group_modify(~{
    
    df <- .x
    
    # remove NA
    df <- df %>%
      filter(!is.na(WBGT))
    
    # skip counties with too few distinct values
    if(n_distinct(df$WBGT) < 3){
      
      df$cluster <- NA_integer_
      df$season  <- NA_character_
      
      return(df)
      
    }
    
    # K-means
    km <- kmeans(df$WBGT,
                 centers = 3,
                 nstart = 100)
    
    df$cluster <- km$cluster
    
    # Order clusters by mean WBGT
    cluster_order <- df %>%
      group_by(cluster) %>%
      summarise(
        mean_wbgt = mean(WBGT),
        .groups = "drop"
      ) %>%
      arrange(mean_wbgt) %>%
      mutate(
        season = c("Cold", "Temperate", "Warm")
      )
    
    df %>%
      left_join(cluster_order, by = "cluster")
    
  }) %>%
  ungroup()
month_cluster_pct <- county_clusters %>%
  group_by(StCoFIPS, month, season) %>%
  summarise(
    n_days = n(),
    .groups = "drop"
  ) %>%
  group_by(StCoFIPS, month) %>%
  mutate(
    total_days = sum(n_days),
    pct_days = n_days / total_days
  ) %>%
  ungroup()

head(month_cluster_pct)
monthly_season <- month_cluster_pct %>%
  group_by(StCoFIPS, month) %>%
  slice_max(
    order_by = pct_days,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(StCoFIPS, month)

head(county_clusters)

write_csv(monthly_season, "C:/Users/velvet/OneDrive - The Pennsylvania State University/Paper/US_heat_cluster/us-heat-clusters/data/processed/season cluster/monthly_season.csv")


#---- County heatwaves#----
head(heat2)

# 95th percentile of WBGT for each county during warm season
heat_warm_95th <- heat2 %>%
  left_join(
    monthly_season %>%
      dplyr::select(StCoFIPS, month, season),
    by = c("StCoFIPS", "month")
  ) %>%
  filter(season == "Warm") %>%
  group_by(StCoFIPS) %>%
  summarise(
    WBGTmin_95th  = quantile(WBGTmin_C, 0.95, na.rm = TRUE),
    WBGTmean_95th = quantile(WBGTmean_C, 0.95, na.rm = TRUE),
    WBGTmax_95th  = quantile(WBGTmax_C, 0.95, na.rm = TRUE),
    .groups = "drop"
  )


# dummy value for 3 time periods * 3 wbgt definition
## 1 day heatwave events
heat_hw <- heat2 %>%
  left_join(monthly_season %>% dplyr::select(StCoFIPS, month, season),
            by = c("StCoFIPS", "month")) %>%
  filter(season == "Warm") %>%
  left_join(heat_warm_95th, by = "StCoFIPS") %>%
  mutate(
    hw1_min  = WBGTmin_C  >= WBGTmin_95th,
    hw1_mean = WBGTmean_C >= WBGTmean_95th,
    hw1_max  = WBGTmax_C  >= WBGTmax_95th
  )

mark_heatwave <- function(x, min_length) {
  
  r <- rle(x)
  
  inverse.rle(list(
    values = r$values & r$lengths >= min_length,
    lengths = r$lengths
  ))
}
## 3- and 5- consecutive days heatwave events
heat_hw <- heat2 %>%
  left_join(
    monthly_season %>%
      dplyr::select(StCoFIPS, month, season),
    by = c("StCoFIPS", "month")
  ) %>%
  filter(season == "Warm") %>%
  left_join(heat_warm_95th, by = "StCoFIPS") %>%
  mutate(
    hw1_min  = WBGTmin_C  >= WBGTmin_95th,
    hw1_mean = WBGTmean_C >= WBGTmean_95th,
    hw1_max  = WBGTmax_C  >= WBGTmax_95th
  ) %>%
  arrange(StCoFIPS, Date) %>%
  group_by(StCoFIPS) %>%
  mutate(
    hw3_min  = mark_heatwave(hw1_min, 3),
    hw3_mean = mark_heatwave(hw1_mean, 3),
    hw3_max  = mark_heatwave(hw1_max, 3),
    
    hw5_min  = mark_heatwave(hw1_min, 5),
    hw5_mean = mark_heatwave(hw1_mean, 5),
    hw5_max  = mark_heatwave(hw1_max, 5)
  ) %>%
  ungroup() %>%
  dplyr::select(
    StCoFIPS, Date,
    WBGTmin_C, WBGTmean_C, WBGTmax_C,
    WBGTmin_95th, WBGTmean_95th, WBGTmax_95th,
    hw1_min, hw1_mean, hw1_max,
    hw3_min, hw3_mean, hw3_max,
    hw5_min, hw5_mean, hw5_max,
    doy, month, season
  )

write_csv(heat_hw, "C:/Users/velvet/OneDrive - The Pennsylvania State University/Paper/US_heat_cluster/us-heat-clusters/data/processed/season cluster/heat_wave.csv")
head(heat_hw)


#---- Heatwave characteristics#----

library(dplyr)
calc_heatwave_stats <- function(flag, value){
  
  r <- rle(flag)
  
  event_lengths <- r$lengths[r$values]
  
  frequency <- length(event_lengths)
  
  duration <- if(frequency == 0) NA else mean(event_lengths)
  
  event_id <- rep(seq_along(r$lengths), r$lengths)
  
  intensity <- tapply(value[flag], event_id[flag], mean)
  
  intensity <- if(frequency == 0) NA else mean(intensity)
  
  c(
    frequency = frequency,
    duration = duration,
    intensity = intensity
  )
}


heatwave_characteristics <-
  
  heat_hw %>%
  group_split(StCoFIPS) %>%
  lapply(function(df){
    
    data.frame(
      
      StCoFIPS = df$StCoFIPS[1],
      
      t(calc_heatwave_stats(df$hw1_min,  df$WBGTmin_C)),
      t(calc_heatwave_stats(df$hw1_mean, df$WBGTmean_C)),
      t(calc_heatwave_stats(df$hw1_max,  df$WBGTmax_C)),
      
      t(calc_heatwave_stats(df$hw3_min,  df$WBGTmin_C)),
      t(calc_heatwave_stats(df$hw3_mean, df$WBGTmean_C)),
      t(calc_heatwave_stats(df$hw3_max,  df$WBGTmax_C)),
      
      t(calc_heatwave_stats(df$hw5_min,  df$WBGTmin_C)),
      t(calc_heatwave_stats(df$hw5_mean, df$WBGTmean_C)),
      t(calc_heatwave_stats(df$hw5_max,  df$WBGTmax_C))
      
    )
    
  }) %>%
  bind_rows()

names(heatwave_characteristics) <-
  c(
    "StCoFIPS",
    
    "hw1_min_frequency","hw1_min_duration","hw1_min_intensity",
    "hw1_mean_frequency","hw1_mean_duration","hw1_mean_intensity",
    "hw1_max_frequency","hw1_max_duration","hw1_max_intensity",
    
    "hw3_min_frequency","hw3_min_duration","hw3_min_intensity",
    "hw3_mean_frequency","hw3_mean_duration","hw3_mean_intensity",
    "hw3_max_frequency","hw3_max_duration","hw3_max_intensity",
    
    "hw5_min_frequency","hw5_min_duration","hw5_min_intensity",
    "hw5_mean_frequency","hw5_mean_duration","hw5_mean_intensity",
    "hw5_max_frequency","hw5_max_duration","hw5_max_intensity"
  )

head(heatwave_characteristics)
write_csv(heat_hw, "C:/Users/velvet/OneDrive - The Pennsylvania State University/Paper/US_heat_cluster/us-heat-clusters/data/processed/season cluster/heat_wave_characteristics.csv")















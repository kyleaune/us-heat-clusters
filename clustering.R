#---- Script Metadata #----
# Title: clustering
# Author: Hang Li
# Date: 07/23/2026
#---- 1. Set up #----
rm(list = ls())
pkgs <- c(  "dplyr",
            "tidyr",
            "readr",
            "tibble",
            "sf",
            "tigris",
            "psych",
            "corrplot",
            "ggplot2",
            "cluster",
            "RColorBrewer",
            "NbClust",
            "tmap")
invisible(lapply(pkgs, library, character.only = TRUE))

setwd("/Users/hangli/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Aune, Kyle T's files - 260204 - US Heat Clusters/Data/")
tract_master <- readRDS("tract_master_2020.rds")
heat_char <- read_csv("processed/season cluster/heat_wave_characteristics.csv")
#---- 2. Heat Vulnerability #----

pca_data <- tract_master %>%
  st_drop_geometry() %>%
  filter(valid_ct) %>%
  select(
    GEOID,
    water_pct,
    pop_density,
    pct_age65,
    pct_age5,
    median_income,
    poverty,
    energy_burden,
    electricity_burden,
    AC_central,
    AC_none,
    tree_canopy,
    albedo,
    ndvi,
    dev_open_pct,
    dev_low_pct,
    dev_medium_pct,
    dev_high_pct,
    building_coverage,
    impervious
  )
vars <- pca_data %>%
  mutate(
    water_pct     = -water_pct,
    median_income = -median_income,
    AC_central    = -AC_central,
    tree_canopy   = -tree_canopy,
    albedo        = -albedo,
    ndvi          = -ndvi
  ) %>%
  na.omit()

vars_scaled <- scale(vars %>% select(-GEOID))
pc <- principal(
  vars_scaled,
  nfactors = 7,
  rotate = "varimax",
  scores = TRUE
)

pca_scores <- vars %>%
  select(GEOID) %>%
  bind_cols(
    as.data.frame(pc$scores)
  )

names(pca_scores)[2:8] <- c(
  "PC1_urban",
  "PC2_socioeconomic",
  "PC3_age_structure",
  "PC4_AC",
  "PC5_greenness",
  "PC6_low_density",
  "PC7_water"
)
tract_master <- tract_master %>%
  left_join(pca_scores, by = "GEOID")

tract_pca <- tract_master %>%
  st_drop_geometry() %>%
  select(
    GEOID,
    PC1_urban,
    PC2_socioeconomic,
    PC3_age_structure,
    PC4_AC,
    PC5_greenness,
    PC6_low_density,
    PC7_water
  )
# write.csv(tract_pca, "processed/season cluster/tract_pca.csv")

#---- 3. Comparison #----
tract_pca<- read.csv("processed/season cluster/tract_pca.csv")
county_cluster <- read.csv("processed/season cluster/county_clustered.csv")

tract_pca <- tract_pca %>%
  mutate(
    GEOID = sprintf("%011.0f", GEOID),
    StCoFIPS = as.integer(substr(GEOID, 1, 5))
  )

tract_pca <- tract_pca %>%
  left_join(
    county_cluster,
    by = "StCoFIPS"
  )
names(tract_pca)

library(dplyr)

pc_summary <- tract_pca %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    across(
      starts_with("PC"),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  )

pc_summary

tract_pca %>%
  pivot_longer(
    starts_with("PC"),
    names_to = "Component",
    values_to = "Score"
  ) %>%
  ggplot(aes(x = factor(cluster), y = Score, fill = factor(cluster))) +
  geom_boxplot(outlier.size = 0.3) +
  facet_wrap(~Component, scales = "free_y") +
  labs(
    x = "Heat Cluster",
    y = "PCA Score"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_blank()
  )
#---- 3. Stack #----
library(dplyr)

pc_vars <- c(
  "PC1_urban",
  "PC2_socioeconomic",
  "PC3_age_structure",
  "PC4_AC",
  "PC5_greenness",
  "PC6_low_density",
  "PC7_water"
)

tract_pca_scored <- tract_pca %>%
  mutate(
    across(
      all_of(pc_vars),
      ~ {
        pc_mean <- mean(.x, na.rm = TRUE)
        pc_sd   <- sd(.x, na.rm = TRUE)
        
        case_when(
          is.na(.x)                    ~ NA_integer_,
          .x <= pc_mean - 2 * pc_sd    ~ 1L,
          .x <= pc_mean - 1 * pc_sd    ~ 2L,
          .x <  pc_mean                ~ 3L,
          .x <  pc_mean + 1 * pc_sd    ~ 4L,
          .x <  pc_mean + 2 * pc_sd    ~ 5L,
          .x >= pc_mean + 2 * pc_sd    ~ 6L
        )
      },
      .names = "{.col}_score"
    )
  )
tract_pca_scored %>%
  select(ends_with("_score")) %>%
  summarise(
    across(
      everything(),
      ~ paste(names(table(.x, useNA = "ifany")),
              table(.x, useNA = "ifany"),
              sep = ": ",
              collapse = "; ")
    )
  ) %>%
  pivot_longer(
    everything(),
    names_to = "PC",
    values_to = "distribution"
  )
pc_cluster_score <- tract_pca_scored %>%
  filter(!is.na(cluster)) %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    across(
      ends_with("_score"),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

pc_cluster_score
library(tidyr)
library(ggplot2)

pc_cluster_long <- pc_cluster_score %>%
  select(-n) %>%
  pivot_longer(
    cols = ends_with("_score"),
    names_to = "component",
    values_to = "mean_score"
  ) %>%
  mutate(
    component = factor(
      component,
      levels = paste0(pc_vars, "_score"),
      labels = c(
        "Urbanization",
        "Socioeconomic disadvantage",
        "Age structure",
        "AC vulnerability",
        "Greenness and reflectance",
        "Low-density development",
        "Surface water"
      )
    ),
    cluster = factor(cluster)
  )

p_pc_stack <- ggplot(
  pc_cluster_long,
  aes(
    x = cluster,
    y = mean_score,
    fill = component
  )
) +
  geom_col(
    width = 0.75,
    color = "white",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = sprintf("%.2f", mean_score)),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  coord_flip() +
  scale_fill_brewer(
    palette = "Set2",
    name = "PCA component"
  ) +
  labs(
    x = "Heat cluster",
    y = "Sum of mean assigned PCA scores",
    title = "PCA profiles by heat cluster",
    subtitle = "Original PCA scores categorized from 1 to 6 using standard-deviation thresholds"
  ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

p_pc_stack
#---- 3. Clustering heatwave characteristics #----

heat_char[is.na(heat_char)] <- 0
# heat_char <- heat_char %>%
#   dplyr::select(
#     StCoFIPS,
# 
#     hw3_mean_frequency, hw3_mean_duration, hw3_mean_intensity,
#     hw3_max_frequency, hw3_max_duration, hw3_max_intensity,
#     hw3_min_frequency, hw3_min_duration, hw3_min_intensity
#   )
# 
# heat_char <- heat_char %>%
#   dplyr::select(
#     StCoFIPS,
# 
#     hw3_min_frequency, hw3_min_duration, hw3_min_intensity,
#     hw3_mean_frequency, hw3_mean_duration, hw3_mean_intensity,
#     hw3_max_frequency, hw3_max_duration, hw3_max_intensity
#   )

heat_scaled <- heat_char %>%
  column_to_rownames("StCoFIPS") %>%
  scale()
ev <- eigen(cor(heat_scaled))
ev$values
# eigen>1; 6 clusters
pc_heat <- principal(
  heat_scaled,
  nfactors = 6,
  rotate = "varimax",
  scores = TRUE
)
print(
  pc_heat,
  digits = 3,
  cutoff = 0.30,
  sort = TRUE
)

scores <- as.data.frame(pc_heat$scores)

cor(scores)

# Elbow
set.seed(123)

wss <- sapply(1:10, function(k){
  
  kmeans(
    scores,
    centers = k,
    nstart = 100
  )$tot.withinss
  
})

elbow <- data.frame(
  k = 1:10,
  wss = wss
)
ggplot(elbow, aes(k, wss)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = 1:10) +
  theme_bw()

# Silhouette analysis
sil <- sapply(2:10, function(k){
  km <- kmeans(
    scores,
    centers = k,
    nstart = 100
  )
  mean(
    silhouette(km$cluster,
               dist(scores))[,3]
  )
})
data.frame(
  k = 2:10,
  silhouette = sil
) |>
  ggplot(aes(k, silhouette)) +
  geom_line() +
  geom_point(size=3) +
  theme_bw()

# NbClust

set.seed(123)

nb <- NbClust(
  data = scores,
  distance = "euclidean",
  min.nc = 2,
  max.nc = 10,
  method = "kmeans",
  index = "all"
)
nb$Best.nc


km <- kmeans(
  scores,
  centers = 4,
  nstart = 500
)

heat_char$cluster <- factor(km$cluster)

cluster_profile <- scores |>
  mutate(cluster = factor(km$cluster)) |>
  group_by(cluster) |>
  summarise(across(everything(), mean))
cluster_profile

library(pheatmap)
centers <- cluster_profile %>%
  column_to_rownames("cluster")

pheatmap(
  centers,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  color = colorRampPalette(c("blue","white","red"))(100)
)





options(tigris_use_cache = TRUE)

counties_sf <- counties(cb = TRUE, year = 2020) %>%
  st_transform(5070)

conus <- counties_sf %>%
  filter(
    !STATEFP %in% c("02", "15", "60", "66", "69", "72", "78")
  )

heat_char <- heat_char %>%
  mutate(
    StCoFIPS = sprintf("%05s", StCoFIPS)
  )

conus_cluster <- conus %>%
  left_join(
    heat_char,
    by = c("GEOID" = "StCoFIPS")
  )

ggplot(conus_cluster) +
  geom_sf(aes(fill = factor(cluster)), color = NA) +
  scale_fill_brewer(palette = "Set2", name = "Cluster") +
  labs(
    title = "Heatwave Pattern Clusters Across U.S. Counties"
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold")
  )

hw3_summary <- heat_char %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    
    frequency_mean = mean(hw3_mean_frequency, na.rm = TRUE),
    frequency_sd   = sd(hw3_mean_frequency, na.rm = TRUE),
    
    duration_mean  = mean(hw3_mean_duration, na.rm = TRUE),
    duration_sd    = sd(hw3_mean_duration, na.rm = TRUE),
    
    intensity_mean = mean(hw3_mean_intensity, na.rm = TRUE),
    intensity_sd   = sd(hw3_mean_intensity, na.rm = TRUE),
    
    .groups = "drop"
  )

hw3_summary
#---- 4. Percentile classification #----
heat_scaled <- heat_char %>%
  mutate(across(-StCoFIPS, scale))

heat_scores <- heat_scaled %>%
  transmute(
    StCoFIPS,
    
    frequency_score = rowMeans(
      dplyr::select(., ends_with("frequency")),
      na.rm = TRUE
    ),
    
    duration_score = rowMeans(
      dplyr::select(., ends_with("duration")),
      na.rm = TRUE
    ),
    
    intensity_score = rowMeans(
      dplyr::select(., ends_with("intensity")),
      na.rm = TRUE
    )
  )
##---- 3.1 Mapping scores #----
options(tigris_use_cache = TRUE)

counties <- counties(
  cb = FALSE,
  year = 2020,
  class = "sf"
)

counties <- counties %>%
  filter(
    !STATEFP %in% c("02", "15", "60", "66", "69", "72", "78")
  ) %>%
  mutate(
    StCoFIPS = GEOID
  )

counties_heat <- counties %>%
  left_join(heat_scores, by = "StCoFIPS")



set.seed(1234)

counties_heat <- counties_heat %>%
  filter(
    !is.na(frequency_score),
    !is.na(duration_score),
    !is.na(intensity_score),
    is.finite(frequency_score),
    is.finite(duration_score),
    is.finite(intensity_score)
  )
cluster3 <- function(x){
  
  km <- kmeans(
    x = matrix(x, ncol = 1),
    centers = 3,
    nstart = 100
  )
  
  centers <- km$centers[,1]
  
  ord <- order(centers)
  
  lab <- rep(NA_character_, 3)
  lab[ord] <- c("Low", "Medium", "High")
  
  factor(
    lab[km$cluster],
    levels = c("Low","Medium","High")
  )
}


counties_heat <- counties_heat %>%
  mutate(
    freq_group = cluster3(frequency_score),
    dur_group  = cluster3(duration_score),
    int_group  = cluster3(intensity_score)
  )

counties_heat <- counties_heat %>%
  mutate(
    heat_profile = paste(
      int_group,
      freq_group,
      dur_group,
      
      sep = "_"
    )
  )


heat_profile_summary <- counties_heat %>%
  st_drop_geometry() %>%
  count(
    int_group,
    freq_group,
    dur_group,
    
    heat_profile,
    sort = TRUE
  )

print(heat_profile_summary)


table(counties_heat$heat_profile)

with(
  st_drop_geometry(counties_heat),
  table(int_group, freq_group, dur_group)
)

hue_vals <- c(
  Low = 220,      # blue
  Medium = 40,    # orange
  High = 0        # red
)

chroma_vals <- c(
  Low = 30,
  Medium = 55,
  High = 80
)

luminance_vals <- c(
  Low = 90,
  Medium = 70,
  High = 50
)

palette27 <-
  expand.grid(
    int_group  = c("Low","Medium","High"),
    freq_group = c("Low","Medium","High"),
    dur_group  = c("Low","Medium","High")
  ) %>%
  mutate(
    heat_profile = paste(int_group,freq_group, dur_group,  sep="_"),
    color = hcl(
      h = hue_vals[int_group],
      c = chroma_vals[freq_group],
      l = luminance_vals[dur_group]
    )
  )


counties_heat <- counties_heat %>%
  left_join(
    palette27 %>%
      dplyr::select(heat_profile, color),
    by="heat_profile"
  )
tmap_mode("plot")


map_heat <- tm_shape(counties_heat) +
  tm_polygons(
    col = "heat_profile",
    palette = setNames(
      palette27$color,
      palette27$heat_profile
    ),
    border.col = NA
  ) +
  tm_layout(
    frame = FALSE,
    legend.outside = TRUE,
    legend.outside.position = "right"
  )

tmap_save(
  map_heat,
  filename = "Output/Figures/Heatwave_Profile_Map.jpg",
  width = 12,
  height = 7,
  units = "in",
  dpi = 600
)

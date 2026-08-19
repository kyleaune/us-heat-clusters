#---- Script Metadata #----
# Title: Heatwave clustering method comparison
# Author: Hang Li
# Purpose: Compare four ways of grouping counties by heatwave characteristics:
#   1) ordinary k-means (non-spatial, same approach as heatwave_characteristics.R)
#   2) Ward-like hierarchical clustering with a spatial-contiguity constraint,
#      via ClustGeo::hclustgeo()
#   3) ordinary hierarchical clustering (Ward, no spatial constraint) - a
#      second, independent non-spatial baseline
#   4) IECC 2021 building-energy climate zones - an existing, non-data-driven
#      classification, used as an external reference rather than fit here
# All three fitted methods cluster on the SAME variable portfolio and the SAME
# county set, so differences between them isolate the effect of the spatial
# constraint (method 2) and of the algorithm itself (method 3 vs 1).
#
# Input: heat_characteristics_v2.csv, already computed by
# heatwave_characteristics.R's "#---- Heatwave characteristics #----" section.
# Nothing here re-reads the 1.1 GB raw WBGT file or recomputes heatwave season/
# spell logic - that is done once, upstream, in that script.
#---------------------------------------------------------------------------


#===========================================================================
# 1. FILE INPUT AND SETUP
#===========================================================================
rm(list = ls())
pkgs <- c("dplyr", "readr", "tibble", "tidyr", "ggplot2",
          "cluster", "NbClust", "parallel", "doParallel", "foreach",
          "sf", "tigris", "spdep", "ClustGeo", "fpc", "clevr",
          "terra", "exactextractr", "officer", "flextable")
install.packages(setdiff(pkgs, rownames(installed.packages())))
invisible(lapply(pkgs, library, character.only = TRUE))
library(patchwork)
data_dir   <- "/Users/hangli/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Aune, Kyle T's files - 260204 - US Heat Clusters/Data"
season_dir <- file.path(data_dir, "processed", "season cluster")
fig_dir    <- "/Users/hangli/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Aune, Kyle T's files - 260204 - US Heat Clusters/Output/Figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

options(tigris_use_cache = TRUE)

## --- heatwave characteristics, already computed ----------------------------
county_heat_characteristics <- read_csv(
  file.path(season_dir, "heat_characteristics_v2.csv"),
  col_types = cols(StCoFIPS = col_character(), .default = col_double())
)

print(setdiff(names(county_heat_characteristics), "StCoFIPS"))

## --- county geometry, needed only by section 3's spatial constraint --------

counties_sf <- suppressMessages(
  counties(cb = TRUE, year = 2020, class = "sf", progress_bar = FALSE)
) %>%
  mutate(GEOID = as.character(GEOID)) %>%
  filter(!STATEFP %in% c("02", "15", "60", "66", "69", "72", "78")) %>%
  # AK, HI, American Samoa, Guam, N. Mariana Is., Puerto Rico, US Virgin Is.
  st_transform(5070)

cat("\nCONUS counties with geometry:", nrow(counties_sf), "\n")

# Queen contiguity, built once and reused by section 3's k-selection and by
# the ClustGeo/k-means comparison - poly2nb is the most expensive step in the
# fragmentation diagnostic below, so it should not be rebuilt on every call.
nb          <- poly2nb(counties_sf, queen = TRUE)
geoid_order <- counties_sf$GEOID

# Given a cluster label per county, sever every adjacency edge that crosses a
# cluster boundary and count the surviving connected components. This is the
# spatial-coherence half of evaluating a clustering - silhouette only sees the
# feature space and cannot tell a compact region from one scattered across the
# country.
count_fragments <- function(county_ids, cluster_vec) {
  cl   <- cluster_vec[match(geoid_order, county_ids)]
  keep <- !is.na(cl)

  nb_sub <- subset(nb, keep)
  cl_sub <- cl[keep]

  nb_cut <- nb_sub
  for (i in seq_along(nb_cut)) {
    if (identical(nb_cut[[i]], 0L)) next
    same <- nb_cut[[i]][cl_sub[nb_cut[[i]]] == cl_sub[i]]
    nb_cut[[i]] <- if (length(same) == 0) 0L else same
  }

  comp <- n.comp.nb(nb_cut)
  list(nc = comp$nc, comp_id = comp$comp.id, cluster = cl_sub, n_county = sum(keep))
}

# Same anchor-color ramp used for cluster maps elsewhere in this project
# (heatwave_characteristics.R's mapping section), duplicated here as a small,
# stable utility rather than sourcing that whole (heavy) script just for one
# helper function.
severity_anchors <- c("#FFFFB2", "#FECC5C", "#FD8D3C", "#F03B20",
                      "#BD0026", "#7A0177", "#49006A")
build_severity_pal <- function(k) {
  stats::setNames(grDevices::colorRampPalette(severity_anchors, space = "Lab")(k),
                  as.character(seq_len(k)))
}

# Ranks cluster labels 1 (coolest) .. k (hottest) by mean rel_intensity_tmax,
# the same single-variable rule used in heatwave_characteristics.R step 5b.
# Shared by every method below so all of them are labelled the same way and
# read directly comparably.
relabel_by_tmax <- function(county_ids, raw_cluster, heat_dat, k) {
  key <- tibble(StCoFIPS = county_ids, cluster_raw = raw_cluster) %>%
    left_join(heat_dat %>% dplyr::select(StCoFIPS, rel_intensity_tmax), by = "StCoFIPS") %>%
    group_by(cluster_raw) %>%
    summarise(mean_tmax = mean(rel_intensity_tmax, na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_tmax) %>%
    transmute(cluster_raw, cluster = row_number())
  stopifnot(nrow(key) == k)
  cat("  crosswalk (raw id -> coolest..hottest by mean rel_intensity_tmax):\n")
  print(as.data.frame(key), row.names = FALSE)
  factor(key$cluster[match(raw_cluster, key$cluster_raw)], levels = seq_len(k))
}

# CH (Calinski-Harabasz) index for any grouping vector, backed by
# fpc::calinhara() rather than a hand-rolled formula (the hand-rolled version
# used earlier in this file's development matched fpc::calinhara() to float
# precision on synthetic data, so there is no re-derivation risk in switching).
#
# fpc::calinhara(x, clustering, cn = max(clustering)) is only correct if group
# labels are consecutive integers 1..k - confirmed by testing: with
# non-consecutive labels (e.g. 1, 3, 5) the default cn silently returns the
# WRONG value (it computed 194.9 instead of the correct 490.3 in that test).
# Converting to a factor and taking as.integer() first guarantees consecutive
# labels regardless of the input's original coding, and cn is always passed
# explicitly rather than relying on the default.
ch_index <- function(x, group) {
  group <- droplevels(as.factor(group))
  fpc::calinhara(x, as.integer(group), cn = nlevels(group))
}


#===========================================================================
# 2. K-MEANS METHOD (non-spatial baseline, mirrors heatwave_characteristics.R)
#===========================================================================

##---- 2.1 variable portfolio #----
# Kept identical to cluster_vars in heatwave_characteristics.R so the two
# scripts' k-means results match. If that portfolio changes there, update it
# here too.
cluster_vars <- c(
  "abs_hwdays_per_yr_28",
  "rel_intensity_tmax",
  "rel_duration_max",
  "rel_hw_span"
)
cluster_input <- county_heat_characteristics %>%
  dplyr::select(StCoFIPS, dplyr::all_of(cluster_vars)) %>%
  na.omit()

cluster_scaled <- scale(cluster_input[-1])

##---- 2.2 profile: elbow, silhouette, CH index #----
# Three views of the same k-means fits, one kmeans() call per k rather than
# one loop per metric. WSS and betweenss come straight off the kmeans object;
# CH now goes through the shared ch_index() helper (section 1) instead of its
# own inline formula, so every method in this file uses one CH implementation.
K_RANGE <- 2:10
d_full  <- dist(cluster_scaled)
n_obs   <- nrow(cluster_scaled)

profile_km <- t(sapply(K_RANGE, function(k) {
  km_k <- kmeans(cluster_scaled, centers = k, nstart = 50)
  c(wss        = km_k$tot.withinss,
    silhouette = mean(silhouette(km_k$cluster, d_full)[, 3]),
    ch         = ch_index(cluster_scaled, km_k$cluster))
}))
profile_km <- data.frame(k = K_RANGE, profile_km)

cat("\nk-means profile (", length(cluster_vars), "-variable portfolio):\n", sep = "")
print(round(profile_km, 3), row.names = FALSE)
cat("silhouette peaks at k =", K_RANGE[which.max(profile_km$silhouette)],
    "| CH peaks at k =", K_RANGE[which.max(profile_km$ch)], "\n")

p1 <- ggplot(profile_km, aes(k, wss)) +
  geom_line() + geom_point() +
  labs(x = "k", y = "Total within-cluster SS", title = "(a) Elbow") +
  theme_minimal(base_family = "Arial", base_size = 8)

p2 <- ggplot(profile_km, aes(k, silhouette)) +
  geom_line() + geom_point() +
  labs(x = "k", y = "Mean silhouette width", title = "(b) Silhouette") +
  theme_minimal(base_family = "Arial", base_size = 8)

p3 <- ggplot(profile_km, aes(k, ch)) +
  geom_line() + geom_point() +
  labs(x = "k", y = "Calinski-Harabasz index", title = "(c) CH index") +
  theme_minimal(base_family = "Arial", base_size = 8)

combined <- p1 + p2 + p3
combined
ggsave(file.path(fig_dir,"k_selection_diagnostics.png"), combined, width = 6.5, height = 3, dpi = 300)

##---- 2.3 final k-means #----
optimal_k <- 7

set.seed(123)
km <- kmeans(cluster_scaled, centers = optimal_k, nstart = 100)

sil_km <- silhouette(km$cluster, d_full)
cat("\nk-means, k =", optimal_k, "- mean silhouette width:",
    round(mean(sil_km[, 3]), 4), "\n")

##---- 2.4 relabel clusters by ascending tmax #----
# k-means cluster IDs are arbitrary. Same rule as heatwave_characteristics.R
# step 5b: rank by mean rel_intensity_tmax alone (ascending) and relabel
# 1 = coolest ... optimal_k = hottest.
cat("\nk-means cluster relabelling:\n")
cluster_kmeans <- cluster_input %>%
  dplyr::select(StCoFIPS) %>%
  mutate(
    cluster_raw = km$cluster,
    cluster     = relabel_by_tmax(StCoFIPS, km$cluster, county_heat_characteristics, optimal_k)
  )

# write_csv(cluster_kmeans, file.path(season_dir, "clusterCompare_kmeans.csv"))


#===========================================================================
# 3. WARD-LIKE SPATIALLY-CONSTRAINED CLUSTERING (ClustGeo)
#===========================================================================
# hclustgeo() minimizes a convex combination of two Ward-type criteria: one on
# a "feature" dissimilarity D0 (same as k-means clusters on) and one on a
# "constraint" dissimilarity D1 (here, geographic distance between county
# centroids). alpha in [0,1] sets the mix: alpha = 0 ignores geography
# entirely (Ward clustering on D0 alone); alpha = 1 ignores the heatwave
# variables and clusters purely on geographic proximity.
#
# D0 and D1 must describe the SAME counties in the SAME order - built from one
# joined table below rather than two separately filtered objects, so there is
# no chance of them silently going out of alignment.
#
# Order of steps: D0/D1 -> choose alpha -> FIT the tree -> choose k on that
# fitted tree -> cut at the chosen k -> relabel. k is chosen on an
# already-fitted tree (cutree() is nearly free once hclustgeo() has run), so
# k-selection sits between fitting and the final cut, not after it.

##---- 3.1 build D0 (feature space) and D1 (geographic space) #----
county_pts <- counties_sf %>%
  st_point_on_surface() %>%     # guaranteed inside the polygon, unlike centroid
  st_coordinates() %>%
  as.data.frame() %>%
  mutate(StCoFIPS = counties_sf$GEOID) %>%
  rename(x = X, y = Y)

geo_input <- cluster_input %>%
  inner_join(county_pts, by = "StCoFIPS")

n_dropped <- nrow(cluster_input) - nrow(geo_input)
if (n_dropped > 0) {
  message(n_dropped, " counties in cluster_input have no matching geometry ",
          "and are dropped from the ClustGeo comparison")
}

D0 <- dist(scale(geo_input[, cluster_vars]))
D1 <- dist(geo_input[, c("x", "y")])

stopifnot(attr(D0, "Size") == attr(D1, "Size"),
          attr(D0, "Size") == nrow(geo_input))

##---- 3.2 choose alpha #----
# D1 is in metres (up to ~4,000,000 across CONUS) and D0 is in z-score units,
# so scale = TRUE (the default) is essential here - without it D1 would
# dominate the combination at every alpha. choicealpha() rescales both to
# [0,1] before mixing and plots how much of each criterion's inertia survives
# a K-cluster partition across the alpha grid; the usual choice is the alpha
# just past where the D1 (spatial) curve stops rising steeply, without giving
# up much of the D0 (feature) curve. K here is a reference value only (matched
# to k-means' optimal_k), used solely to evaluate the alpha tradeoff - the
# actual k for the final ClustGeo partition is chosen in 3.4, after the tree
# below is fit.
range_alpha <- seq(0, 1, 0.1)

alpha_choice <- choicealpha(D0, D1, range.alpha = range_alpha,
                            K = optimal_k, scale = TRUE, graph = FALSE)

cat("\nchoicealpha: proportion of explained inertia (Q) by alpha, reference K =", optimal_k, "\n")
print(data.frame(alpha = range_alpha,
                 Q0_feature = round(alpha_choice$Q[, 1], 3),
                 Q1_spatial = round(alpha_choice$Q[, 2], 3)))

plot(alpha_choice, main = paste("choicealpha, reference K =", optimal_k))

# Pick alpha from the plot/table above, then re-run from here down.
alpha_geo <- 0.2

##---- 3.3 fit the tree #----
# Fitting only - no cutree() yet. hclustgeo() builds the full dendrogram
# regardless of k, so k is chosen (3.4) on this one fitted tree before the
# final cut (3.5).
tree_geo <- hclustgeo(D0, D1, alpha = alpha_geo, scale = TRUE)

##---- 3.4 choosing k for ClustGeo #----
# tree_geo already holds the full dendrogram at alpha_geo - unlike k-means,
# trying a different k costs nothing further; it is just cutree() at a
# different height. Four diagnostics, read together rather than trusting one:
#
#   a) fusion-height elbow    cheapest - uses the tree as already fit
#   b) D0 silhouette across k feature-side fit, same lens as section 2.2
#   c) fragmentation across k spatial-side coherence, which (a)/(b) cannot see
#                              at all: silhouette on D0 cannot tell a compact
#                              region from one scattered across the country
#   d) alpha = 0 profile      hclustgeo(D0, alpha = 0) IS ordinary Ward on D0
#                              (no geography at all), so this is a free
#                              cross-check against k-means' own silhouette
#                              profile in section 2.2 - if they diverge a lot,
#                              Ward and k-means disagree on structure before
#                              geography even enters the picture. (Section 4
#                              below fits the same ordinary-Ward idea as its
#                              own full method, on the complete county set.)
#
# alpha was chosen (3.2) at a fixed reference K to match k-means for
# comparability. If the diagnostics below argue for a different k, that alpha
# choice technically ought to be revisited too - choicealpha()'s Q0/Q1
# tradeoff can shift with K. Re-running 3.2 with K set to the new k is the
# lightweight check; a full alpha-by-k grid is the heavier, rarely-necessary
# version of that check.

# (a) fusion-height elbow: last N merge heights, most recent merge first
n_show  <- 15
heights <- rev(tail(tree_geo$height, n_show))
plot(seq_len(n_show), heights, type = "b", pch = 19,
     xlab = "merges before the final partition (1 = last)", ylab = "fusion height",
     main = sprintf("ClustGeo fusion heights, alpha = %.1f", alpha_geo))

# (b) + (c): same k range, cut from the one already-fit tree
K_RANGE_GEO <- 2:10

sil_profile_geo <- sapply(K_RANGE_GEO, function(k) {
  mean(silhouette(cutree(tree_geo, k = k), D0)[, 3])
})

frag_profile_geo <- sapply(K_RANGE_GEO, function(k) {
  count_fragments(geo_input$StCoFIPS, cutree(tree_geo, k = k))$nc
})

# (d) ordinary Ward on D0 alone, for comparison with k-means' silhouette curve
tree_ward0 <- hclustgeo(D0, alpha = 0)
sil_profile_ward0 <- sapply(K_RANGE_GEO, function(k) {
  mean(silhouette(cutree(tree_ward0, k = k), D0)[, 3])
})

k_profile_geo <- data.frame(
  k                  = K_RANGE_GEO,
  silhouette_alpha   = round(sil_profile_geo, 4),
  silhouette_ward0   = round(sil_profile_ward0, 4),
  pieces_total       = frag_profile_geo,
  pieces_per_cluster = round(frag_profile_geo / K_RANGE_GEO, 1)
)
cat("\nClustGeo k profile (alpha =", alpha_geo, "):\n")
print(k_profile_geo, row.names = FALSE)

par(mfrow = c(1, 2))
plot(K_RANGE_GEO, sil_profile_geo, type = "b", pch = 19, col = "firebrick",
     ylim = range(c(sil_profile_geo, sil_profile_ward0)),
     xlab = "k", ylab = "mean silhouette (on D0)",
     main = sprintf("Silhouette: alpha=%.1f vs ordinary Ward", alpha_geo))
lines(K_RANGE_GEO, sil_profile_ward0, type = "b", pch = 17, col = "grey40", lty = 2)
legend("topright", legend = c(sprintf("alpha=%.1f", alpha_geo), "alpha=0 (ordinary Ward)"),
       col = c("firebrick", "grey40"), pch = c(19, 17), lty = c(1, 2),
       bty = "n", cex = 0.8)
plot(K_RANGE_GEO, k_profile_geo$pieces_per_cluster, type = "b", pch = 19,
     xlab = "k", ylab = "fragments per cluster",
     main = sprintf("Spatial fragmentation, alpha=%.1f", alpha_geo))
par(mfrow = c(1, 1))

#write_csv(k_profile_geo, file.path(season_dir, "clustgeo_k_profile.csv"))

# Pick k for the final ClustGeo partition from 3.4 above. Defaults to
# optimal_k (matching k-means) so the three methods stay directly comparable
# unless the diagnostics argue otherwise.

optimal_k_geo <- optimal_k

##---- 3.5 cut the tree and relabel by ascending tmax #----
raw_cut_geo <- cutree(tree_geo, k = optimal_k_geo)

cat("\nClustGeo cluster relabelling:\n")
cluster_clustgeo <- geo_input %>%
  dplyr::select(StCoFIPS) %>%
  mutate(
    cluster_raw = raw_cut_geo,
    cluster     = relabel_by_tmax(StCoFIPS, raw_cut_geo, county_heat_characteristics, optimal_k_geo)
  )

cat("\nClustGeo, alpha =", alpha_geo, ", k =", optimal_k_geo, "- cluster sizes:\n")
print(table(cluster_clustgeo$cluster))

# write_csv(cluster_clustgeo, file.path(season_dir, "clusterCompare_clustgeo.csv"))


#===========================================================================
# 4. ORDINARY HIERARCHICAL CLUSTERING (Ward, no spatial constraint)
#===========================================================================
# A second non-spatial baseline, distinct in implementation from k-means
# (partition-based, refit per k) and from ClustGeo (spatial constraint). This
# is the same "ordinary Ward on D0" idea already used as a cross-check inside
# section 3.4(d), but fit here as its own full method: base R hclust(), on
# the COMPLETE cluster_scaled county set (not geo_input, which is a few
# counties smaller wherever county geometry was unavailable) - no coordinates
# are needed for this method, so there is no reason to accept that loss here.

##---- 4.1 fit the tree #----
D_hc      <- dist(cluster_scaled)
tree_hc   <- hclust(D_hc, method = "ward.D2")

##---- 4.2 choosing k #----
# Same three-part logic as sections 2.2 (k-means) and 3.4 (ClustGeo): fusion
# height, silhouette, CH index, all read from the one already-fitted tree.
K_RANGE_HC <- 2:10

sil_profile_hc <- sapply(K_RANGE_HC, function(k) {
  mean(silhouette(cutree(tree_hc, k = k), D_hc)[, 3])
})
ch_profile_hc <- sapply(K_RANGE_HC, function(k) {
  ch_index(cluster_scaled, cutree(tree_hc, k = k))
})

k_profile_hc <- data.frame(
  k          = K_RANGE_HC,
  silhouette = round(sil_profile_hc, 4),
  ch         = round(ch_profile_hc, 1)
)
cat("\nhierarchical (Ward) k profile:\n")
print(k_profile_hc, row.names = FALSE)
cat("silhouette peaks at k =", K_RANGE_HC[which.max(sil_profile_hc)],
    "| CH peaks at k =", K_RANGE_HC[which.max(ch_profile_hc)], "\n")

n_show_hc  <- 15
heights_hc <- rev(tail(tree_hc$height, n_show_hc))

par(mfrow = c(1, 3))
plot(seq_len(n_show_hc), heights_hc, type = "b", pch = 19,
     xlab = "merges before the final partition (1 = last)", ylab = "fusion height",
     main = "1. Elbow (fusion height)")
plot(K_RANGE_HC, sil_profile_hc, type = "b", pch = 19,
     xlab = "k", ylab = "Mean silhouette width", main = "2. Silhouette")
plot(K_RANGE_HC, ch_profile_hc, type = "b", pch = 19,
     xlab = "k", ylab = "Calinski-Harabasz index", main = "3. CH index")
par(mfrow = c(1, 1))

optimal_k_hc <- 7

##---- 4.3 cut the tree and relabel by ascending tmax #----
raw_cut_hc <- cutree(tree_hc, k = optimal_k_hc)

cat("\nHierarchical (Ward) cluster relabelling:\n")
cluster_hclust <- cluster_input %>%
  dplyr::select(StCoFIPS) %>%
  mutate(
    cluster_raw = raw_cut_hc,
    cluster     = relabel_by_tmax(StCoFIPS, raw_cut_hc, county_heat_characteristics, optimal_k_hc)
  )

cat("\nHierarchical (Ward), k =", optimal_k_hc, "- cluster sizes:\n")
print(table(cluster_hclust$cluster))

# write_csv(cluster_hclust, file.path(season_dir, "clusterCompare_hclust.csv"))


#===========================================================================
# 5. COMPARISON WITH EXTERNAL CLIMATE/REGION CLASSIFICATIONS
#===========================================================================
# How much does each data-driven clustering (k-means, ClustGeo, ordinary
# hierarchical) agree with FOUR existing, non-data-driven reference systems,
# and does any of them separate the heat-metric feature space better than
# those references do?
#   1) IECC 2021        building-energy climate zones (county, shapefile)
#   2) Koppen-Geiger     1991-2020 climate classification (global raster,
#                        majority class per county)
#   3) IPCC WGI          AR6 reference regions (global polygons, county
#                        centroid membership)
#   4) NOAA 9            NCEI's 9 climate regions (state-based, not spatial -
#                        every county in a state gets that state's region)
#
#   V-measure  - agreement between a clustering and a reference's labels
#                (harmonic mean of homogeneity and completeness; 0 = no
#                shared structure, 1 = identical partitions)
#   CH index   - cluster quality of EACH labelling (3 methods + 4 references)
#                on the SAME feature space, so a higher CH for one of our
#                clusterings than for a reference means the four heatwave
#                metrics separate better under that clustering than under
#                that reference system.
#
# clevr::v_measure() is used rather than a hand-rolled version - the
# hand-rolled entropy formula was checked against it on synthetic data and
# matched to float precision, so there is no re-derivation risk here, and
# using the package function is less code to maintain.

##---- 5.1 load IECC21, Koppen-Geiger, IPCC WGI regions, NOAA9 and build one common county set #----

# --- IECC21 (county, shapefile) ---------------------------------------------
# GEOID in this file is "G" + 5-digit FIPS (e.g. "G39057"), not the bare FIPS
# used everywhere else in this script - stripped here at the point of entry so
# every join below uses the same StCoFIPS convention as cluster_kmeans etc.
climate_zones <- st_read(file.path(data_dir, "ClimateZoneDataFiles/ClimateZones.shp"),
                         quiet = TRUE) %>%
  st_drop_geometry() %>%
  transmute(StCoFIPS = sub("^G", "", GEOID), iecc21 = IECC21) %>%
  filter(!is.na(iecc21))

# --- Koppen-Geiger (global raster, majority class per county) --------------
# Raster is EPSG:4326 (lon/lat); counties_sf is already projected to 5070 for
# the ClustGeo section, so reprojecting the (smaller) county polygons to match
# the raster is cheaper than reprojecting the raster itself.
kg_rast     <- terra::rast(file.path(
  data_dir, "ClimateZoneDataFiles/koppen_geiger_tif/1991_2020/koppen_geiger_0p00833333.tif"
))
counties_ll <- st_transform(counties_sf, sf::st_crs(kg_rast))
kg_mode     <- exactextractr::exact_extract(kg_rast, counties_ll, fun = "mode")

koppen_zones <- tibble(StCoFIPS = counties_ll$GEOID, koppen = kg_mode) %>%
  filter(!is.na(koppen))

# --- IPCC WGI AR6 reference regions (global polygons, county centroid) -----
ipcc_regions <- st_read(file.path(data_dir, "ClimateZoneDataFiles/IPCC-WGI-reference-regions-v4.geojson"),
                        quiet = TRUE) %>%
  st_transform(5070) %>%
  dplyr::select(Acronym)

ipcc_zones <- counties_sf %>%
  st_point_on_surface() %>%
  dplyr::select(GEOID) %>%
  st_join(ipcc_regions, join = st_within, left = FALSE) %>%
  st_drop_geometry() %>%
  transmute(StCoFIPS = GEOID, ipcc = Acronym)

# --- NOAA 9 climate regions (state-based, not spatial) ----------------------
# The state index numbers in NOAA's own region tables (e.g. "Connecticut (6)")
# are NOT Census state FIPS codes - they are the 48 CONUS states in
# ALPHABETICAL order (Alabama = 1 ... Wyoming = 48, skipping AK/HI), which is
# how NOAA indexes states in its regional-climate tables, not how Census
# TIGER data (used everywhere else in this project, via STATEFP) does. The
# join below therefore goes through STATE NAME, not that index number, to the
# standard Census STATEFP that counties_sf already carries.
noaa9_states <- tibble::tribble(
  ~state_name,      ~noaa9_region,
  "Connecticut",    "Northeast",
  "Delaware",       "Northeast",
  "Maine",          "Northeast",
  "Maryland",       "Northeast",
  "Massachusetts",  "Northeast",
  "New Hampshire",  "Northeast",
  "New Jersey",     "Northeast",
  "New York",       "Northeast",
  "Pennsylvania",   "Northeast",
  "Rhode Island",   "Northeast",
  "Vermont",        "Northeast",
  "Iowa",           "Upper Midwest",
  "Michigan",       "Upper Midwest",
  "Minnesota",      "Upper Midwest",
  "Wisconsin",      "Upper Midwest",
  "Illinois",       "Ohio Valley",
  "Indiana",        "Ohio Valley",
  "Kentucky",       "Ohio Valley",
  "Missouri",       "Ohio Valley",
  "Ohio",           "Ohio Valley",
  "Tennessee",      "Ohio Valley",
  "West Virginia",  "Ohio Valley",
  "Alabama",        "Southeast",
  "Florida",        "Southeast",
  "Georgia",        "Southeast",
  "North Carolina", "Southeast",
  "South Carolina", "Southeast",
  "Virginia",       "Southeast",
  "Montana",        "Northern Rockies and Plains",
  "Nebraska",       "Northern Rockies and Plains",
  "North Dakota",   "Northern Rockies and Plains",
  "South Dakota",   "Northern Rockies and Plains",
  "Wyoming",        "Northern Rockies and Plains",
  "Arkansas",       "South",
  "Kansas",         "South",
  "Louisiana",      "South",
  "Mississippi",    "South",
  "Oklahoma",       "South",
  "Texas",          "South",
  "Arizona",        "Southwest",
  "Colorado",       "Southwest",
  "New Mexico",     "Southwest",
  "Utah",           "Southwest",
  "Idaho",          "Northwest",
  "Oregon",         "Northwest",
  "Washington",     "Northwest",
  "California",     "West",
  "Nevada",         "West"
)

# Standard Census state FIPS (STATEFP) - stable public reference, distinct
# from NOAA's own alphabetical state index above.
state_fips <- tibble::tribble(
  ~state_name,      ~STATEFP,
  "Alabama",        "01", "Arizona",         "04", "Arkansas",       "05",
  "California",     "06", "Colorado",        "08", "Connecticut",    "09",
  "Delaware",        "10", "Florida",         "12", "Georgia",        "13",
  "Idaho",           "16", "Illinois",        "17", "Indiana",        "18",
  "Iowa",            "19", "Kansas",          "20", "Kentucky",       "21",
  "Louisiana",       "22", "Maine",           "23", "Maryland",       "24",
  "Massachusetts",   "25", "Michigan",        "26", "Minnesota",      "27",
  "Mississippi",     "28", "Missouri",        "29", "Montana",        "30",
  "Nebraska",        "31", "Nevada",          "32", "New Hampshire",  "33",
  "New Jersey",      "34", "New Mexico",      "35", "New York",       "36",
  "North Carolina",  "37", "North Dakota",    "38", "Ohio",           "39",
  "Oklahoma",        "40", "Oregon",          "41", "Pennsylvania",   "42",
  "Rhode Island",    "44", "South Carolina",  "45", "South Dakota",   "46",
  "Tennessee",       "47", "Texas",           "48", "Utah",           "49",
  "Vermont",         "50", "Virginia",        "51", "Washington",     "53",
  "West Virginia",   "54", "Wisconsin",       "55", "Wyoming",        "56"
)

noaa9_lookup <- noaa9_states %>% left_join(state_fips, by = "state_name")
stopifnot(!any(is.na(noaa9_lookup$STATEFP)))  # every NOAA9 state name should match a known FIPS

noaa9_zones <- counties_sf %>%
  st_drop_geometry() %>%
  transmute(StCoFIPS = GEOID, STATEFP) %>%
  inner_join(noaa9_lookup %>% dplyr::select(STATEFP, noaa9 = noaa9_region), by = "STATEFP") %>%
  dplyr::select(StCoFIPS, noaa9)

# --- common county set -------------------------------------------------------
# Inner join across all 3 methods + 4 references so V-measure and the CH
# comparison run on an IDENTICAL county set - ClustGeo already drops a handful
# of counties with no matching geometry (section 3.1), so without this
# common-set step the group counts would not be directly comparable.
common <- cluster_kmeans %>%
  dplyr::select(StCoFIPS, cluster_km = cluster) %>%
  inner_join(cluster_clustgeo %>% dplyr::select(StCoFIPS, cluster_geo = cluster), by = "StCoFIPS") %>%
  inner_join(cluster_hclust %>% dplyr::select(StCoFIPS, cluster_hc = cluster), by = "StCoFIPS") %>%
  inner_join(climate_zones, by = "StCoFIPS") %>%
  inner_join(koppen_zones,  by = "StCoFIPS") %>%
  inner_join(ipcc_zones,    by = "StCoFIPS") %>%
  inner_join(noaa9_zones,   by = "StCoFIPS") %>%
  inner_join(county_heat_characteristics %>% dplyr::select(StCoFIPS, dplyr::all_of(cluster_vars)),
            by = "StCoFIPS") %>%
  mutate(iecc21 = factor(iecc21), koppen = factor(koppen), ipcc = factor(ipcc), noaa9 = factor(noaa9))

cat("\ncommon county set (k-means, ClustGeo, hierarchical, IECC21, Koppen, IPCC, NOAA9, all non-missing):",
    nrow(common), "\n")
cat("IECC21 zones:", nlevels(common$iecc21), "| Koppen classes:", nlevels(common$koppen),
    "| IPCC regions:", nlevels(common$ipcc), "| NOAA9 regions:", nlevels(common$noaa9), "\n")

##---- 5.2 V-measure: each clustering vs each reference #----
method_cols   <- c("k-means" = "cluster_km", "ClustGeo" = "cluster_geo", "Hierarchical" = "cluster_hc")
external_cols <- c("IECC21" = "iecc21", "Koppen-Geiger" = "koppen",
                   "IPCC WGI region" = "ipcc", "NOAA 9 climate region" = "noaa9")

vmeasure_grid <- expand.grid(method = names(method_cols), external = names(external_cols),
                             stringsAsFactors = FALSE)

vmeasure_tab <- vmeasure_grid %>%
  dplyr::mutate(
    n            = nrow(common),
    homogeneity  = mapply(function(m, e) homogeneity(common[[external_cols[e]]], common[[method_cols[m]]]),
                          method, external),
    completeness = mapply(function(m, e) completeness(common[[external_cols[e]]], common[[method_cols[m]]]),
                          method, external),
    v_measure    = mapply(function(m, e) v_measure(common[[external_cols[e]]], common[[method_cols[m]]]),
                          method, external)
  ) %>%
  dplyr::arrange(external, method) %>%
  dplyr::transmute(Method = method, `External system` = external, N = n,
                   Homogeneity = round(homogeneity, 4), Completeness = round(completeness, 4),
                   `V-measure` = round(v_measure, 4))

cat("\n=== V-measure vs 4 external systems (0 = no agreement, 1 = identical partitions) ===\n")
print(as.data.frame(vmeasure_tab), row.names = FALSE)

##---- 5.3 CH index: all 7 labellings on the SAME feature space #----
# Same scaling convention as section 2.1 (scale() on the portfolio columns),
# recomputed on the common county set so all 7 CH values share one feature
# space. Uses the shared ch_index() helper from section 1 (fpc::calinhara()
# under the hood).
common_scaled <- scale(common[, cluster_vars])

group_cols <- c(method_cols, external_cols)

ch_tab <- tibble::tibble(
  Method   = names(group_cols),
  n_groups = sapply(group_cols, function(col) nlevels(droplevels(as.factor(common[[col]])))),
  ch_index = sapply(group_cols, function(col) ch_index(common_scaled, common[[col]]))
) %>%
  dplyr::mutate(ch_index = round(ch_index, 1)) %>%
  dplyr::rename(`CH index` = ch_index, `N groups` = n_groups)

cat("\n=== CH index on the heatwave-metric feature space (higher = better separated) ===\n")
print(as.data.frame(ch_tab), row.names = FALSE)

write_csv(vmeasure_tab, file.path(season_dir, "external_vmeasure_comparison.csv"))
write_csv(ch_tab,       file.path(season_dir, "external_ch_comparison.csv"))

##---- 5.4 export to Word #----
doc <- officer::read_docx()
doc <- officer::body_add_par(doc, "Table 1. V-measure - each clustering method vs each external system",
                             style = "heading 1")
doc <- officer::body_add_flextable(doc, flextable::autofit(flextable::flextable(vmeasure_tab)))
doc <- officer::body_add_par(doc, "Table 2. CH index - all methods and external systems, same feature space",
                             style = "heading 1")
doc <- officer::body_add_flextable(doc, flextable::autofit(flextable::flextable(ch_tab)))

docx_path <- file.path(fig_dir, "external_climate_comparison.docx")
print(doc, target = docx_path)

message("done - exported: external_vmeasure_comparison.csv, external_ch_comparison.csv, ",
       basename(docx_path))


#===========================================================================
# 6. CLUSTER DESCRIPTIVES: COUNTIES, POPULATION, % US POPULATION
#===========================================================================
# For the primary (k-means) 7 mild -> severe exposure clusters: number of
# counties, total population, and % of total US population. Population is
# tract-level (tract_master_2020.rds's pop column) summed up to county via
# StCoFIPS = substr(GEOID, 1, 5) - cached the same way (tract_pop.csv)
# domain_exposure_heatmap.R already does, so this reuses that cache if it
# already exists rather than re-reading the 542 MB source file.

pop_cache <- file.path(season_dir, "tract_pop.csv")
if (!file.exists(pop_cache)) {
  message("caching tract population from tract_master_2020.rds ...")
  tm <- readRDS(file.path(data_dir, "tract_master_2020.rds"))
  tm %>% sf::st_drop_geometry() %>%
    dplyr::select(GEOID, pop) %>%
    write_csv(pop_cache)
  rm(tm); invisible(gc())
}

tract_pop <- read_csv(pop_cache, col_types = cols(GEOID = col_character(), pop = col_double()))

county_pop <- tract_pop %>%
  dplyr::mutate(StCoFIPS = substr(GEOID, 1, 5)) %>%
  dplyr::group_by(StCoFIPS) %>%
  dplyr::summarise(pop = sum(pop, na.rm = TRUE), .groups = "drop")

cluster_describe <- cluster_kmeans %>%
  dplyr::select(StCoFIPS, cluster) %>%
  dplyr::left_join(county_pop, by = "StCoFIPS")

n_missing_pop <- sum(is.na(cluster_describe$pop))
if (n_missing_pop > 0) {
  message(n_missing_pop, " counties have no population match - excluded from the population totals below")
}

total_us_pop <- sum(cluster_describe$pop, na.rm = TRUE)

cluster_summary <- cluster_describe %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(n_counties = dplyr::n(), population = sum(pop, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(pct_us_population = round(100 * population / total_us_pop, 2)) %>%
  dplyr::arrange(cluster)

cat("\n=== cluster descriptives: counties, population, % US population ===\n")
print(as.data.frame(cluster_summary), row.names = FALSE)

write_csv(cluster_summary, file.path(season_dir, "cluster_descriptives.csv"))


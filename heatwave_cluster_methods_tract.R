#---- Script Metadata #----
# Title: Heatwave clustering method comparison - TRACT level
# Author: Hang Li
# Purpose: Same comparison as heatwave_cluster_methods.R (k-means vs Ward
# hierarchical, vs IECC 2021 climate zones), but on TRACT-level heat
# characteristics instead of county-level, and simplified to a diagnostic-only
# check: no ClustGeo (no spatial constraint, no tract geometry needed at all),
# no file exports. Just: profile (elbow/silhouette/CH) -> pick k -> fit ->
# print silhouette, for both methods, then compare method vs method and
# method vs IECC21.
#
# Input: heat_characteristics_tract_2001_2020.csv - GEOID-keyed, 14 metrics
# (7 per {rel, abs} heatwave definition). Column names differ slightly from
# the county-level file: "rel_intensity_max"/"abs_intensity_max" here, not
# "..._tmax" - substituted below (flagged, not silent) since cluster_vars as
# given named the county file's column, which does not exist in this one.
#
# SCALE WARNING: this file has ~149,000 tracts, not ~3,100 counties.
#   - stats::hclust()/dist() need O(n^2) memory - infeasible here (would need
#     tens of TB). fastcluster::hclust.vector(method="ward") computes the
#     identical Ward tree (verified: merge matrix, heights, and cutree() all
#     matched stats::hclust(dist(x), "ward.D2") exactly on a test dataset) via
#     the O(n) memory NN-chain algorithm instead - the only reason hierarchical
#     clustering is feasible at all at this n. One fit still took ~5 minutes
#     on this machine for n = 149,037 - expect the hierarchical section to be
#     the slow part of this script.
#   - cluster::silhouette() also needs an O(n^2) pairwise-distance object, so
#     it is computed on a random SUBSAMPLE of tracts (SIL_N below), not all
#     149k - standard practice for large-n silhouette estimation. k-means
#     fitting, CH index, and the elbow diagnostics do NOT need pairwise
#     distances and run on the FULL tract set.
#---------------------------------------------------------------------------

#===========================================================================
# 1. SETUP
#===========================================================================
rm(list = ls())
pkgs <- c("dplyr", "readr", "tibble", "ggplot2", "cluster", "fastcluster", "fpc", "sf", "clevr")
install.packages(setdiff(pkgs, rownames(installed.packages())))
invisible(lapply(pkgs, library, character.only = TRUE))

data_dir   <- "/Users/hangli/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Aune, Kyle T's files - 260204 - US Heat Clusters/Data"
season_dir <- file.path(data_dir, "processed", "season cluster")

tract_heat <- read_csv(
  file.path(season_dir, "heat_characteristics_tract_2001_2020.csv"),
  col_types = cols(GEOID = col_character(), .default = col_double())
)

cat("tract_heat:", nrow(tract_heat), "tracts x", ncol(tract_heat) - 1, "metrics\n")
cat("available metrics:\n")
print(setdiff(names(tract_heat), "GEOID"))

# CH (Calinski-Harabasz) index, backed by fpc::calinhara() rather than a
# hand-rolled formula. fpc::calinhara(x, clustering, cn = max(clustering)) is
# only correct when group labels are consecutive integers 1..k - confirmed
# separately that non-consecutive labels silently give the WRONG value with
# the default cn. as.integer(factor(...)) + explicit cn = nlevels(...)
# guarantees consecutive labels regardless of the input's original coding.
ch_index <- function(x, group) {
  group <- droplevels(as.factor(group))
  fpc::calinhara(x, as.integer(group), cn = nlevels(group))
}

#===========================================================================
# 2. VARIABLE PORTFOLIO AND SCALING
#===========================================================================
# ---------------------------------------------------------------- EDIT THIS --
cluster_vars <- c(
  "abs_hwdays_per_yr",
  "rel_intensity_max",   # tract-file equivalent of the county file's rel_intensity_tmax
  "rel_duration_max",
  "rel_hw_span"
)
# -----------------------------------------------------------------------------

cluster_input <- tract_heat %>%
  dplyr::select(GEOID, dplyr::all_of(cluster_vars)) %>%
  na.omit()

cluster_scaled <- scale(cluster_input[-1])
n_obs <- nrow(cluster_input)

cat("\nportfolio:", paste(cluster_vars, collapse = ", "), "\n")
cat("tracts:", n_obs, "of", nrow(tract_heat),
    sprintf("(%d dropped for NA)\n", nrow(tract_heat) - n_obs))

# Subsample used for every silhouette computation below (k-means and
# hierarchical alike), so the two methods' silhouette numbers stay directly
# comparable (same tracts evaluated both times).
SIL_N <- 5000
set.seed(123)
sil_idx <- sample.int(n_obs, min(SIL_N, n_obs))
d_sil   <- dist(cluster_scaled[sil_idx, ])
cat("silhouette subsample size:", length(sil_idx), "of", n_obs, "tracts\n")

K_RANGE <- 2:10

#===========================================================================
# 3. K-MEANS
#===========================================================================

##---- 3.1 profile: elbow, silhouette, CH index #----
profile_km <- t(sapply(K_RANGE, function(k) {
  km_k <- kmeans(cluster_scaled, centers = k, nstart = 25)
  c(wss        = km_k$tot.withinss,
    silhouette = mean(silhouette(km_k$cluster[sil_idx], d_sil)[, 3]),
    ch         = ch_index(cluster_scaled, km_k$cluster))
}))
profile_km <- data.frame(k = K_RANGE, profile_km)

cat("\nk-means profile (", length(cluster_vars), "-variable portfolio, tract level):\n", sep = "")
print(round(profile_km, 3), row.names = FALSE)
cat("silhouette peaks at k =", K_RANGE[which.max(profile_km$silhouette)],
    "| CH peaks at k =", K_RANGE[which.max(profile_km$ch)], "\n")

par(mfrow = c(1, 3))
plot(K_RANGE, profile_km$wss, type = "b", pch = 19,
     xlab = "k", ylab = "Total within-cluster SS", main = "1. Elbow")
plot(K_RANGE, profile_km$silhouette, type = "b", pch = 19,
     xlab = "k", ylab = "Mean silhouette width (subsample)", main = "2. Silhouette")
plot(K_RANGE, profile_km$ch, type = "b", pch = 19,
     xlab = "k", ylab = "Calinski-Harabasz index", main = "3. CH index")
par(mfrow = c(1, 1))

##---- 3.2 select k and fit #----
# ---------------------------------------------------------------- EDIT THIS --
optimal_k <- 7
# -----------------------------------------------------------------------------

set.seed(123)
km <- kmeans(cluster_scaled, centers = optimal_k, nstart = 50)

cat("\nk-means, k =", optimal_k, "- cluster sizes:\n")
print(table(km$cluster))
cat("k-means, k =", optimal_k, "- mean silhouette width (subsample):",
    round(mean(silhouette(km$cluster[sil_idx], d_sil)[, 3]), 4), "\n")
cat("k-means, k =", optimal_k, "- CH index (full data):",
    round(ch_index(cluster_scaled, km$cluster), 1), "\n")

#===========================================================================
# 4. WARD HIERARCHICAL CLUSTERING (no spatial constraint)
#===========================================================================
# fastcluster::hclust.vector(method = "ward") - see the scale warning in the
# metadata header for why stats::hclust(dist(...)) cannot be used at this n.
# Fit ONCE; every k below is a cutree() on this one tree, which is cheap.

cat("\nfitting Ward hierarchical tree on", n_obs, "tracts (fastcluster::hclust.vector) ...\n")
t0 <- Sys.time()
tree_hc <- fastcluster::hclust.vector(cluster_scaled, method = "ward")
cat("done in", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "sec\n")

##---- 4.1 profile: fusion-height elbow, silhouette, CH index #----
n_show  <- 15
heights <- rev(tail(tree_hc$height, n_show))

sil_profile_hc <- sapply(K_RANGE, function(k) {
  mean(silhouette(cutree(tree_hc, k = k)[sil_idx], d_sil)[, 3])
})
ch_profile_hc <- sapply(K_RANGE, function(k) {
  ch_index(cluster_scaled, cutree(tree_hc, k = k))
})

profile_hc <- data.frame(k = K_RANGE, silhouette = round(sil_profile_hc, 4),
                         ch = round(ch_profile_hc, 1))
cat("\nhierarchical (Ward) profile:\n")
print(profile_hc, row.names = FALSE)
cat("silhouette peaks at k =", K_RANGE[which.max(sil_profile_hc)],
    "| CH peaks at k =", K_RANGE[which.max(ch_profile_hc)], "\n")

par(mfrow = c(1, 3))
plot(seq_len(n_show), heights, type = "b", pch = 19,
     xlab = "merges before the final partition (1 = last)", ylab = "fusion height",
     main = "1. Elbow (fusion height)")
plot(K_RANGE, sil_profile_hc, type = "b", pch = 19,
     xlab = "k", ylab = "Mean silhouette width (subsample)", main = "2. Silhouette")
plot(K_RANGE, ch_profile_hc, type = "b", pch = 19,
     xlab = "k", ylab = "Calinski-Harabasz index", main = "3. CH index")
par(mfrow = c(1, 1))

##---- 4.2 select k and cut #----
# ---------------------------------------------------------------- EDIT THIS --
optimal_k_hc <- 7
# -----------------------------------------------------------------------------

raw_cut_hc <- cutree(tree_hc, k = optimal_k_hc)

cat("\nhierarchical (Ward), k =", optimal_k_hc, "- cluster sizes:\n")
print(table(raw_cut_hc))
cat("hierarchical (Ward), k =", optimal_k_hc, "- mean silhouette width (subsample):",
    round(mean(silhouette(raw_cut_hc[sil_idx], d_sil)[, 3]), 4), "\n")
cat("hierarchical (Ward), k =", optimal_k_hc, "- CH index (full data):",
    round(ch_index(cluster_scaled, raw_cut_hc), 1), "\n")

#===========================================================================
# 5. K-MEANS VS HIERARCHICAL - DIAGNOSTIC COMPARISON
#===========================================================================
diag_compare <- tibble::tibble(
  method     = c("k-means", "Hierarchical (Ward)"),
  k          = c(optimal_k, optimal_k_hc),
  silhouette = c(mean(silhouette(km$cluster[sil_idx], d_sil)[, 3]),
                mean(silhouette(raw_cut_hc[sil_idx], d_sil)[, 3])),
  ch_index   = c(ch_index(cluster_scaled, km$cluster),
                ch_index(cluster_scaled, raw_cut_hc))
)
cat("\n=== k-means vs hierarchical (Ward), diagnostics at each method's chosen k ===\n")
print(as.data.frame(diag_compare %>%
                      dplyr::mutate(silhouette = round(silhouette, 4),
                                   ch_index   = round(ch_index, 1))),
      row.names = FALSE)

# Cross-tab: do the two methods agree on which tracts go together, regardless
# of label numbering?
cat("\ncross-tab, k-means cluster (rows) x hierarchical cluster (cols):\n")
print(table(km$cluster, raw_cut_hc))

#===========================================================================
# 6. COMPARISON WITH IECC 2021 CLIMATE ZONES
#===========================================================================
# IECC21 is COUNTY-level; tracts are mapped to their county via the first 5
# digits of GEOID (StCoFIPS), same convention used elsewhere in this project
# (e.g. domain_exposure_heatmap.R). Every tract in a county inherits that
# county's IECC21 zone.

climate_zones <- st_read(file.path(data_dir, "ClimateZoneDataFiles/ClimateZones.shp"),
                         quiet = TRUE) %>%
  st_drop_geometry() %>%
  transmute(StCoFIPS = sub("^G", "", GEOID), iecc21 = IECC21) %>%
  filter(!is.na(iecc21))

common <- cluster_input %>%
  dplyr::select(GEOID) %>%
  dplyr::mutate(
    StCoFIPS   = substr(GEOID, 1, 5),
    cluster_km = km$cluster,
    cluster_hc = raw_cut_hc
  ) %>%
  dplyr::inner_join(climate_zones, by = "StCoFIPS") %>%
  dplyr::mutate(iecc21 = factor(iecc21))

cat("\ncommon tract set (k-means, hierarchical, IECC21, all non-missing):",
    nrow(common), "\n")
cat("IECC21 zones present:", nlevels(common$iecc21), "-",
    paste(levels(common$iecc21), collapse = ", "), "\n")

##---- 6.1 V-measure #----
vmeasure_tab <- tibble::tibble(
  comparison   = c("k-means vs IECC21", "Hierarchical vs IECC21"),
  n            = nrow(common),
  homogeneity  = c(homogeneity(common$iecc21, common$cluster_km),
                   homogeneity(common$iecc21, common$cluster_hc)),
  completeness = c(completeness(common$iecc21, common$cluster_km),
                   completeness(common$iecc21, common$cluster_hc)),
  v_measure    = c(v_measure(common$iecc21, common$cluster_km),
                   v_measure(common$iecc21, common$cluster_hc))
)

cat("\n=== V-measure vs IECC21 (0 = no agreement, 1 = identical partitions) ===\n")
print(as.data.frame(vmeasure_tab %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 4)))),
      row.names = FALSE)

##---- 6.2 CH index on the same feature space #----
common_scaled <- scale(tract_heat[match(common$GEOID, tract_heat$GEOID), cluster_vars])

ch_tab <- tibble::tibble(
  method   = c("k-means", "Hierarchical (Ward)", "IECC21"),
  n_groups = c(nlevels(droplevels(factor(common$cluster_km))),
              nlevels(droplevels(factor(common$cluster_hc))),
              nlevels(common$iecc21)),
  ch_index = c(ch_index(common_scaled, common$cluster_km),
              ch_index(common_scaled, common$cluster_hc),
              ch_index(common_scaled, common$iecc21))
)

cat("\n=== CH index on the heatwave-metric feature space (higher = better separated) ===\n")
print(as.data.frame(ch_tab %>% dplyr::mutate(ch_index = round(ch_index, 1))), row.names = FALSE)

message("done - diagnostic check only, nothing written to disk")

#---- Script Metadata #----
# Title: heat cluster
# Author: Hang Li
# Date: 07/08/2026


#---- Setup #----
rm(list = ls())
pkgs <- c("tidyverse", "sf", "stars", "raster", "tidycensus", "tigris", "exactextractr",
          "parallel", "doParallel", "foreach", "readxl", "httr", "jsonlite", "ggplot2")
lapply(pkgs, library, character.only = TRUE, quietly = TRUE, verbose = FALSE)

setwd("/Users/hangli/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Aune, Kyle T's files - 260204 - US Heat Clusters/Data")

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
wbgt_clim <- wbgt_clim %>%
  filter(!is.na(WBGT))

##---- elbow method to determine clusters#----
max_k <- 8

elbow_all <- wbgt_clim %>%
  group_by(StCoFIPS) %>%
  group_modify(~{
    
    wss <- sapply(1:max_k, function(k){
      
      kmeans(
        .x$WBGT,
        centers = k,
        nstart = 50
      )$tot.withinss
      
    })
    
    tibble(
      k = 1:max_k,
      wss = wss
    )
    
  }) %>%
  ungroup()
elbow_summary <- elbow_all %>%
  group_by(k) %>%
  summarise(
    mean_wss = mean(wss),
    median_wss = median(wss),
    sd_wss = sd(wss),
    se_wss = sd_wss / sqrt(n()),
    .groups = "drop"
  )
elbow <- ggplot(elbow_summary, aes(k, mean_wss)) +
  geom_line(size = 0.5) +
  geom_point(size = 1) +
  geom_errorbar(
    aes(
      ymin = mean_wss - se_wss,
      ymax = mean_wss + se_wss
    ),
    width = 0.15
  ) +
  scale_x_continuous(breaks = 1:max_k) +
  labs(
    x = "Number of clusters (k)",
    y = "Mean within-cluster sum of squares",
    title = "Elbow method across all U.S. counties"
  ) +
  theme_bw(base_size = 12)
ggsave(
  filename = "../Output/Figures/elbow_plot.png",
  plot = elbow,
  width = 7,
  height = 5,
  dpi = 600
)
##---- k-means for cold temp warm#----
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

write_csv(monthly_season, "/Users/hangli/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Paper/US_heat_cluster/us-heat-clusters/data/processed/season cluster/monthly_season.csv")


#---- County heatwaves#----
# Two heatwave definitions, both requiring at least MIN_SPELL consecutive days
# inside the county's own warm season (as classified in the k-means step above):
#
#   A. RELATIVE  WBGTmean_C >= that county's 95th percentile of warm-season
#                WBGTmean_C. Scales to local acclimatisation.
#   B. ABSOLUTE  WBGTmax_C >= 31 C. A fixed physiological threshold, so a
#                northern county can have zero events.
#
# Everything below is restricted to 2000-2020.

STUDY_START    <- as.Date("2000-01-01")
STUDY_END      <- as.Date("2020-12-31")
MIN_SPELL      <- 2L      # consecutive days required to count as a heatwave
ABS_WBGT       <- 31      # degrees C, absolute-threshold definition - drives
                          # every abs_* metric except the sensitivity columns below
ABS_WBGT_EXTRA <- c(28, 29, 30, 32)  # degrees C - abs_hwdays_per_yr ONLY (see
                          # Heatwave characteristics section); every other
                          # abs_* metric stays at the single ABS_WBGT = 31 C
                          # definition

heat_warm <- heat2 %>%
  left_join(
    monthly_season %>% dplyr::select(StCoFIPS, month, season),
    by = c("StCoFIPS", "month")
  ) %>%
  filter(season == "Warm", Date >= STUDY_START, Date <= STUDY_END) %>%
  mutate(year = year(Date))

# County-specific relative threshold, computed on warm-season days only
heat_thresh <- heat_warm %>%
  group_by(StCoFIPS) %>%
  summarise(WBGTmean_95th = quantile(WBGTmean_C, 0.95, na.rm = TRUE),
            .groups = "drop")

# Break a county's series into spells that share a flag value AND are contiguous
# in calendar time. The warm-season filter removes the cold months, so
# consecutive ROWS are not necessarily consecutive DAYS - without the date test
# the last hot day of one summer and the first of the next would be merged into
# a single spell, inflating both frequency and duration.
spell_id <- function(date, flag) {
  n <- length(flag)
  if (n == 0L) return(integer(0))
  if (n == 1L) return(1L)
  brk <- c(TRUE,
           (flag[-1] != flag[-n]) |
             (as.integer(date[-1]) != as.integer(date[-n]) + 1L))
  cumsum(brk)
}

# TRUE on days belonging to a run of >= min_len consecutive flagged days
mark_hw <- function(date, flag, min_len = MIN_SPELL) {
  flag <- !is.na(flag) & flag
  g    <- spell_id(date, flag)
  len  <- ave(g, g, FUN = length)
  flag & len >= min_len
}

# hw_abs_28/29/30/32 are extra absolute-threshold flags used ONLY to build the
# abs_hwdays_per_yr sensitivity columns below - every other abs_* metric comes
# from hw_abs (ABS_WBGT = 31) alone, unchanged. Written out explicitly (not
# built from ABS_WBGT_EXTRA programmatically) because mutate()'s data-masking
# evaluates column expressions lazily inside the pipeline - a list of
# mark_hw() calls built by lapply() BEFORE entering mutate() would try to
# resolve Date/WBGTmax_C in the calling environment instead, where they don't
# exist (confirmed: that approach fails with "object 'WBGTmax_C' not found").
heat_hw <- heat_warm %>%
  left_join(heat_thresh, by = "StCoFIPS") %>%
  arrange(StCoFIPS, Date) %>%
  group_by(StCoFIPS) %>%
  mutate(
    hw_rel    = mark_hw(Date, WBGTmean_C >= WBGTmean_95th),
    hw_abs    = mark_hw(Date, WBGTmax_C  >= ABS_WBGT),
    hw_abs_28 = mark_hw(Date, WBGTmax_C  >= 28),
    hw_abs_29 = mark_hw(Date, WBGTmax_C  >= 29),
    hw_abs_30 = mark_hw(Date, WBGTmax_C  >= 30),
    hw_abs_32 = mark_hw(Date, WBGTmax_C  >= 32)
  ) %>%
  ungroup() %>%
  dplyr::select(
    StCoFIPS, Date, year, doy, month, season,
    WBGTmin_C, WBGTmean_C, WBGTmax_C, WBGTmean_95th,
    hw_rel, hw_abs, dplyr::all_of(paste0("hw_abs_", ABS_WBGT_EXTRA))
  )

cat("warm-season days:", format(nrow(heat_hw), big.mark = ","),
    "| counties:", n_distinct(heat_hw$StCoFIPS),
    "| years:", n_distinct(heat_hw$year), "\n")
cat("heatwave days - relative:", format(sum(heat_hw$hw_rel), big.mark = ","),
    " absolute:", format(sum(heat_hw$hw_abs), big.mark = ","), "\n")

write_csv(heat_hw, "processed/season cluster/heat_wave.csv")


#---- Heatwave characteristics#----
# Eight metrics per definition, in four families:
#   intensity  mean WBGTmax and mean WBGTmin across heatwave days
#   frequency  events per year and heatwave days per year
#   duration   mean event length and longest event length
#   timing     warm-season length, and first-to-last heatwave span within a year
#
# Rates are per year over the observed record rather than totals, so counties
# are comparable even if a few have incomplete years.
#
# Counties with zero events under a definition get 0 for the count-based metrics
# and NA for the ones that are undefined without an event (intensity, duration,
# span). That distinction matters for the absolute threshold, where northern
# counties legitimately never reach 31 C.

hw_summary <- function(dat, hw_col, prefix) {

  hwrows <- dat %>% dplyr::filter(.data[[hw_col]])

  base <- dat %>%
    group_by(StCoFIPS) %>%
    summarise(n_years   = n_distinct(year),
              warm_days = n(),
              .groups   = "drop")

  inten <- hwrows %>%
    group_by(StCoFIPS) %>%
    summarise(intensity_tmax = mean(WBGTmax_C, na.rm = TRUE),
              intensity_tmin = mean(WBGTmin_C, na.rm = TRUE),
              hw_days        = n(),
              .groups        = "drop")

  # Spell ids computed on heatwave days only: the date test still separates
  # events, since non-heatwave days leave a calendar gap.
  ev <- hwrows %>%
    group_by(StCoFIPS) %>%
    mutate(ev = spell_id(Date, rep(TRUE, n()))) %>%
    group_by(StCoFIPS, ev) %>%
    summarise(len = n(), .groups = "drop") %>%
    group_by(StCoFIPS) %>%
    summarise(n_events      = n(),
              duration_mean = mean(len),
              duration_max  = max(len),
              .groups       = "drop")

  # Within-year span from the first to the last heatwave day, averaged over
  # years that had at least one event. CONUS warm seasons sit inside a calendar
  # year, so day-of-year differences do not wrap.
  span <- hwrows %>%
    group_by(StCoFIPS, year) %>%
    summarise(sp = max(doy) - min(doy) + 1, .groups = "drop") %>%
    group_by(StCoFIPS) %>%
    summarise(hw_span = mean(sp), .groups = "drop")

  out <- base %>%
    left_join(inten, by = "StCoFIPS") %>%
    left_join(ev,    by = "StCoFIPS") %>%
    left_join(span,  by = "StCoFIPS") %>%
    transmute(
      StCoFIPS,
      intensity_tmax,
      intensity_tmin,
      events_per_yr = coalesce(n_events, 0L) / n_years,
      hwdays_per_yr = coalesce(hw_days,  0L) / n_years,
      duration_mean,
      duration_max,
      season_days   = warm_days / n_years,
      hw_span
    )

  names(out)[-1] <- paste0(prefix, "_", names(out)[-1])
  out
}

county_heat_characteristics <-
  hw_summary(heat_hw, "hw_rel", "rel") %>%
  full_join(hw_summary(heat_hw, "hw_abs", "abs"), by = "StCoFIPS")

# season_days is a property of the county, not of the definition, so the two
# copies are identical by construction - kept for symmetry between the sets.
stopifnot(all.equal(county_heat_characteristics$rel_season_days,
                    county_heat_characteristics$abs_season_days))

# abs_hwdays_per_yr sensitivity: same metric (heatwave days per year, still
# requiring MIN_SPELL consecutive days), recomputed at every threshold in
# c(ABS_WBGT, ABS_WBGT_EXTRA) = 28-32 C. Every OTHER abs_* metric (intensity,
# events, duration, span) is left at the single ABS_WBGT = 31 C definition
# above - recomputing those per threshold isn't needed for a frequency
# sensitivity check, so only hwdays_per_yr is repeated, via a stripped-down
# version of hw_summary() rather than the full function.
hwdays_only <- function(dat, hw_col) {
  base <- dat %>%
    group_by(StCoFIPS) %>%
    summarise(n_years = n_distinct(year), .groups = "drop")
  hwd <- dat %>%
    dplyr::filter(.data[[hw_col]]) %>%
    group_by(StCoFIPS) %>%
    summarise(hw_days = n(), .groups = "drop")
  base %>%
    left_join(hwd, by = "StCoFIPS") %>%
    transmute(StCoFIPS, hwdays_per_yr = coalesce(hw_days, 0L) / n_years)
}

abs_wbgt_all <- sort(c(ABS_WBGT, ABS_WBGT_EXTRA))
abs_hwdays_multi <- Reduce(
  function(x, y) full_join(x, y, by = "StCoFIPS"),
  lapply(abs_wbgt_all, function(t) {
    hw_col <- if (t == ABS_WBGT) "hw_abs" else paste0("hw_abs_", t)
    hwdays_only(heat_hw, hw_col) %>%
      rename(!!paste0("abs_hwdays_per_yr_", t) := hwdays_per_yr)
  })
)

county_heat_characteristics <- county_heat_characteristics %>%
  full_join(abs_hwdays_multi, by = "StCoFIPS")

cat("\nabs_hwdays_per_yr across thresholds", paste(abs_wbgt_all, collapse = ", "), "C:\n")
print(summary(county_heat_characteristics[paste0("abs_hwdays_per_yr_", abs_wbgt_all)]))

cat("\ncounty_heat_characteristics:", nrow(county_heat_characteristics),
    "counties x", ncol(county_heat_characteristics) - 1, "metrics\n")
cat("counties with no ABSOLUTE-threshold heatwave:",
    sum(county_heat_characteristics$abs_events_per_yr == 0), "\n\n")
print(summary(county_heat_characteristics[-1]))

write_csv(county_heat_characteristics,
          "processed/season cluster/heat_characteristics_v2.csv")

#---- Cluster #----
# Workflow for choosing a clustering:
#   1  pick a variable portfolio          <- EDIT cluster_vars
#   2  build and scale the matrix
#   3  kfx(): NbClust silhouette vote over repeated subsamples
#   4  silhouette profile on the FULL data, for comparison with kfx
#   5  fit the final k-means              <- EDIT optimal_k
#   5b relabel clusters mild -> severe (moved out of the mapping section, so
#      the exported county_clustered.csv already carries the final IDs)
#   6  silhouette diagnostics for the chosen solution
#
# Steps 1-4 are cheap to re-run. Change cluster_vars, source from here down,
# and compare the kfx vote and the silhouette profile before committing to a k.

library(cluster)
library(NbClust)

##---- 1. variable portfolio #----

# What is available (8 metrics x 2 heatwave definitions):
#
#   intensity   *_intensity_tmax   *_intensity_tmin
#   frequency   *_events_per_yr    *_hwdays_per_yr
#   duration    *_duration_mean    *_duration_max
#   timing      *_season_days      *_hw_span
#
# prefix rel_ = relative (county 95th percentile), abs_ = absolute (WBGT >= 31).
cat("\navailable metrics:\n")
print(setdiff(names(county_heat_characteristics), "StCoFIPS"))

cat("\ncompleteness (counties with a non-NA value):\n")
print(colSums(!is.na(county_heat_characteristics[-1])))

# ---------------------------------------------------------------- EDIT THIS --
cluster_vars <- c(
  "abs_hwdays_per_yr_29",
  "rel_intensity_tmax",
  "rel_duration_max",
  "rel_hw_span"
)
# -----------------------------------------------------------------------------

# Short axis/legend labels derived from cluster_vars - reused by the severity
# ranking below (step 5b) and by the radar chart in the mapping section, so
# there is one source of truth instead of a separately maintained label list
# that can drift out of step when the portfolio changes.
radar_vars <- setNames(
  sub("^(rel|abs)_", "", cluster_vars) |>
    sub(pattern = "intensity_tmax", replacement = "Intensity") |>
    sub(pattern = "hwdays_per_yr",    replacement = "Frequency") |>
    sub(pattern = "duration_max",  replacement = "Duration") |>
    sub(pattern = "season_days", replacement = "season") |>
    sub(pattern = "hw_span",     replacement = "Season length"),
  cluster_vars
)

##---- 2. build and scale #----
# Scaling is not optional here: season_days is ~150 and events_per_yr is ~1.5,
# so on raw values Euclidean distance would be almost entirely season_days.
cluster_input <- county_heat_characteristics %>%
  dplyr::select(StCoFIPS, dplyr::all_of(cluster_vars)) %>%
  na.omit()

cluster_scaled <- scale(cluster_input[-1])

cat("\nportfolio:", paste(cluster_vars, collapse = ", "), "\n")
cat("counties:", nrow(cluster_input), "of", nrow(county_heat_characteristics),
    sprintf("(%d dropped for NA)\n",
            nrow(county_heat_characteristics) - nrow(cluster_input)))

# Highly collinear inputs effectively double-weight whatever they measure.
cor_in <- cor(cluster_scaled)
cat("\ncorrelation among portfolio variables:\n")
print(round(cor_in, 2))
hi <- which(abs(cor_in) > 0.9 & upper.tri(cor_in), arr.ind = TRUE)
if (nrow(hi)) {
  cat("note - |r| > 0.9 between:",
      paste(sprintf("%s/%s", rownames(cor_in)[hi[, 1]], colnames(cor_in)[hi[, 2]]),
            collapse = "; "), "\n")
}

##---- 3. kfx: repeated-subsample NbClust vote #----
# Specify dataframe, variable name(s), # of replicants, and the k range to
# search. min.nc/max.nc generalize the original hardcoded 2-5 range: the
# indicator columns (k2, k3, ... kN) and the "best" lookup are built from
# min.nc:max.nc rather than a fixed set of four, so this works for any range.
kfx <- function(dat, var, rep, min.nc = 2, max.nc = 10) {
  require(dplyr)
  require(doParallel)
  require(foreach)
  require(parallel)

  # Setting up core cluster
  if (rep > (detectCores() - 1)) {
    n.cl <- detectCores() - 1}
  else {
    n.cl <- rep}

  cl <- makeCluster(n.cl)
  registerDoParallel(cl)

  k_seq  <- min.nc:max.nc
  k_cols <- paste0("k", k_seq)

  # Repeating test rep times
  set.seed(123)
  out.ls <- foreach(ii = seq_len(rep)) %dopar% {
    require(NbClust)

    dat.sub <- dat[sample(seq_len(nrow(dat)),
                          size = nrow(dat) / 5,
                          replace = FALSE),
                   var]

    tryCatch({
      x <- NbClust(
        data = dat.sub,
        distance = "euclidean",
        min.nc = min.nc,
        max.nc = max.nc,
        method = "kmeans",
        index = "silhouette")

      y <- x$Best.nc[1]

      # one 0/1 indicator per candidate k in min.nc:max.nc, built dynamically
      # rather than one hardcoded column per k value
      ind        <- as.list(as.integer(k_seq == y))
      names(ind) <- k_cols
      out        <- as.data.frame(c(list(rep = ii), ind))

      return(out)
      },
      error = function(e) {
        message(paste0("An error occurred on run ", ii))
        print(e)
        }
    )
    }
  stopCluster(cl)
  out <- bind_rows(out.ls) %>%
    mutate(best = k_cols[max.col(.[, k_cols], ties.method = "first")])

  return(out)
}

# Each replicate uses a fresh subsample, so the spread across replicates is the
# point: a portfolio with real structure votes consistently. (set.seed on the
# master does not control the workers, so the vote shifts slightly run to run;
# wrap in doRNG if you need it bit-reproducible.)
#
# min.nc/max.nc here match K_RANGE in step 4, so the subsample vote and the
# full-data silhouette profile are directly comparable over the same range.

kfx_res <- kfx(dat    = as.data.frame(cluster_scaled),
               var    = colnames(cluster_scaled),
               rep    = 20,
               min.nc = 2,
               max.nc = 10)

cat("\nkfx vote across replicates:\n")
print(table(kfx_res$best))
cat("\nper-replicate detail:\n")
print(kfx_res, row.names = FALSE)

##---- 4. silhouette profile on the full data #----
# kfx votes on subsamples; this is the same criterion on all counties, over a
# wider k range. Read the two together rather than trusting either alone.
K_RANGE <- 2:10

set.seed(123)
d_full <- dist(cluster_scaled)

sil_profile <- sapply(K_RANGE, function(k) {
  km_k <- kmeans(cluster_scaled, centers = k, nstart = 50)
  mean(silhouette(km_k$cluster, d_full)[, 3])
})

wss_profile <- sapply(K_RANGE, function(k) {
  kmeans(cluster_scaled, centers = k, nstart = 50)$tot.withinss
})

profile_tab <- data.frame(k = K_RANGE,
                          mean_silhouette = round(sil_profile, 4),
                          within_ss       = round(wss_profile, 1))
par(mfrow = c(1, 2))

plot(K_RANGE, wss_profile, type = "b", pch = 19,
     xlab = "k", ylab = "Total within-cluster SS",
     main = paste("Elbow -", length(cluster_vars), "variable portfolio"))

plot(K_RANGE, sil_profile, type = "b", pch = 19,
     xlab = "k", ylab = "Mean silhouette width",
     main = paste("Silhouette -", length(cluster_vars), "variable portfolio"))

par(mfrow = c(1, 1))

##---- 5. final k-means #----

optimal_k <- 7

set.seed(123)
km <- kmeans(cluster_scaled, centers = optimal_k, nstart = 100)

# km$cluster is a plain integer vector in the row order of cluster_scaled,
# which came from cluster_input - join back on StCoFIPS rather than assuming
# county_heat_characteristics has the same rows in the same order. Kept as an
# integer (not a factor) until the mild -> severe relabeling in 5b, so the
# intermediate joins do not have to reconcile factor levels.
county_clustered <- cluster_input %>%
  dplyr::mutate(cluster_raw = km$cluster) %>%
  dplyr::select(StCoFIPS, cluster_raw) %>%
  left_join(county_heat_characteristics, by = "StCoFIPS")

##---- 5b. relabel clusters mild -> severe #----
# k-means cluster IDs are arbitrary. Rank them by mean rel_intensity_tmax alone
# (ascending) and relabel 1 = coolest ... optimal_k = hottest. Ranking runs on
# the RAW partition, so it does not depend on the labels being replaced.
#
# Previously this ranked on a composite (sum of z-scored cluster_vars). Now it
# is single-variable, so scaling would not change the order - the mean is
# taken on the raw rel_intensity_tmax directly, no z-scoring needed.
#
# This used to happen later, in the mapping section, and only touched an
# in-memory column - county_clustered.csv was exported with the raw k-means
# IDs. It now happens here, before export, so `cluster` in the exported CSV
# (and everywhere downstream) is already the mild -> severe ID. The original
# k-means label survives as `cluster_raw`, kept for traceability.
cluster_tmax <- county_clustered %>%
  dplyr::group_by(cluster_raw) %>%
  dplyr::summarise(mean_tmax = mean(rel_intensity_tmax, na.rm = TRUE),
                   .groups = "drop")

cluster_key <- cluster_tmax %>%
  dplyr::arrange(mean_tmax) %>%
  dplyr::transmute(cluster_raw, cluster = dplyr::row_number(), mean_tmax)

cat("\ncrosswalk (raw k-means ID -> coolest..hottest ID by mean rel_intensity_tmax):\n")
print(as.data.frame(cluster_key), row.names = FALSE)
# write_csv(cluster_key, "processed/season cluster/cluster_severity_crosswalk.csv")

county_clustered <- county_clustered %>%
  dplyr::left_join(cluster_key %>% dplyr::select(cluster_raw, cluster),
                   by = "cluster_raw") %>%
  dplyr::mutate(cluster = factor(cluster, levels = seq_len(optimal_k))) %>%
  dplyr::relocate(cluster, cluster_raw, .after = StCoFIPS)

stopifnot(!anyNA(county_clustered$cluster))   # every county should have matched

##---- 6. silhouette diagnostics for the chosen solution #----
# silhouette() itself is computed on the RAW partition (km$cluster / d_full
# are in that order and know nothing about the relabeling) - only the printed
# summary is shown in the mild -> severe labelling, via cluster_key, so the
# diagnostics read consistently with everything else in the script.
sil <- silhouette(km$cluster, d_full)

cat("\n=== silhouette, k =", optimal_k, "===\n")
cat("mean silhouette width:", round(mean(sil[, 3]), 4), "\n")
cat("counties with negative silhouette (better fit in another cluster):",
    sum(sil[, 3] < 0), sprintf("(%.1f%%)\n", 100 * mean(sil[, 3] < 0)))

sil_by_cluster <- data.frame(
  cluster_raw = as.integer(sil[, 1]),
  width       = sil[, 3]
) %>%
  dplyr::left_join(cluster_key %>% dplyr::select(cluster_raw, cluster), by = "cluster_raw") %>%
  group_by(cluster) %>%
  summarise(n            = n(),
            mean_width   = round(mean(width), 4),
            pct_negative = round(100 * mean(width < 0), 1),
            .groups      = "drop") %>%
  arrange(cluster)

cat("\nby cluster (mild -> severe order):\n")
print(as.data.frame(sil_by_cluster), row.names = FALSE)
cat("\n(interpretation: >0.50 strong, 0.25-0.50 weak but usable,\n")
cat(" <0.25 means the structure is largely arbitrary)\n")

plot(sil, border = NA,
     main = paste("Silhouette, k =", optimal_k))

# Cluster centres in the ORIGINAL units of the portfolio variables, already in
# mild -> severe order since county_clustered$cluster is relabeled.
cluster_profile_raw <- county_clustered %>%
  group_by(cluster) %>%
  summarise(n_counties = n(),
            across(dplyr::all_of(cluster_vars), ~ round(mean(.x, na.rm = TRUE), 3)),
            .groups = "drop") %>%
  arrange(cluster)

cat("\ncluster centres (original units, mild -> severe order):\n")
print(as.data.frame(cluster_profile_raw), row.names = FALSE)

write_csv(county_clustered, "processed/season cluster/county_clustered.csv")
write_csv(profile_tab,      "processed/season cluster/cluster_k_profile.csv")
write_csv(kfx_res,          "processed/season cluster/cluster_kfx_votes.csv")


#---- mapping #----
# Both figures below key off county_clustered$cluster, which is already the
# mild -> severe ID from the Cluster section (1 = mildest ... optimal_k = most
# severe) - no re-ranking happens here.

library(dplyr)
library(sf)
library(tigris)
library(tmap)
library(tidyr)
library(patchwork)

options(tigris_use_cache = TRUE)

##---- 1. county geometry + join #----
counties_conus <- counties(
  year = 2020,
  cb = TRUE,
  class = "sf"
) %>%
  mutate(GEOID = as.character(GEOID)) %>%
  filter(
    !STATEFP %in% c(
      "02", # Alaska
      "15", # Hawaii
      "60", # American Samoa
      "66", # Guam
      "69", # Northern Mariana Islands
      "72", # Puerto Rico
      "78"  # U.S. Virgin Islands
    )
  )

# Defensive - StCoFIPS is already a 5-character string from the source data,
# but this is a no-op in that case and cheap insurance if it ever is not.
county_clustered <- county_clustered %>%
  mutate(StCoFIPS = sprintf("%05s", as.character(StCoFIPS)))

counties_clustered_sf <- counties_conus %>%
  left_join(
    county_clustered %>%
      dplyr::select(StCoFIPS, cluster, cluster_raw, dplyr::all_of(cluster_vars)),
    by = c("GEOID" = "StCoFIPS")
  ) %>%
  st_transform(5070)

match_check <- counties_clustered_sf %>%
  st_drop_geometry() %>%
  summarise(
    total_counties     = n(),
    matched_counties   = sum(!is.na(cluster)),
    unmatched_counties = sum(is.na(cluster))
  )
cat("\ncounty geometry match:\n")
print(as.data.frame(match_check), row.names = FALSE)

##---- 2. severity palette #----
# Sequential yellow -> orange -> red -> purple ramp (ColorBrewer's YlOrRd
# extended with a purple tail), in the spirit of the NWS HeatRisk scale, which
# uses yellow/orange/red/magenta for increasing heat danger. Kept off a
# blue-to-red diverging ramp on purpose: every cluster here is a flavor of
# "heatwave", not a mix of hot and cold categories, so blue at the mild end
# would misleadingly read as "cold" rather than "less severe heat".
severity_anchors <- c(
  "#FFFFB2",  # pale yellow   - mildest
  "#FECC5C",
  "#FD8D3C",
  "#F03B20",
  "#BD0026",  # dark red
  "#7A0177",
  "#49006A"   # deep purple   - most severe
)

# Interpolated in Lab space (perceptually smoother than default RGB) to any k,
# so the palette is not tied to exactly 7 clusters.
build_severity_pal <- function(k) {
  stats::setNames(
    grDevices::colorRampPalette(severity_anchors, space = "Lab")(k),
    as.character(seq_len(k))
  )
}

k_final <- dplyr::n_distinct(county_clustered$cluster)
severity_pal <- build_severity_pal(k_final)

cat("\nseverity palette (", k_final, "clusters):\n")
print(severity_pal)

##---- 3. static county map (tmap) #----

cluster_map <- tm_shape(counties_clustered_sf) +
  tm_polygons(
    col        = "cluster",
    title      = sprintf("Heatwave cluster\n(1 = mildest, %d = most severe)", k_final),
    palette    = severity_pal,
    border.col = NA,
    colorNA    = "grey90",
    textNA     = "No data"
  ) +
  tm_layout(
    main.title          = "Heatwave Clusters across CONUS Counties",
    main.title.position = "center",
    legend.outside      = TRUE,
    frame               = FALSE
  )

cluster_map

tmap_save(
  cluster_map,
  filename = "../Output/Figures/county_clusters_6.jpg",
  width    = 12,
  height   = 8,
  units    = "in",
  dpi      = 300
)

##---- 4. combined map + radar profiles (ggplot) #----
# Same cluster map rebuilt in ggplot (so it composes with patchwork), stacked
# over small-multiple radars of each cluster's mean profile on the portfolio
# variables chosen in the Cluster section.

cluster_mean_profile <- county_clustered %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    n_county = dplyr::n(),
    dplyr::across(dplyr::all_of(cluster_vars), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::arrange(cluster)

print(as.data.frame(cluster_mean_profile), row.names = FALSE)
write_csv(cluster_mean_profile,
          "processed/season cluster/cluster_mean_profile.csv")

# coord_polar() draws connecting lines as ARCS, which turns a radar into a
# crescent. Overriding is_linear makes the segments straight.
coord_radar <- function(theta = "x", start = 0, direction = 1) {
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") "y" else "x"
  ggproto("CoordRadar", CoordPolar,
          theta = theta, r = r, start = start,
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}

# Radial scale in z-units computed across COUNTIES, not across the cluster
# means, so axes with very different natural ranges (e.g. a ~1-2 day metric
# next to a ~10-30 day one) stay visually comparable.
cluster_z_plot <- county_clustered %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(cluster_vars), ~ as.numeric(scale(.x)))) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(cluster_vars), ~ mean(.x, na.rm = TRUE)),
                   .groups = "drop")

radar_dat <- cluster_z_plot %>%
  tidyr::pivot_longer(dplyr::all_of(cluster_vars), names_to = "dim", values_to = "z") %>%
  dplyr::left_join(cluster_mean_profile %>% dplyr::select(cluster, n_county), by = "cluster") %>%
  dplyr::mutate(
    dim       = factor(dim, levels = names(radar_vars), labels = radar_vars),
    facet_lab = factor(sprintf("Cluster %s\n(n = %d)", cluster, n_county),
                       levels = sprintf("Cluster %s\n(n = %d)",
                                        cluster_mean_profile$cluster,
                                        cluster_mean_profile$n_county))
  )

# Reference contour at the national mean. geom_hline would be drawn as a
# straight chord under the is_linear override, so trace it on the axes instead
# - that is the true zero contour in this coordinate system.
zero_ref <- radar_dat %>%
  dplyr::distinct(facet_lab, dim) %>%
  dplyr::mutate(z = 0)

p_cluster_map <- ggplot(counties_clustered_sf) +
  geom_sf(aes(fill = cluster), color = NA) +
  scale_fill_manual(values = severity_pal, na.value = "grey90",
                    name = "Heatwave cluster",
                    labels = function(x) ifelse(is.na(x), "No data", x),
                    drop = FALSE) +
  labs(title = "Heatwave clusters across CONUS counties") +
  # Type sized for a 6.5 in wide figure at 300 dpi. ggplot text sizes are
  # absolute points, so they do NOT rescale with the canvas.
  theme_void(base_size = 8) +
  theme(plot.background   = element_rect(fill = "white", color = NA),
        plot.title        = element_text(face = "bold", size = 10, hjust = 0.5),
        legend.position   = "right",
        legend.title      = element_text(size = 7.5, face = "bold"),
        legend.text       = element_text(size = 6.5),
        legend.key.height = unit(0.30, "cm"),
        legend.key.width  = unit(0.30, "cm"))

rad_lim <- c(min(radar_dat$z) - 0.45, max(radar_dat$z) + 0.25)

# Grid shape adapts to how many clusters were actually fit, rather than a
# hardcoded nrow = 2 that only looked right for exactly 7 panels.
facet_nrow <- 2

p_cluster_radar <- ggplot(radar_dat, aes(x = dim, y = z, group = cluster)) +
  # Stroke widths and point sizes are absolute (mm), scaled down alongside the
  # type for the 6.5 in canvas.
  geom_polygon(data = zero_ref, aes(x = dim, y = z, group = facet_lab),
               inherit.aes = FALSE, fill = NA, color = "grey55",
               linewidth = 0.25, linetype = "dashed") +
  geom_polygon(aes(fill = cluster, color = cluster), alpha = 0.4, linewidth = 0.45) +
  geom_point(aes(color = cluster), size = 0.9) +
  coord_radar() +
  facet_wrap(~ facet_lab, nrow = facet_nrow) +
  scale_fill_manual(values = severity_pal, guide = "none") +
  scale_color_manual(values = severity_pal, guide = "none") +
  scale_y_continuous(limits = rad_lim, breaks = c(-1, 0, 1)) +
  labs(
    title = "Cluster profiles across the portfolio variables",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 7) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    axis.text.x      = element_text(size = 5.5, face = "bold"),
    axis.text.y      = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.spacing    = unit(0.15, "cm"),
    strip.text       = element_text(face = "bold", size = 6.5, lineheight = 1.0),
    plot.title       = element_text(face = "bold", size = 9)
  )

cluster_map_profiles <- p_cluster_map / p_cluster_radar +
  plot_layout(heights = c(1.35, 1))

cluster_map_profiles

# 6.5 in = full text width of a US Letter page with 1 in margins.
ggsave(
  filename = "../Output/Figures/cluster_map_with_profiles_v3.png",
  plot     = cluster_map_profiles,
  width    = 6.5,
  height   = 8,
  units    = "in",
  dpi      = 600,
  bg       = "white"
)


#---- Descriptive analysis: 4 heatwave characteristics #----
# STANDALONE - reads heat_characteristics_v2.csv directly rather than reusing
# county_heat_characteristics/counties_conus/severity_pal etc. from the
# sections above, so this can be run on its own without re-running the whole
# pipeline (the 1.1 GB raw WBGT read, the per-county k-means season
# classification, the heatwave-spell detection). Only needs Setup (library
# loads + setwd()) to have already run, same as every other section's
# relative file paths in this script.
#   1. Four characteristic maps (choropleth, continuous scale, one per metric)
#   2. Pairwise correlation panel (Spearman - these metrics are right-skewed,
#      same reasoning as domain_exposure_heatmap.R section 7) - the point is
#      to show the four are not redundant, not just re-measuring one thing
#   3. Table of descriptive statistics
# The correlation panel and both tables are exported to one Word document
# (officer/flextable, same pattern already used in heat_vulnerability.R);
# maps are saved as PNG, matching this script's own map-export convention.

library(dplyr)
library(readr)
library(tidyr)
library(sf)
library(tigris)
library(ggplot2)
library(patchwork)
library(officer)
library(flextable)

options(tigris_use_cache = TRUE)

# Built via code points rather than typed literally, same reasoning as the
# rest of this project: Rscript in a shell with LANG unset parses the source
# as US-ASCII and would mangle a literal degree sign.
deg_sym <- intToUtf8(0x00B0)
rho_sym <- intToUtf8(0x03C1)

desc_vars <- c(
  "abs_hwdays_per_yr",
  "rel_intensity_tmax",
  "rel_duration_max",
  "rel_hw_span"
)
desc_labs <- c(
  abs_hwdays_per_yr  = paste0("Heatwave-day frequency (days/yr)"),
  rel_intensity_tmax = paste0("Max WBGT intensity (", deg_sym, "C)"),
  rel_duration_max   = paste0("Longest event duration (days)"),
  rel_hw_span        = paste0("Warm-season heatwave span (days)")
)

county_heat <- read_csv(
  "processed/season cluster/heat_characteristics_v2.csv",
  col_types = cols(StCoFIPS = col_character(), .default = col_double())
)
stopifnot(all(desc_vars %in% names(county_heat)))

cat("\ndescriptive section: county_heat characteristics loaded -",
    nrow(county_heat), "counties\n")

##---- 1. four characteristic maps #----
counties_desc_sf <- suppressMessages(
  counties(year = 2020, cb = TRUE, class = "sf")
) %>%
  mutate(GEOID = as.character(GEOID)) %>%
  filter(!STATEFP %in% c("02", "15", "60", "66", "69", "72", "78")) %>%
  st_transform(5070) %>%
  left_join(county_heat %>% dplyr::select(StCoFIPS, dplyr::all_of(desc_vars)),
           by = c("GEOID" = "StCoFIPS"))

# scale_fill_viridis_c() - a continuous, perceptually-uniform scale, kept
# deliberately distinct from severity_pal (categorical, cluster membership)
# used elsewhere in this script: these are raw continuous metrics, not
# cluster IDs, so they should not visually read as "which cluster".
desc_map <- function(var) {
  ggplot(counties_desc_sf) +
    geom_sf(aes(fill = .data[[var]]), color = NA) +
    scale_fill_viridis_c(option = "inferno", na.value = "grey85", name = NULL) +
    labs(title = desc_labs[var]) +
    theme_void(base_size = 9) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.title       = element_text(face = "bold", size = 9, hjust = 0.5),
      legend.key.width = unit(0.3, "cm"),
      legend.text      = element_text(size = 6)
    )
}

desc_maps <- lapply(desc_vars, desc_map)

p_desc_maps <- patchwork::wrap_plots(desc_maps, ncol = 2, nrow = 2) +
  patchwork::plot_annotation(
    title = "Heatwave characteristics across CONUS counties",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

p_desc_maps

ggsave(
  filename = "../Output/Figures/heat_characteristic_maps.png",
  plot     = p_desc_maps,
  width    = 9,
  height   = 7,
  dpi      = 600,
  bg       = "white"
)

##---- 2. pairwise correlation panel (Spearman) #----
# Just the correlation matrix, exported to Word (section 4 below) - showing
# whether the four characteristics are redundant with each other.
cor_mat <- cor(county_heat[desc_vars], method = "spearman", use = "complete.obs")

cat("\nSpearman correlation matrix (4 characteristics):\n")
print(round(cor_mat, 3))

##---- 3. table of descriptive statistics #----
desc_table <- county_heat %>%
  dplyr::select(dplyr::all_of(desc_vars)) %>%
  tidyr::pivot_longer(dplyr::everything(), names_to = "variable", values_to = "value") %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(
    N      = dplyr::n(),
    Mean   = mean(value),
    SD     = sd(value),
    Min    = min(value),
    Q1     = quantile(value, 0.25),
    Median = median(value),
    Q3     = quantile(value, 0.75),
    Max    = max(value),
    .groups = "drop"
  ) %>%
  dplyr::mutate(variable = desc_labs[variable]) %>%
  dplyr::arrange(match(variable, desc_labs[desc_vars])) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 2))) %>%
  dplyr::rename(Variable = variable)

cat("\ndescriptive statistics (4 characteristics):\n")
print(as.data.frame(desc_table), row.names = FALSE)

write_csv(desc_table, "processed/season cluster/heat_characteristic_descriptives.csv")
write_csv(as.data.frame(cor_mat) %>% tibble::rownames_to_column("variable"),
          "processed/season cluster/heat_characteristic_correlations.csv")

##---- 4. export tables + correlation panel to Word #----
# COLUMN NAMES (not values) containing the intToUtf8()-built degree sign make
# flextable() fail with "undefined columns selected" under this session's C
# locale - confirmed in isolation: the identical data.frame flextable()s fine
# when that character is only a CELL VALUE (as in desc_table's Variable
# column below), but errors as soon as it is part of a column NAME. ASCII-only
# column headers sidestep it; the row labels can keep the nicer degree sign
# since those are just values.
desc_labs_ascii <- gsub(deg_sym, "", desc_labs, fixed = TRUE)

cor_mat_tab <- as.data.frame(round(cor_mat, 3)) %>%
  tibble::rownames_to_column("Variable") %>%
  dplyr::mutate(Variable = unname(desc_labs[desc_vars]))
names(cor_mat_tab)[-1] <- unname(desc_labs_ascii[desc_vars])

doc <- read_docx()

doc <- body_add_par(doc, "Table 1. Descriptive statistics - 4 heatwave characteristics",
                    style = "heading 1")
doc <- body_add_flextable(doc, autofit(flextable(desc_table)))

doc <- body_add_par(doc, "Table 2. Spearman correlation matrix", style = "heading 1")
doc <- body_add_flextable(doc, autofit(flextable(cor_mat_tab)))

docx_path <- "../Output/Figures/heat_characteristics_descriptive_report.docx"
print(doc, target = docx_path)

cat("\nexported:\n")
cat("  ../Output/Figures/heat_characteristic_maps.png\n")
cat("  processed/season cluster/heat_characteristic_descriptives.csv\n")
cat("  processed/season cluster/heat_characteristic_correlations.csv\n")
cat("  ", docx_path, "\n", sep = "")

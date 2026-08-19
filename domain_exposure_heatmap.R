#---- Script Metadata #----
# Title: Domain x Exposure heatmap
# Author: Hang Li
# Purpose: Replace the stacked-bar PCA profile figure with a domain x exposure
#          matrix. The stack sums the 7 components, which cancels their opposing
#          relationships with heat exposure (cumulative index eta2 = 2.7% vs
#          13.5% for greenness and 6.6% for socioeconomic alone). This figure
#          shows each domain separately so the sign reversals stay visible.
#---------------------------------------------------------------------------

#---- 1. Setup #----
rm(list = ls())

pkgs <- c("dplyr", "readr", "tidyr", "ggplot2", "ragg")
invisible(lapply(pkgs, library, character.only = TRUE))

# The default macOS png device drops UTF-8 glyphs (eta-squared, degree sign).
# ragg::agg_png renders them correctly.
png_dev <- ragg::agg_png

# Built via code points rather than typed literally: Rscript launched from a
# shell with LANG unset parses the source as US-ASCII and would mangle literal
# UTF-8 bytes. This keeps the source pure ASCII and is locale-independent.
eta2_sym  <- intToUtf8(c(0x03B7, 0x00B2))   # eta-squared
deg_sym   <- intToUtf8(0x00B0)              # degree sign
arrow_sym <- intToUtf8(0x2192)              # rightwards arrow
rho_sym   <- intToUtf8(0x03C1)              # greek rho
dash_sym  <- intToUtf8(0x2014)              # em dash

data_dir <- "/Users/hangli/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Aune, Kyle T's files - 260204 - US Heat Clusters/Data"
fig_dir  <- "/Users/hangli/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Aune, Kyle T's files - 260204 - US Heat Clusters/Output/Figures"

season_dir <- file.path(data_dir, "processed", "season cluster")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


#---- 2. Data #----
# Per-tract component scores from heat_vulnerability.R section 9. That solution
# is MSA- and communality-screened (16 of 19 variables, 5 components, 78.3% of
# variance) and supersedes the older tract_pca.csv, which came from an
# unscreened 19-variable / 7-component PCA whose 7th component was a single
# variable and whose file carried unscored rows.
pca <- read_csv(
  file.path(season_dir, "pca_tract_scores.csv"),
  col_types = cols(GEOID = col_character(), .default = col_double())
) %>%
  mutate(StCoFIPS = substr(GEOID, 1, 5))

county <- read_csv(
  file.path(season_dir, "county_clustered.csv"),
  col_types = cols(StCoFIPS = col_character(),
                   cluster  = col_integer(),
                   .default = col_double())
)

# county_clustered.csv now carries `cluster` already relabeled 1 (mildest) ...
# 7 (most severe) - that relabeling happens in heatwave_characteristics.R's
# "#---- Cluster #----" section (step 5b), not here, so no crosswalk join is
# needed any more: the file used to export the RAW k-means id and rely on a
# separate cluster_severity_crosswalk.csv to relabel it, which could (and did)
# go stale relative to county_clustered.csv after a re-run. `cluster_raw`
# survives in the file for traceability but is not used below.
#
# heat_characteristics_v2.csv was reworked into 8 metrics x {rel, abs} =
# 16 columns (frequency/intensity/duration/timing under two heatwave
# definitions); the old 4-column set (frequency_pct, mean_Tmax_intensity,
# mean_Tmin_intensity, longest_duration) no longer exists.
#
# Continuous per-county hazard composite (z-scored across counties), used for
# sections 9-11's tercile map / LISA where hazard has to vary continuously,
# not just take 7 discrete cluster values. Built from the SAME four variables
# as heatwave_characteristics.R's `cluster_vars` (the portfolio that drove the
# k-means fit and the mild -> severe ranking), so this script's hazard and
# that script's cluster order describe the same thing. If cluster_vars changes
# there, update hazard_vars here to match - the check right below will catch a
# mismatch (loudly) rather than silently drifting apart again.
hazard_vars <- c("abs_hwdays_per_yr", "rel_intensity_tmax",
                 "rel_duration_max", "rel_hw_span")
stopifnot(all(hazard_vars %in% names(county)))

county <- county %>%
  mutate(
    hazard = as.numeric(scale(abs_hwdays_per_yr)) +
             as.numeric(scale(rel_intensity_tmax)) +
             as.numeric(scale(rel_duration_max)) +
             as.numeric(scale(rel_hw_span))
  )

hazard_check <- county %>%
  group_by(cluster) %>%
  summarise(mean_hazard = mean(hazard), .groups = "drop") %>%
  arrange(cluster)
if (!all(diff(hazard_check$mean_hazard) > 0)) {
  warning("mean hazard is not monotonically increasing across cluster 1..",
          max(county$cluster), " - hazard_vars here may no longer match ",
          "cluster_vars in heatwave_characteristics.R")
}

d <- pca %>%
  inner_join(county, by = "StCoFIPS")

message("tracts matched: ", nrow(d),
        " across ", n_distinct(d$StCoFIPS), " counties")

# Component order comes from the exported variance table rather than being
# hardcoded: psych names components by extraction order but reports them by
# variance, so the columns are RC1, RC2, RC3, RC5, RC4 - not sequential.
pca_var <- read_csv(file.path(season_dir, "pca_variance_explained.csv"),
                    col_types = cols(component = col_character(),
                                     .default = col_double()))
pcs <- pca_var$component
stopifnot(all(pcs %in% names(pca)))

pca_load <- read_csv(file.path(season_dir, "pca_loadings_contributions.csv"),
                     col_types = cols(Variable = col_character(),
                                      .default = col_double()))

L <- as.matrix(pca_load[, paste0(pcs, "_loading")])
rownames(L) <- pca_load$Variable
colnames(L) <- pcs

# Readable labels, with an automatic fallback to the two leading variables so a
# changed variable set degrades to accurate-but-ugly rather than to wrong.
#
# Re-derived from the current pca_loadings_contributions.csv (17 variables:
# the screened land-cover/socioeconomic set plus the new pct_disability,
# pct_nonwhite, pct_living_alone from heat_vulnerability.R). Component
# meanings shifted from the previous run - dropping albedo removed the
# separate "greenness vs albedo" axis, so tree_canopy/ndvi now anchor RC1
# together with mid-density development instead; the extra development/
# population-density variables split into a second, separate "dense urban
# core" component (RC5) that did not exist as its own axis before:
#   RC2 (20.2%) income, poverty, energy/electricity burden, disability
#   RC1 (19.6%) low tree canopy/NDVI + medium-density development
#   RC5 (17.4%) high-density development + population density
#   RC3 ( 9.6%) AC_central + AC_none, essentially nothing else
#   RC4 ( 9.3%) age 65+, disability, living alone (vs. LOWER pct_nonwhite -
#               bipolar, see the check below; label describes the majority
#               of the salient loadings, not all of them)
pc_labels_manual <- c(
  RC1 = "Low greenness",
  RC2 = "Socioeconomic disadvantage",
  RC3 = "Low AC access",
  RC4 = "Demographic susceptibility",
  RC5 = "Built density"
)
pc_labels <- vapply(pcs, function(p) {
  if (p %in% names(pc_labels_manual)) return(unname(pc_labels_manual[p]))
  paste(rownames(L)[order(abs(L[, p]), decreasing = TRUE)[1:2]], collapse = " + ")
}, character(1))
names(pc_labels) <- pcs


#---- 3. Order exposure clusters along the hazard gradient #----
# `cluster` already IS the mild -> severe order (1..7) coming out of
# heatwave_characteristics.R, so this just arranges by it and builds display
# labels - no re-ranking happens here. Column stats shown are the same two
# variables (frequency, max intensity) that headline the radar chart there,
# for a consistent read across both scripts' figures.
cluster_info <- d %>%
  group_by(cluster) %>%
  summarise(
    hazard    = mean(hazard),
    hwdays    = mean(abs_hwdays_per_yr),
    tmax      = mean(rel_intensity_tmax),
    dur_max   = mean(rel_duration_max),
    hw_span   = mean(rel_hw_span),
    n_county  = n_distinct(StCoFIPS),
    n_tract   = n(),
    .groups   = "drop"
  ) %>%
  arrange(cluster) %>%
  mutate(
    rank  = cluster,
    label = sprintf("Cluster %d\n%.1f d/yr, %.1f%sC\n(%d co.)",
                    cluster, hwdays, tmax, deg_sym, n_county)
  )

print(as.data.frame(cluster_info), row.names = FALSE)


#---- 4. Effect size per domain #----
# eta-squared = share of tract-level PC variance explained by the exposure
# typology. With n ~ 84k every p-value is < 1e-16, so report effect size instead.
eta_sq <- function(x, g) {
  ok <- !is.na(x); x <- x[ok]; g <- g[ok]
  gm <- mean(x)
  ss_between <- sum(tapply(x, g, function(v) length(v) * (mean(v) - gm)^2))
  ss_between / sum((x - gm)^2)
}

effect <- tibble(
  domain = pcs,
  eta2   = sapply(pcs, function(p) eta_sq(d[[p]], d$cluster))
) %>%
  arrange(desc(eta2))


#---- 5. Build the matrix #----
# PC scores from psych::principal are already unit-variance, so these cell
# values are directly interpretable as SD units relative to the national mean.
mat <- d %>%
  dplyr::select(cluster, all_of(pcs)) %>%
  pivot_longer(all_of(pcs), names_to = "domain", values_to = "score") %>%
  group_by(cluster, domain) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    se         = sd(score, na.rm = TRUE) / sqrt(sum(!is.na(score))),
    .groups    = "drop"
  ) %>%
  left_join(effect, by = "domain") %>%
  left_join(cluster_info %>% dplyr::select(cluster, rank, label), by = "cluster") %>%
  mutate(
    # rows ordered by effect size, strongest at top
    domain_lab  = factor(
      sprintf("%s\n(%s = %.1f%%)", pc_labels[domain], eta2_sym, 100 * eta2),
      levels = rev(sprintf("%s\n(%s = %.1f%%)",
                           pc_labels[effect$domain], eta2_sym, 100 * effect$eta2))
    ),
    cluster_lab = factor(label, levels = cluster_info$label)
  )

write_csv(
  mat %>% dplyr::select(cluster, rank, domain, mean_score, se, eta2),
  file.path(season_dir, "domain_exposure_matrix.csv")
)


#---- 6. Plot #----
lim <- max(abs(mat$mean_score))

brks <- c(-0.45, -0.30, -0.15, 0, 0.15, 0.30, 0.45)   # 7 breaks -> 8 classes

p_heat <- ggplot(mat, aes(x = cluster_lab, y = domain_lab, fill = mean_score)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(abs(round(mean_score, 2)) < 0.005,
                               "0.00", sprintf("%+.2f", mean_score)),
                color = abs(mean_score) >= 0.45),
            size = 3.1, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")) +
  scale_fill_fermenter(
    palette   = "RdBu",
    direction = -1,
    breaks    = brks,
    limits    = c(-0.6, 0.8),
    name      = "Mean PCA score\n(SD units)",
    guide     = guide_coloursteps(
      barheight     = unit(55, "mm"),
      barwidth      = unit(4.5, "mm"),
      show.limits   = TRUE,
      even.steps    = TRUE,
      title.position = "top"
    )
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x        = paste0("Heatwave exposure cluster"),
    y        = NULL,
    title    = "Heat vulnerability domains across heatwave exposure clusters"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid       = element_blank(),
    axis.text.x.top  = element_text(size = 8.5, lineheight = 1.05),
    axis.text.y      = element_text(size = 9, lineheight = 1.05, hjust = 1),
    axis.title.x     = element_text(margin = margin(t = 10), face = "bold"),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 8.7,
                                    margin = margin(b = 10)),
    plot.caption     = element_text(size = 7.6, hjust = 0,
                                    margin = margin(t = 10)),
    legend.position  = "right",
    legend.key.height = unit(1.1, "cm")
  )

if (interactive()) print(p_heat)

ggsave(
  filename = file.path(fig_dir, "domain_exposure_heatmap_6.png"),
  plot     = p_heat,
  width    = 10,
  height   = 5,
  dpi      = 300,
  device   = png_dev
)


#---- 7. Companion: continuous exposure vs domain (no binning) #----
# Same story without collapsing counties into 7 bins. Spearman is used because
# the exposure metrics are right-skewed. Uses the same four hazard_vars as the
# cluster composite above (not the full 16-metric set), so this panel is the
# continuous-variable version of exactly what produced the clusters.
expo_vars <- c(abs_hwdays_per_yr  = "Heatwave-day frequency",
               rel_intensity_tmax = "Max WBGT intensity",
               rel_duration_max   = "Longest event duration",
               rel_hw_span        = "Warm-season heatwave span")

corr_tab <- expand_grid(domain = pcs, expo = names(expo_vars)) %>%
  rowwise() %>%
  mutate(rho = cor(d[[domain]], d[[expo]],
                   method = "spearman", use = "complete.obs")) %>%
  ungroup() %>%
  mutate(
    domain_lab = factor(pc_labels[domain], levels = rev(pc_labels[effect$domain])),
    expo_lab   = factor(expo_vars[expo],   levels = expo_vars)
  )

p_corr <- ggplot(corr_tab, aes(x = rho, y = domain_lab, fill = rho)) +
  geom_vline(xintercept = 0, color = "grey55", linewidth = 0.4) +
  geom_col(width = 0.65) +
  facet_wrap(~ expo_lab, nrow = 1) +
  scale_fill_gradient2(low = "#2166AC", mid = "grey90", high = "#B2182B",
                       midpoint = 0, guide = "none") +
  labs(
    x        = paste0("Spearman ", rho_sym, "  (tract level)"),
    y        = NULL,
    title    = "Vulnerability domains vs. continuous heatwave exposure",
    subtitle = paste0("Sign reversals across domains ", dash_sym,
                      " the effect the cumulative index averages away")
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold")
  )

if (interactive()) print(p_corr)

ggsave(
  filename = file.path(fig_dir, "domain_exposure_correlations.png"),
  plot     = p_corr,
  width    = 11,
  height   = 4.2,
  dpi      = 600,
  device   = png_dev
)

write_csv(corr_tab %>% dplyr::select(domain, expo, rho),
          file.path(season_dir, "domain_exposure_correlations.csv"))

# --- regime-stratified correlations: All vs each cluster ---------------------
# Same Spearman correlations as above, computed separately WITHIN each
# heatwave exposure cluster ("regime"), not just pooled across all of them.
# A domain x exposure relationship that looks weak (or a particular sign) in
# the all-sample result could still be strong, weak, or reversed inside a
# specific regime - averaging across clusters can mask that the same
# relationship. cor.test() adds a p-value alongside rho, since regime sample
# sizes are much smaller than the full tract set and worth knowing whether a
# given regime's correlation is actually distinguishable from zero.
regime_corr <- function(dat, regime_lab) {
  expand_grid(domain = pcs, expo = names(expo_vars)) %>%
    rowwise() %>%
    mutate(
      n     = sum(complete.cases(dat[[domain]], dat[[expo]])),
      rho   = cor(dat[[domain]], dat[[expo]], method = "spearman", use = "complete.obs"),
      p_val = suppressWarnings(
        cor.test(dat[[domain]], dat[[expo]], method = "spearman")$p.value
      )
    ) %>%
    ungroup() %>%
    mutate(regime = regime_lab)
}

corr_all_tab <- regime_corr(d, "All")

corr_regime_tab <- lapply(sort(unique(d$cluster)), function(cl) {
  regime_corr(d %>% dplyr::filter(cluster == cl), as.character(cl))
}) %>%
  dplyr::bind_rows()

corr_compare <- dplyr::bind_rows(corr_all_tab, corr_regime_tab) %>%
  dplyr::mutate(
    domain_lab = factor(pc_labels[domain], levels = pc_labels[effect$domain]),
    expo_lab   = factor(expo_vars[expo], levels = expo_vars),
    regime     = factor(regime, levels = c("All", as.character(sort(unique(d$cluster)))))
  ) %>%
  dplyr::arrange(domain_lab, expo_lab, regime)

cat("\nregime-stratified Spearman correlations (domain x exposure), All + 7 clusters:\n")
print(as.data.frame(corr_compare %>%
                      dplyr::select(domain_lab, expo_lab, regime, n, rho, p_val) %>%
                      dplyr::mutate(rho = round(rho, 3), p_val = signif(p_val, 3))),
      row.names = FALSE)

write_csv(corr_compare %>% dplyr::select(domain, expo, regime, n, rho, p_val),
          file.path(season_dir, "domain_exposure_correlations_by_regime.csv"))

# Wide comparison table - one row per domain x exposure pair, one column per
# regime (All, 1..7) - is the actual "compare all-sample vs regime-specific"
# view, read left-to-right across a row.
corr_wide_rho <- corr_compare %>%
  dplyr::select(domain_lab, expo_lab, regime, rho) %>%
  tidyr::pivot_wider(names_from = regime, values_from = rho) %>%
  dplyr::arrange(domain_lab, expo_lab) %>%
  dplyr::rename(Domain = domain_lab, Exposure = expo_lab) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x, 3)))

cat("\nwide comparison table (Spearman rho), All vs regime 1-7:\n")
print(as.data.frame(corr_wide_rho), row.names = FALSE)

# --- export to Word (same officer/flextable pattern already used in --------
# heat_vulnerability.R's PCA export, for consistency across this project's
# scripts).
suppressMessages(library(officer))
suppressMessages(library(flextable))

corr_long_tab <- corr_compare %>%
  dplyr::transmute(Domain = domain_lab, Exposure = expo_lab, Regime = regime,
                   N = n, Rho = round(rho, 3), P = signif(p_val, 3))

doc <- read_docx()
doc <- body_add_par(
  doc, "Table 1. Domain x exposure Spearman correlations - all samples vs regime-stratified (wide)",
  style = "heading 1"
)
doc <- body_add_flextable(doc, autofit(flextable(corr_wide_rho)))
doc <- body_add_par(
  doc, "Table 2. Domain x exposure Spearman correlations - all samples vs regime-stratified (long, with n and p)",
  style = "heading 1"
)
doc <- body_add_flextable(doc, autofit(flextable(corr_long_tab)))

docx_path <- file.path(fig_dir, "domain_exposure_correlations_by_regime.docx")
print(doc, target = docx_path)
cat("\nexported:", docx_path, "\n")


#---- 7b. Domain score distributions across clusters (violin) #----
# Section 6's heatmap collapses each cluster x domain cell to one mean score,
# which hides shape: a tight, unimodal distribution and a wide, bimodal one
# can produce the same mean. This shows the full tract-level distribution
# instead - one panel PER DOMAIN (5 panels, 3+2 over 2 rows), 7 violins per
# panel (the mild -> severe exposure clusters), so within a domain you can
# compare the clusters directly. Each panel has its own y-axis: the domains
# differ enough in spread that a shared axis would flatten the narrower ones.

violin_dat <- d %>%
  dplyr::select(cluster, dplyr::all_of(pcs)) %>%
  pivot_longer(dplyr::all_of(pcs), names_to = "domain", values_to = "score") %>%
  dplyr::filter(!is.na(score)) %>%
  left_join(cluster_info %>% dplyr::select(cluster, rank), by = "cluster") %>%
  mutate(
    # same effect-size order as the heatmap rows, so the two figures agree
    domain_lab  = factor(pc_labels[domain], levels = pc_labels[effect$domain]),
    cluster_ord = factor(rank, levels = seq_len(nrow(cluster_info)))
  )

# Same mild (pale yellow) -> severe (deep purple) ramp used for the cluster
# maps in heatwave_characteristics.R and heatwave_cluster_methods.R, so a
# cluster's color means the same thing in every figure in this project. Not
# defined elsewhere in this script, so it is declared here rather than
# assumed to already exist.
severity_anchors <- c(
  "#FFFFB2",  # pale yellow   - mildest
  "#FECC5C",
  "#FD8D3C",
  "#F03B20",
  "#BD0026",  # dark red
  "#7A0177",
  "#49006A"   # deep purple   - most severe
)
build_severity_pal <- function(k) {
  stats::setNames(grDevices::colorRampPalette(severity_anchors, space = "Lab")(k),
                  as.character(seq_len(k)))
}
severity_pal <- build_severity_pal(nrow(cluster_info))

p_violin <- ggplot(violin_dat, aes(x = cluster_ord, y = score, fill = cluster_ord)) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.2, alpha = 0.9) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white",
              alpha = 0.7, linewidth = 0.3) +
  facet_wrap(~ domain_lab, nrow = 2, scales = "free_y") +
  scale_fill_manual(values = severity_pal, guide = "none") +
  labs(
    x        = paste0("Heatwave exposure cluster (1 = mildest ", arrow_sym, " 7 = most severe)"),
    y        = "Tract-level PCA score (SD units)",
    title    = "Vulnerability domain distributions across heatwave exposure clusters",
    subtitle = paste0(
      "One panel per domain (ordered by ", eta2_sym, "), 7 clusters per panel, ",
      "colored mild ", arrow_sym, " severe as in the cluster maps.\n",
      "Y-axes are independent per panel - only compare shape/spread WITHIN a ",
      "panel, not the absolute scale ACROSS panels."
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor    = element_blank(),
    strip.text          = element_text(face = "bold", size = 9),
    plot.title           = element_text(face = "bold", size = 12),
    plot.subtitle        = element_text(size = 8.3, color = "grey30"),
    axis.text.x          = element_text(size = 8)
  )

if (interactive()) print(p_violin)

ggsave(
  filename = file.path(fig_dir, "domain_score_violins_by_domain.png"),
  plot     = p_violin,
  width    = 11,
  height   = 8.5,
  dpi      = 600,
  device   = png_dev
)

violin_summary <- violin_dat %>%
  dplyr::group_by(domain, cluster_ord) %>%
  dplyr::summarise(
    n      = dplyr::n(),
    mean   = mean(score),
    median = median(score),
    sd     = sd(score),
    iqr    = IQR(score),
    .groups = "drop"
  ) %>%
  dplyr::arrange(domain, cluster_ord)

write_csv(violin_summary, file.path(season_dir, "domain_score_violin_summary.csv"))


#---- 7c. Composite HVI vs domain distributions across clusters (boxplot) #----
# Composite heat vulnerability index (HVI: equal-weight sum of the 5 oriented
# domain z-scores, then re-standardized) compared against its own 5 component
# domains, both across the 7 heatwave exposure clusters. Same construction as
# section 8's vuln_tract$hvi, computed here directly from d's own pcs columns
# so this section stays self-contained and does not need section 8 to have
# run first (matches 7b/7c's established independence from later sections).
#
# Boxplot only (no violin), same severity_pal color ramp as 7b/7d/7e/7f, no
# significance testing. Two stacked panels in one figure:
#   top:    composite HVI, 1 panel, 7 clusters
#   bottom: the 5 domains that make it up, one panel each (same order as 7b),
#           7 clusters per panel - a boxplot-only version of 7b

# patchwork's `/` stacking operator only works once the package is attached
# (namespacing individual calls with patchwork:: is not enough) - loaded here
# since this is the first section that needs it in file execution order.
suppressMessages(library(patchwork))

hvi_vec <- as.numeric(scale(rowSums(as.matrix(d[pcs]))))

hvi_dat <- d %>%
  dplyr::select(cluster) %>%
  dplyr::mutate(hvi = hvi_vec) %>%
  dplyr::filter(!is.na(hvi)) %>%
  left_join(cluster_info %>% dplyr::select(cluster, rank), by = "cluster") %>%
  mutate(cluster_ord = factor(rank, levels = seq_len(nrow(cluster_info))))

domain_box_dat <- d %>%
  dplyr::select(cluster, dplyr::all_of(pcs)) %>%
  pivot_longer(dplyr::all_of(pcs), names_to = "domain", values_to = "score") %>%
  dplyr::filter(!is.na(score)) %>%
  left_join(cluster_info %>% dplyr::select(cluster, rank), by = "cluster") %>%
  mutate(
    # same effect-size order as 7b/section 6's heatmap, so all three figures agree
    domain_lab  = factor(pc_labels[domain], levels = pc_labels[effect$domain]),
    cluster_ord = factor(rank, levels = seq_len(nrow(cluster_info)))
  )

cat("\ncomposite HVI, tracts by cluster:\n")
print(table(hvi_dat$cluster_ord))

# y-axis zoomed to the Tukey whisker range (Q1-1.5*IQR to Q3+1.5*IQR) plus 8%
# padding, not the full data range. outlier.shape = NA hides outlier POINTS
# but does NOT shrink ggplot2's computed axis range on its own (confirmed by
# testing directly with ggplot_build() in 7f) - a handful of extreme tracts
# were stretching every panel's axis far past where its box+whiskers actually
# sit, flattening the boxes down to slivers. Same fix as 7f's vuln_ylim(),
# duplicated here rather than shared, since 7f runs later in this file and
# this section is meant to stay self-contained.
box_ylim <- function(x) {
  s   <- boxplot.stats(x)$stats   # Tukey 5-number: lo whisker..hi whisker
  pad <- 0.08 * diff(range(s))
  if (pad == 0) pad <- 1e-6 * max(abs(s), 1)  # guard against a degenerate/constant panel
  c(s[1] - pad, s[5] + pad)
}

p_hvi_box <- ggplot(hvi_dat, aes(x = cluster_ord, y = hvi, fill = cluster_ord)) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3) +
  coord_cartesian(ylim = box_ylim(hvi_dat$hvi)) +
  scale_fill_manual(values = severity_pal, guide = "none") +
  labs(
    x     = NULL,
    y     = "Composite HVI (SD units)",
    title = "Composite HVI"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title         = element_text(face = "bold", size = 11),
    axis.text.x        = element_text(size = 8)
  )

# facet_wrap(scales = "free_y") cannot give each panel its OWN
# coord_cartesian(ylim = ...) - that argument applies to the whole plot, not
# per facet - so each domain is its own small ggplot (whisker-zoomed to that
# domain's own range) instead, combined with patchwork in the same 3+2 grid
# 7b uses for its 5 domain panels.
domain_box_plot <- function(lab) {
  dat_d <- domain_box_dat %>% dplyr::filter(domain_lab == lab)
  ggplot(dat_d, aes(x = cluster_ord, y = score, fill = cluster_ord)) +
    geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.3) +
    coord_cartesian(ylim = box_ylim(dat_d$score)) +
    scale_fill_manual(values = severity_pal, guide = "none") +
    labs(x = NULL, y = NULL, title = lab) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title         = element_text(face = "bold", size = 9),
      axis.text.x        = element_text(size = 7.5)
    )
}

domain_box_plots <- lapply(levels(domain_box_dat$domain_lab), domain_box_plot)

p_domain_box <- patchwork::wrap_plots(domain_box_plots, ncol = 3, nrow = 2) +
  patchwork::plot_annotation(
    title = "5 vulnerability domains",
    theme = theme(plot.title = element_text(face = "bold", size = 11))
  )

p_hvi_vs_domains <- (p_hvi_box / p_domain_box) +
  patchwork::plot_layout(heights = c(1, 2)) +
  patchwork::plot_annotation(
    title    = "Composite HVI vs. its 5 component domains across heatwave exposure clusters",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 8.3, color = "grey30")
    )
  )

if (interactive()) print(p_hvi_vs_domains)

ggsave(
  filename = file.path(fig_dir, "hvi_vs_domains_boxplot_by_cluster.png"),
  plot     = p_hvi_vs_domains,
  width    = 11,
  height   = 10,
  dpi      = 600,
  device   = png_dev
)

# One-way ANOVA per variable (composite HVI + its 5 component domains),
# regime membership (heatwave exposure cluster) as the only predictor.
# Eta-squared = SS(cluster_ord) / SS(total), read directly off the ANOVA
# table - the proportion of that variable's total tract-level variance
# explained by which cluster a tract falls in. With n in the tens of
# thousands, p-values are essentially always < 2.2e-16 regardless of how
# small the effect is (same point already made for the domain-level version
# of this in section 4's `effect` table) - eta-squared is the number worth
# reading, not significance.
anova_dat <- dplyr::bind_rows(
  hvi_dat        %>% dplyr::transmute(variable = "Composite HVI", value = hvi, cluster_ord),
  domain_box_dat %>% dplyr::transmute(variable = as.character(domain_lab), value = score, cluster_ord)
)

regime_anova <- anova_dat %>%
  dplyr::group_by(variable) %>%
  dplyr::group_modify(~{
    fit       <- aov(value ~ cluster_ord, data = .x)
    tab       <- summary(fit)[[1]]
    ss_regime <- tab["cluster_ord", "Sum Sq"]
    ss_resid  <- tab["Residuals",   "Sum Sq"]
    tibble::tibble(
      df_regime = tab["cluster_ord", "Df"],
      df_resid  = tab["Residuals",   "Df"],
      f_value   = tab["cluster_ord", "F value"],
      p_value   = tab["cluster_ord", "Pr(>F)"],
      eta_sq    = ss_regime / (ss_regime + ss_resid)
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(dplyr::desc(eta_sq))

cat("\n=== one-way ANOVA: value ~ cluster_ord, composite HVI + 5 domains ===\n")
print(as.data.frame(regime_anova %>%
                      dplyr::mutate(f_value = round(f_value, 1),
                                   eta_sq  = round(eta_sq, 4),
                                   p_value = signif(p_value, 3))),
      row.names = FALSE)

write_csv(regime_anova, file.path(season_dir, "hvi_domain_regime_anova.csv"))

hvi_summary <- hvi_dat %>%
  dplyr::group_by(cluster_ord) %>%
  dplyr::summarise(n = dplyr::n(), mean = mean(hvi), median = median(hvi),
                   sd = sd(hvi), iqr = IQR(hvi), .groups = "drop")

domain_box_summary <- domain_box_dat %>%
  dplyr::group_by(domain, cluster_ord) %>%
  dplyr::summarise(n = dplyr::n(), mean = mean(score), median = median(score),
                   sd = sd(score), iqr = IQR(score), .groups = "drop") %>%
  dplyr::arrange(domain, cluster_ord)

write_csv(hvi_summary, file.path(season_dir, "hvi_boxplot_summary.csv"))
write_csv(domain_box_summary, file.path(season_dir, "domain_boxplot_summary.csv"))


#---- 7d. Domain radar profiles by cluster #----
# Radar of each cluster's mean profile across the 5 vulnerability domains,
# same technique as heatwave_characteristics.R's "##---- 4. combined map +
# radar profiles (ggplot) #----" section: a custom coord_radar() so polygon
# edges are straight (coord_polar() alone draws them as arcs), z-scored axes
# so domains with different natural spread stay visually comparable, one
# small-multiple panel per cluster with a dashed zero-reference contour.
#
# TRACT level (uses `d`, one row per tract; pcs = the 5 oriented PCA domain
# scores already built in section 2) - this section now profiles
# VULNERABILITY shape by cluster, not the heat-exposure metrics that used to
# live here as a violin; those are 7f's territory now.

# coord_polar() draws connecting lines as ARCS, which turns a radar into a
# crescent. Overriding is_linear makes the segments straight. Duplicated from
# heatwave_characteristics.R (a standalone ggproto helper, not from a
# package) rather than sourcing that whole script for one function.
coord_radar <- function(theta = "x", start = 0, direction = 1) {
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") "y" else "x"
  ggproto("CoordRadar", CoordPolar,
          theta = theta, r = r, start = start,
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}

cluster_domain_profile <- d %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    n_tract = dplyr::n(),
    dplyr::across(dplyr::all_of(pcs), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::arrange(cluster)

print(as.data.frame(cluster_domain_profile), row.names = FALSE)
# write_csv(cluster_domain_profile, file.path(season_dir, "cluster_domain_profile.csv"))

# Radial scale in z-units computed across TRACTS, not across the cluster
# means, so all 5 domains stay visually comparable regardless of each one's
# natural spread (pcs are already close to unit variance out of the PCA, but
# not necessarily exactly - re-scaling here rather than assuming it matches
# the same explicit-over-assumed convention used for the 7e score bins).
cluster_z_plot <- d %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(pcs), ~ as.numeric(scale(.x)))) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(pcs), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

radar_dat <- cluster_z_plot %>%
  tidyr::pivot_longer(dplyr::all_of(pcs), names_to = "dim", values_to = "z") %>%
  dplyr::left_join(cluster_domain_profile %>% dplyr::select(cluster, n_tract), by = "cluster") %>%
  dplyr::mutate(
    # Wrapped onto 2 lines - "Socioeconomic disadvantage" and "Demographic
    # susceptibility" run wide enough as one line that adjacent axis labels
    # collided/clipped in the unwrapped version (4 panels per row in an
    # 11 in figure), especially the "Low greenness" label sitting at the
    # right edge of each panel.
    dim       = factor(dim, levels = pcs,
                       labels = vapply(pc_labels[pcs], function(x)
                         paste(strwrap(x, width = 12), collapse = "\n"), character(1))),
    facet_lab = factor(sprintf("Cluster %s\n(n = %s tracts)", cluster,
                               format(n_tract, big.mark = ",")),
                       levels = sprintf("Cluster %s\n(n = %s tracts)",
                                        cluster_domain_profile$cluster,
                                        format(cluster_domain_profile$n_tract, big.mark = ",")))
  )

# Reference contour at the tract-level mean (z = 0). geom_hline would be drawn
# as a straight chord under the is_linear override, so trace it on the axes
# instead - that is the true zero contour in this coordinate system.
zero_ref <- radar_dat %>%
  dplyr::distinct(facet_lab, dim) %>%
  dplyr::mutate(z = 0)

rad_lim <- c(min(radar_dat$z) - 0.45, max(radar_dat$z) + 0.25)

p_domain_radar <- ggplot(radar_dat, aes(x = dim, y = z, group = cluster)) +
  geom_polygon(data = zero_ref, aes(x = dim, y = z, group = facet_lab),
              inherit.aes = FALSE, fill = NA, color = "grey55",
              linewidth = 0.25, linetype = "dashed") +
  geom_polygon(aes(fill = factor(cluster), color = factor(cluster)),
              alpha = 0.4, linewidth = 0.45) +
  geom_point(aes(color = factor(cluster)), size = 0.9) +
  coord_radar() +
  facet_wrap(~ facet_lab, nrow = 2) +
  scale_fill_manual(values = severity_pal, guide = "none") +
  scale_color_manual(values = severity_pal, guide = "none") +
  scale_y_continuous(limits = rad_lim, breaks = c(-1, 0, 1)) +
  labs(
    title    = "Vulnerability domain profiles by heatwave exposure cluster",
    subtitle = paste0(
      "Mean z-scored domain score per cluster, 5 domains per radar. ",
      "Dashed contour = tract-level mean (z = 0)."
    ),
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 9) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    axis.text.x      = element_text(size = 6.2, face = "bold", lineheight = 0.85),
    axis.text.y      = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.spacing    = unit(0.4, "cm"),
    strip.text       = element_text(face = "bold", size = 7.5, lineheight = 1.0),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 8.3, color = "grey30")
  )

if (interactive()) print(p_domain_radar)

ggsave(
  filename = file.path(fig_dir, "domain_radar_by_cluster.png"),
  plot     = p_domain_radar,
  width    = 13,
  height   = 8,
  dpi      = 300,
  device   = png_dev
)


#---- 7e. Within-cluster inequality (Gini coefficient) by domain, radar #----
# For each of the 5 vulnerability domains, how UNEQUALLY is that domain's
# score distributed among the tracts inside each cluster? Same radar
# technique as 7d (custom coord_radar(), small-multiple per cluster), but the
# axis value is the Gini coefficient of the domain score computed WITHIN each
# cluster's own tracts, not a mean.
#
# The classic Gini coefficient (area between the Lorenz curve and the line of
# equality) is only defined for non-negative values - pcs are PCA component
# scores (z-scored, signed), so DescTools::Gini() returns NA on them directly
# (confirmed by testing). Each domain is therefore min-max rescaled to [0, 1]
# using that domain's GLOBAL min/max across ALL tracts (not each cluster's own
# min/max) before computing Gini - using a per-cluster rescale would apply a
# DIFFERENT linear transform to every cluster and make the resulting Gini
# values incomparable across clusters; a single global rescale per domain
# preserves relative dispersion, so a more spread-out cluster still reads as
# more unequal than a tighter one.
suppressMessages(library(DescTools))

domain_lab_wrapped <- vapply(pc_labels[pcs], function(x)
  paste(strwrap(x, width = 12), collapse = "\n"), character(1))

gini_input <- d %>%
  dplyr::select(cluster, dplyr::all_of(pcs)) %>%
  dplyr::mutate(dplyr::across(
    dplyr::all_of(pcs),
    ~ (.x - min(.x, na.rm = TRUE)) / (max(.x, na.rm = TRUE) - min(.x, na.rm = TRUE))
  )) %>%
  pivot_longer(dplyr::all_of(pcs), names_to = "domain", values_to = "value01") %>%
  dplyr::filter(!is.na(value01)) %>%
  dplyr::mutate(domain_lab = factor(domain_lab_wrapped[domain], levels = domain_lab_wrapped[pcs]))

gini_profile <- gini_input %>%
  dplyr::group_by(cluster, domain_lab) %>%
  dplyr::summarise(gini = DescTools::Gini(value01, na.rm = TRUE), .groups = "drop") %>%
  dplyr::left_join(cluster_info %>% dplyr::select(cluster, n_tract), by = "cluster")

cat("\nwithin-cluster Gini coefficient by domain:\n")
print(as.data.frame(gini_profile %>% tidyr::pivot_wider(
  id_cols = cluster, names_from = domain_lab, values_from = gini)), row.names = FALSE)
write_csv(gini_profile, file.path(season_dir, "cluster_domain_gini.csv"))

radar_gini_dat <- gini_profile %>%
  dplyr::mutate(
    facet_lab = factor(sprintf("Cluster %s\n(n = %s tracts)", cluster,
                               format(n_tract, big.mark = ",")),
                       levels = sprintf("Cluster %s\n(n = %s tracts)",
                                        cluster_info$cluster,
                                        format(cluster_info$n_tract, big.mark = ",")))
  )

# Reference contour at the POPULATION Gini per domain (all tracts pooled,
# ignoring cluster) - lets a reader see at a glance whether a cluster is MORE
# or LESS internally unequal than CONUS as a whole for that domain. Traced on
# the axes rather than geom_hline, since geom_hline would be drawn as a
# straight chord under the is_linear override in coord_radar().
pop_gini <- gini_input %>%
  dplyr::group_by(domain_lab) %>%
  dplyr::summarise(gini = DescTools::Gini(value01, na.rm = TRUE), .groups = "drop")

zero_ref_gini <- radar_gini_dat %>%
  dplyr::distinct(facet_lab, domain_lab) %>%
  dplyr::left_join(pop_gini, by = "domain_lab")

rad_lim_gini <- c(0, max(radar_gini_dat$gini, zero_ref_gini$gini, na.rm = TRUE) * 1.15)

p_gini_radar <- ggplot(radar_gini_dat, aes(x = domain_lab, y = gini, group = cluster)) +
  geom_polygon(data = zero_ref_gini, aes(x = domain_lab, y = gini, group = facet_lab),
              inherit.aes = FALSE, fill = NA, color = "grey55",
              linewidth = 0.25, linetype = "dashed") +
  geom_polygon(aes(fill = factor(cluster), color = factor(cluster)),
              alpha = 0.4, linewidth = 0.45) +
  geom_point(aes(color = factor(cluster)), size = 0.9) +
  coord_radar() +
  facet_wrap(~ facet_lab, nrow = 2) +
  scale_fill_manual(values = severity_pal, guide = "none") +
  scale_color_manual(values = severity_pal, guide = "none") +
  scale_y_continuous(limits = rad_lim_gini, breaks = scales::pretty_breaks(n = 3)) +
  labs(
    title    = "Within-cluster inequality (Gini coefficient) by domain and heatwave exposure cluster",
    subtitle = paste0(
      "Gini coefficient of each domain's score WITHIN each cluster (domains min-max rescaled to ",
      "[0,1] on their global range first, so values are comparable across clusters).\n",
      "Dashed contour = population Gini per domain (all tracts pooled)."
    ),
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 9) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    axis.text.x      = element_text(size = 6.2, face = "bold", lineheight = 0.85),
    axis.text.y      = element_text(size = 5.5),
    panel.grid.major = element_line(linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.spacing    = unit(0.4, "cm"),
    strip.text       = element_text(face = "bold", size = 7.5, lineheight = 1.0),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 8.3, color = "grey30")
  )

if (interactive()) print(p_gini_radar)

ggsave(
  filename = file.path(fig_dir, "domain_gini_radar_by_cluster.png"),
  plot     = p_gini_radar,
  width    = 13,
  height   = 8,
  dpi      = 300,
  device   = png_dev
)


#---- 7f. descriptive figures on heat vulnerability variables (violin + boxplot) #----
# Same severity_pal color scheme as 7b-7e, one panel per variable, no
# significance testing - descriptive only. Covers the 17 raw variables that
# feed heat_vulnerability.R's PCA (its screened pca_vars), in ORIGINAL units -
# not the z-scored component space used in 7b/7c/7e.
#
# TRACT level, not county level (unlike 7d): these ARE tract-level
# characteristics (one value per tract), so this uses `d` (one row per tract,
# already carrying `cluster` from the county join in section 2) rather than
# `county`.
#
# Data prep (vuln_vars ... vuln_dom) is built ONCE here and shared by both
# figures below - violin (independent y-axis, panel width scales with row
# length) and boxplot (y-axis zoomed to the whisker range, every panel the
# same size) - since both plots are two different views of the exact same
# 17-variable, tract-level dataset; there is no reason to rebuild it twice.
vuln_vars <- c(
  "pop_density", "pct_age65", "median_income",
  "poverty", "energy_burden", "electricity_burden", "AC_central", "AC_none",
  "tree_canopy", "ndvi",
  "dev_medium_pct", "dev_high_pct", "building_coverage", "impervious",
  "pct_disability", "pct_nonwhite", "pct_living_alone"
)

# Units taken from how each variable is actually built in heat_vulnerability.R
# (energy/electricity burden are raw spend/income ratios, never multiplied by
# 100 there; AC_central/AC_none, tree_canopy, dev_*_pct, building_coverage,
# impervious, and the three ACS shares are all already 0-100 percentages;
# ndvi is a unitless index).
vuln_units <- c(
  pop_density         = paste0("people/km", intToUtf8(0x00B2)),
  pct_age65           = "%",
  median_income       = "$",
  poverty             = "%",
  energy_burden       = "ratio",
  electricity_burden  = "ratio",
  AC_central          = "%",
  AC_none             = "%",
  tree_canopy         = "%",
  ndvi                = "index",
  dev_medium_pct      = "%",
  dev_high_pct        = "%",
  building_coverage   = "%",
  impervious          = "%",
  pct_disability      = "%",
  pct_nonwhite        = "%",
  pct_living_alone    = "%"
)

vuln_var_labs <- c(
  pop_density         = "Population density",
  pct_age65           = "Age 65+",
  median_income       = "Median household income",
  poverty             = "Poverty rate",
  energy_burden       = "Energy burden",
  electricity_burden  = "Electricity burden",
  AC_central          = "Central AC access",
  AC_none             = "No AC access",
  tree_canopy         = "Tree canopy cover",
  ndvi                = "NDVI (greenness)",
  dev_medium_pct      = "Medium-intensity development",
  dev_high_pct        = "High-intensity development",
  building_coverage   = "Building coverage",
  impervious          = "Impervious surface",
  pct_disability      = "Population with a disability",
  pct_nonwhite        = "Non-white population",
  pct_living_alone    = "Living alone"
)

vuln_metric_labs <- sprintf("%s (%s)", vuln_var_labs[vuln_vars], vuln_units[vuln_vars])

# Cache: only GEOID + these 17 raw columns are needed from tract_master, so
# cache them rather than re-reading the 542 MB file. Independent of section
# 8's pop-only cache below (tract_pop.csv) - this section runs first and needs
# a different column set, so there is no shared state to keep in sync.
vuln_raw_cache <- file.path(season_dir, "tract_vuln_raw_vars.csv")
if (!file.exists(vuln_raw_cache)) {
  message("caching raw vulnerability variables from tract_master_2020.rds ...")
  suppressMessages(library(sf))
  tm <- readRDS(file.path(data_dir, "tract_master_2020.rds"))
  tm %>% sf::st_drop_geometry() %>%
    dplyr::select(GEOID, dplyr::all_of(vuln_vars)) %>%
    write_csv(vuln_raw_cache)
  rm(tm); invisible(gc())
}

vuln_raw <- read_csv(vuln_raw_cache,
                     col_types = cols(GEOID = col_character(), .default = col_double()))

violin_dat_vuln <- d %>%
  dplyr::select(GEOID, cluster) %>%
  dplyr::inner_join(vuln_raw, by = "GEOID") %>%
  dplyr::filter(!is.na(cluster)) %>%
  tidyr::pivot_longer(dplyr::all_of(vuln_vars), names_to = "metric", values_to = "value") %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::mutate(
    metric_lab  = factor(sprintf("%s (%s)", vuln_var_labs[metric], vuln_units[metric]),
                         levels = vuln_metric_labs),
    cluster_ord = factor(cluster, levels = seq_len(nrow(cluster_info)))
  )

cat("\ntracts by cluster (vulnerability-variable violins):\n")
print(table(violin_dat_vuln$cluster_ord[violin_dat_vuln$metric == vuln_vars[1]]))

# --- group variables by the PCA domain they load most strongly on ----------
# One row per domain (5 rows), so a reader can see at a glance which
# vulnerability domain each raw variable stands in for. The assignment is
# read straight off pca_loadings_contributions.csv (`L`, `pc_labels`, `pcs`,
# already built in section 2) - the domain with the largest |loading| for
# that variable - rather than hardcoded, so it stays correct if the PCA
# screening set or its loadings change.
stopifnot(all(vuln_vars %in% rownames(L)))
vuln_dom <- colnames(L)[apply(abs(L[vuln_vars, , drop = FALSE]), 1, which.max)]
names(vuln_dom) <- vuln_vars

cat("\nvariable -> dominant PCA domain:\n")
print(data.frame(variable = vuln_vars, domain = unname(pc_labels[vuln_dom])),
      row.names = FALSE)

# Widest domain (5, Socioeconomic disadvantage) - shared by both figures
# below: both the violin and boxplot grids need to know how many
# variable-columns to lay out, and both pad shorter rows with plot_spacer()
# up to this width so every panel ends up the same size regardless of row.
n_cols_row <- max(table(vuln_dom))

suppressMessages(library(patchwork))

# A row's domain name as its own narrow label panel in column 1, rather than
# a nested plot_annotation(title=...) per row - confirmed by rendering that a
# title on a row's own wrap_plots() does not survive being nested inside the
# outer wrap_plots() below, so it never appeared in the saved figure. An
# explicit ggplot (blank canvas + one rotated text annotation) is a real
# panel, sized like any other, and reliably renders through the same nesting.
# Shared by both 7f.1 and 7f.2.
vuln_row_label <- function(domain_name) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = domain_name, angle = 90,
            fontface = "bold", size = 3.3, color = "grey20", hjust = 0.5, vjust = 0.5) +
    xlim(-1, 1) + ylim(-1, 1) +
    theme_void()
}

# --- 7f.1: violin version - independent y-axis per panel, but EVERY panel ---
# the same size regardless of row length (matches 7f.2's boxplot grid below,
# not facet_wrap - facet_wrap(nrow=1) per row would stretch a 2-variable row
# like "Low AC access" to fill the same width as a 5-variable row, which is
# exactly the stretching this section is avoiding). Same row-label +
# plot_spacer()-padded patchwork grid as 7f.2; widths = c(label column, then
# one column per variable slot) is passed identically to every row's
# wrap_plots() call, so the label column and each variable column are the
# same width in every row, not just within a row.
vuln_violin_plot <- function(var) {
  dat_v <- violin_dat_vuln %>% dplyr::filter(metric == var)
  ggplot(dat_v, aes(x = cluster_ord, y = value, fill = cluster_ord)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.15, alpha = 0.9) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white",
                alpha = 0.7, linewidth = 0.25) +
    scale_fill_manual(values = severity_pal, guide = "none") +
    labs(title = sprintf("%s (%s)", vuln_var_labs[var], vuln_units[var]), x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title         = element_text(face = "bold", size = 7.5, hjust = 0.5),
      axis.text.x        = element_text(size = 6.5),
      axis.text.y        = element_text(size = 6.5)
    )
}

vuln_row_plots <- lapply(pcs, function(p) {
  vars_p <- vuln_vars[vuln_dom[vuln_vars] == p]
  panels <- c(list(vuln_row_label(pc_labels[p])),
             lapply(vars_p, vuln_violin_plot),
             rep(list(patchwork::plot_spacer()), n_cols_row - length(vars_p)))
  patchwork::wrap_plots(panels, ncol = n_cols_row + 1, widths = c(0.3, rep(1, n_cols_row)))
})

p_violin_vuln <- patchwork::wrap_plots(vuln_row_plots, ncol = 1) +
  patchwork::plot_annotation(
    title    = "Tract heat-vulnerability variables by county exposure cluster",
    subtitle = paste0(
      "One row per PCA domain, one panel per variable (original units), all panels the same ",
      "size. Independent y-axis per panel. Cluster 1 = mildest ", arrow_sym,
      " 7 = most severe heatwave exposure."
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "grey30")
    )
  )

if (interactive()) print(p_violin_vuln)

ggsave(
  filename = file.path(fig_dir, "vuln_variable_violins_by_cluster.png"),
  plot     = p_violin_vuln,
  width    = 14,
  height   = 13,
  dpi      = 300,
  device   = png_dev
)

violin_summary_vuln <- violin_dat_vuln %>%
  dplyr::group_by(metric, cluster_ord) %>%
  dplyr::summarise(
    n      = dplyr::n(),
    mean   = mean(value),
    median = median(value),
    sd     = sd(value),
    iqr    = IQR(value),
    min    = min(value),
    max    = max(value),
    .groups = "drop"
  ) %>%
  dplyr::arrange(metric, cluster_ord)

write_csv(violin_summary_vuln,
          file.path(season_dir, "vuln_variable_violin_summary.csv"))


# --- 7f.2: boxplot version - every panel the same size, y-axis reflects ----
# the DISTRIBUTION rather than the full data range.
#
#   1) A violin's kernel density (and a boxplot's default range) both extend
#      to the actual min/max, so a handful of extreme tracts (e.g. one very
#      dense tract) stretch the axis and flatten everything else - confirmed
#      this is real, not just a display quirk: outlier.shape = NA alone does
#      NOT shrink ggplot2's computed y-range here (tested directly with
#      ggplot_build() - the outlier is still included in panel_scales_y even
#      when its point isn't drawn). Instead, each panel's y-axis is set via
#      coord_cartesian(ylim = ...) to the Tukey whisker range (Q1-1.5*IQR to
#      Q3+1.5*IQR, from boxplot.stats()) plus 8% padding - coord_cartesian
#      zooms the viewport without dropping data, so the box/whiskers/outlier
#      points themselves are still computed from the full distribution.
#
#   2) 7f.1's facet_wrap(nrow=1) per row stretches panels to fill the row
#      width, so the 2-variable "Low AC access" row ends up with much wider
#      panels than the 5-variable rows. Here every row is its own patchwork
#      grid with the SAME ncol (n_cols_row, the widest domain) and
#      plot_spacer() padding for shorter rows, so a panel's width is always
#      (figure width) / n_cols_row no matter which row it is in.
vuln_ylim <- function(x) {
  s   <- boxplot.stats(x)$stats   # Tukey 5-number: lo whisker..hi whisker
  pad <- 0.08 * diff(range(s))
  if (pad == 0) pad <- 1e-6 * max(abs(s), 1)  # guard against a degenerate/constant panel
  c(s[1] - pad, s[5] + pad)
}

vuln_box_plot <- function(var) {
  dat_v <- violin_dat_vuln %>% dplyr::filter(metric == var)
  ggplot(dat_v, aes(x = cluster_ord, y = value, fill = cluster_ord)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.3) +
    coord_cartesian(ylim = vuln_ylim(dat_v$value)) +
    scale_fill_manual(values = severity_pal, guide = "none") +
    labs(title = sprintf("%s (%s)", vuln_var_labs[var], vuln_units[var]), x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title         = element_text(face = "bold", size = 7.5, hjust = 0.5),
      axis.text.x        = element_text(size = 6.5),
      axis.text.y        = element_text(size = 6.5)
    )
}

# Same row-label + plot_spacer()-padded grid as 7f.1 (vuln_row_label,
# n_cols_row both built above), so panels are the same size in both figures.
vuln_row_plots_box <- lapply(pcs, function(p) {
  vars_p <- vuln_vars[vuln_dom[vuln_vars] == p]
  panels <- c(list(vuln_row_label(pc_labels[p])),
             lapply(vars_p, vuln_box_plot),
             rep(list(patchwork::plot_spacer()), n_cols_row - length(vars_p)))
  patchwork::wrap_plots(panels, ncol = n_cols_row + 1, widths = c(0.3, rep(1, n_cols_row)))
})

p_box_vuln <- patchwork::wrap_plots(vuln_row_plots_box, ncol = 1) +
  patchwork::plot_annotation(
    title    = "Tract heat-vulnerability variables by county exposure cluster (boxplots)",
    subtitle = paste0(
      "One row per PCA domain, one panel per variable (original units), all panels the same ",
      "size. Boxplot only; y-axis zoomed to the Tukey whisker range (Q1-1.5", intToUtf8(0x00D7),
      "IQR to Q3+1.5", intToUtf8(0x00D7), "IQR) rather than the full data range. ",
      "Cluster 1 = mildest ", arrow_sym, " 7 = most severe heatwave exposure."
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "grey30")
    )
  )

if (interactive()) print(p_box_vuln)

ggsave(
  filename = file.path(fig_dir, "vuln_variable_boxplots_by_cluster_6.png"),
  plot     = p_box_vuln,
  width    = 14,
  height   = 13,
  dpi      = 300,
  device   = png_dev
)




#===========================================================================
# 8. VULNERABILITY INDEX FROM THE SCREENED COMPONENT SCORES
#===========================================================================
# Previously this block re-ran its own 19-variable PCA from tract_master. That
# duplicated - and disagreed with - heat_vulnerability.R section 9, which now
# exports screened per-tract scores. This block consumes those instead: no
# second PCA, no 542 MB read, one definition of the components.
#
# What still has to happen here is ORIENTATION. heat_vulnerability.R
# deliberately leaves the components unoriented, because PCA structure and
# k-means are both sign-invariant. Summing them into a cumulative index is not:
# a component pointing toward LOWER vulnerability would subtract from the total.

# sf is needed by sections 10-11 regardless of whether the population cache has
# to be rebuilt below - loading it only inside that branch would work on the
# first run and fail on every later one.
suppressMessages(library(sf))

# --- verify the scores are already in vulnerability space ------------------
# heat_vulnerability.R section 9 now negates the protective variables BEFORE the
# PCA and sign-checks each component after it, so the exported scores already
# point toward higher vulnerability. This block therefore VERIFIES that rather
# than re-deriving it: applying a variable-level orientation here as well would
# flip the protective variables a second time and undo the correction upstream.
comp_sum <- colSums(L)

cat("\n=== 8. Component orientation (verifying upstream) ===\n")
print(data.frame(
  component = pcs,
  label     = unname(pc_labels[pcs]),
  sum_load  = round(comp_sum, 3),
  status    = ifelse(comp_sum >= 0, "vulnerability-positive", "WRONG SIGN"),
  row.names = NULL
))

if (any(comp_sum < 0)) {
  stop("component(s) ", paste(pcs[comp_sum < 0], collapse = ", "),
       " are not vulnerability-positive. The exported scores are not oriented - ",
       "re-run heat_vulnerability.R section 9 before continuing.")
}

L_vuln <- L   # already oriented upstream; kept as a name for the checks below

# Bipolar check: if a component's SALIENT loadings disagree in sign after
# reorientation, it is a contrast between two things rather than a gradient of
# one, and no single direction means "more vulnerable". Summing such a
# component into a cumulative index is ambiguous regardless of which way it is
# flipped - flag it rather than let the arithmetic hide it.
SALIENT <- 0.40
bipolar <- vapply(pcs, function(p) {
  s <- L_vuln[abs(L_vuln[, p]) >= SALIENT, p]
  length(s) > 1 && any(s > 0) && any(s < 0)
}, logical(1))

if (any(bipolar)) {
  cat("\n!! bipolar component(s):", paste(pcs[bipolar], collapse = ", "), "\n")
  for (p in pcs[bipolar]) {
    s <- L_vuln[abs(L_vuln[, p]) >= SALIENT, p]
    cat(sprintf("   %s (%s): %s\n", p, pc_labels[p],
                paste(sprintf("%s %+.2f", names(s), s), collapse = ", ")))
  }
  cat("   Both poles carry vulnerability, so a signed sum partially cancels\n")
  cat("   them. Treated as-is below; consider reporting these separately.\n")
}

# --- build the index -------------------------------------------------------
scores_oriented <- as.matrix(pca[pcs])   # already oriented upstream
colnames(scores_oriented) <- pcs

vuln_tract <- pca %>%
  dplyr::select(GEOID, StCoFIPS) %>%
  bind_cols(as.data.frame(scores_oriented)) %>%
  dplyr::mutate(hvi = as.numeric(scale(rowSums(scores_oriented))))

# Variance-weighted alternative, kept as a sensitivity check: equal weighting
# lets a component explaining 9.5% count as much as one explaining 30.4%.
w <- pca_var$prop_var_pct / sum(pca_var$prop_var_pct)
vuln_tract$hvi_wvar <- as.numeric(scale(as.numeric(scores_oriented %*% w)))
cat(sprintf("\nequal-weight vs variance-weighted index: r = %+.3f\n",
            cor(vuln_tract$hvi, vuln_tract$hvi_wvar)))

# --- population, for weighting the county aggregation ----------------------
# Only GEOID and pop are needed from tract_master, so cache those two columns
# rather than re-reading 542 MB on every run.
pop_cache <- file.path(season_dir, "tract_pop.csv")
if (!file.exists(pop_cache)) {
  message("caching tract population from tract_master_2020.rds ...")
  suppressMessages(library(sf))
  tm <- readRDS(file.path(data_dir, "tract_master_2020.rds"))
  tm %>% sf::st_drop_geometry() %>%
    dplyr::select(GEOID, pop) %>%
    write_csv(pop_cache)
  rm(tm); invisible(gc())
}

tract_pop <- read_csv(pop_cache, col_types = cols(GEOID = col_character(),
                                                  pop = col_double()))

vuln_tract <- vuln_tract %>% left_join(tract_pop, by = "GEOID")
n_nopop <- sum(is.na(vuln_tract$pop))
if (n_nopop > 0) message("note: ", n_nopop, " tracts have no population value")

write_csv(vuln_tract, file.path(season_dir, "tract_vulnerability_index.csv"))

message("vulnerability index: ", format(nrow(vuln_tract), big.mark = ","),
        " tracts, ", length(pcs), " components")

#===========================================================================
# 9. COUNTY AGGREGATION - the two scalars for the bivariate map and LISA
#===========================================================================
# The LISA must run at ONE scale. Hazard is constant within a county, so pushing
# it down to tracts would give 83k observations carrying only ~3.1k independent
# hazard values - pseudo-replication that badly inflates the permutation test.
# Vulnerability is therefore aggregated UP, population-weighted.

county_vuln <- vuln_tract %>%
  dplyr::mutate(StCoFIPS = substr(GEOID, 1, 5)) %>%
  dplyr::group_by(StCoFIPS) %>%
  dplyr::summarise(
    vuln       = weighted.mean(hvi, pop, na.rm = TRUE),
    vuln_iqr   = IQR(hvi, na.rm = TRUE),      # within-county spread, kept for
    n_tract    = dplyr::n(),                  # the within/between argument
    pop        = sum(pop, na.rm = TRUE),
    .groups    = "drop"
  )

bv_dat <- county %>%
  dplyr::select(StCoFIPS, hazard, cluster, dplyr::all_of(hazard_vars)) %>%
  inner_join(county_vuln, by = "StCoFIPS") %>%
  dplyr::filter(is.finite(hazard), is.finite(vuln))

message(sprintf("counties with both hazard and vulnerability: %s of %s",
                format(nrow(bv_dat), big.mark = ","),
                format(nrow(county),  big.mark = ",")))

cat("\ncorrelation between county hazard and population-weighted vulnerability:\n")
cat(sprintf("  Pearson  r   = %+.3f\n", cor(bv_dat$hazard, bv_dat$vuln)))
cat(sprintf("  Spearman rho = %+.3f\n",
            cor(bv_dat$hazard, bv_dat$vuln, method = "spearman")))


#===========================================================================
# 10. BIVARIATE MAP - hazard tercile x vulnerability tercile
#===========================================================================
options(tigris_use_cache = TRUE)
suppressMessages(library(tigris))

counties_sf <- suppressMessages(
  counties(cb = TRUE, year = 2020, class = "sf", progress_bar = FALSE)
) %>%
  dplyr::filter(!STATEFP %in% c("02", "15", "60", "66", "69", "72", "78")) %>%
  st_transform(5070)

terc <- function(x) cut(x, breaks = quantile(x, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                        labels = 1:3, include.lowest = TRUE)

bv_dat <- bv_dat %>%
  dplyr::mutate(
    haz_t  = terc(hazard),
    vuln_t = terc(vuln),
    biv    = paste(haz_t, vuln_t, sep = "-")
  )

# Stevens-style bivariate palette, keyed "hazard-vulnerability".
# Pure hazard -> red, pure vulnerability -> teal, both high -> dark.
# "no data" must not collide with the palette's lightest cell (#e8e8e8), or
# unscored areas read as low-hazard / low-vulnerability.
no_data_col <- "grey80"

biv_pal <- c(
  "1-1" = "#e8e8e8", "2-1" = "#e4acac", "3-1" = "#c85a5a",
  "1-2" = "#b0d5df", "2-2" = "#ad9ea5", "3-2" = "#985356",
  "1-3" = "#64acbe", "2-3" = "#627f8c", "3-3" = "#574249"
)

bv_sf <- counties_sf %>%
  left_join(bv_dat, by = c("GEOID" = "StCoFIPS"))

# Population living in each bivariate class - the headline numbers.
biv_summary <- bv_dat %>%
  dplyr::group_by(haz_t, vuln_t, biv) %>%
  dplyr::summarise(n_county = dplyr::n(),
                   pop_mil  = sum(pop, na.rm = TRUE) / 1e6,
                   .groups  = "drop") %>%
  dplyr::arrange(desc(pop_mil))

cat("\n=== population by hazard x vulnerability tercile ===\n")
print(as.data.frame(biv_summary), row.names = FALSE)
write_csv(biv_summary, file.path(season_dir, "bivariate_summary.csv"))

p_biv <- ggplot(bv_sf) +
  geom_sf(aes(fill = biv), color = NA) +
  scale_fill_manual(values = biv_pal, na.value = no_data_col, guide = "none") +
  labs(
    title    = "Heat hazard and population vulnerability across CONUS counties",
    subtitle = paste0("Terciles of county heatwave hazard (frequency + intensity + duration + span) ",
                      "crossed with\npopulation-weighted heat vulnerability index. ",
                      format(nrow(bv_dat), big.mark = ","), " counties.")
  ) +
  theme_void(base_size = 11) +
  # theme_void() leaves the canvas transparent, which renders as black in most
  # viewers and hides the black title/legend text. Force an opaque white canvas.
  theme(plot.background = element_rect(fill = "white", color = NA),
        plot.title      = element_text(face = "bold", size = 13),
        plot.subtitle   = element_text(size = 9, color = "grey30"))

# 3x3 legend, drawn as its own plot and inset into the map.
leg_df <- expand_grid(haz_t = 1:3, vuln_t = 1:3) %>%
  dplyr::mutate(biv = paste(haz_t, vuln_t, sep = "-"))

p_leg <- ggplot(leg_df, aes(x = haz_t, y = vuln_t, fill = biv)) +
  geom_tile() +
  scale_fill_manual(values = biv_pal, guide = "none") +
  labs(x = paste0("Hazard ", arrow_sym),
       y = paste0("Vulnerability ", arrow_sym)) +
  coord_fixed() +
  theme_minimal(base_size = 8) +
  theme(axis.text    = element_blank(),
        axis.title   = element_text(size = 7.5),
        panel.grid   = element_blank(),
        plot.background = element_rect(fill = "white", color = NA))

suppressMessages(library(patchwork))
p_biv_full <- p_biv + inset_element(p_leg, left = 0.02, bottom = 0.02,
                                    right = 0.20, top = 0.30)

if (interactive()) print(p_biv_full)

ggsave(file.path(fig_dir, "bivariate_hazard_vulnerability_map.png"),
       p_biv_full, width = 11, height = 7.2, dpi = 600, device = png_dev)


#===========================================================================
# 11. BIVARIATE LISA - is the mismatch itself spatially clustered?
#===========================================================================
# x = own vulnerability, y = neighbours' hazard. Note I(x,y) != I(y,x); this
# orientation answers "are vulnerable counties surrounded by hazard?", so the
# off-diagonal quadrants read as the two mismatch types.
suppressMessages(library(spdep))

lisa_sf <- bv_sf %>% dplyr::filter(!is.na(vuln), !is.na(hazard))

nb <- poly2nb(lisa_sf, queen = TRUE)

# Island counties have no queen contiguity, and moran_bv / lee.mc require every
# unit to have at least one neighbour. Rather than silently dropping them, link
# each isolate to its single nearest county (symmetrically, so the weights stay
# well formed). Patching only the isolates leaves every other county's
# neighbour set exactly as queen contiguity defined it.
iso <- which(card(nb) == 0)
if (length(iso) > 0) {
  message("patching ", length(iso), " isolated county/counties to nearest neighbour: ",
          paste(lisa_sf$NAME[iso], collapse = ", "))
  cent <- st_coordinates(suppressWarnings(st_point_on_surface(st_geometry(lisa_sf))))
  for (i in iso) {
    dd <- sqrt((cent[, 1] - cent[i, 1])^2 + (cent[, 2] - cent[i, 2])^2)
    dd[i] <- Inf
    j <- as.integer(which.min(dd))
    nb[[i]] <- j
    nb[[j]] <- sort(unique(c(nb[[j]][nb[[j]] != 0L], as.integer(i))))
  }
  attr(nb, "sym") <- is.symmetric.nb(nb)
}
stopifnot(all(card(nb) > 0))

lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

set.seed(123)

# Global measures
gm <- moran_bv(lisa_sf$vuln, lisa_sf$hazard, lw, nsim = 999)
cat("\n=== global bivariate Moran's I (vulnerability vs lagged hazard) ===\n")
cat(sprintf("  I = %+.4f   pseudo-p = %.4f\n", gm$t0,
            (sum(abs(gm$t) >= abs(gm$t0)) + 1) / (length(gm$t) + 1)))

# Lee's L - integrates Pearson r with spatial smoothing. Reported alongside
# because bivariate Moran's I is known to largely reflect the plain in-situ
# correlation between x and y rather than genuine spatial association.
lee_res <- lee.mc(lisa_sf$vuln, lisa_sf$hazard, lw, nsim = 999,
                  zero.policy = TRUE, alternative = "two.sided")
cat(sprintf("  Lee's L = %+.4f   pseudo-p = %.4f\n",
            lee_res$statistic, lee_res$p.value))

# Local
bv_lisa <- localmoran_bv(lisa_sf$vuln, lisa_sf$hazard, lw, nsim = 999)

# spdep returns both an analytical and a permutation p-value. Use the folded
# permutation one ("Pr(folded) Sim") - that is what nsim was spent on, and it
# is the convention GeoDa uses for LISA significance.
i_col <- grep("^(Ibvi|Ii)$", colnames(bv_lisa), value = TRUE)[1]
p_col <- if ("Pr(folded) Sim" %in% colnames(bv_lisa)) {
  "Pr(folded) Sim"
} else {
  grep("Sim$", colnames(bv_lisa), value = TRUE)[1]
}
message("using statistic column: ", i_col, " | p-value column: ", p_col)

zx      <- as.numeric(scale(lisa_sf$vuln))
lag_zy  <- lag.listw(lw, as.numeric(scale(lisa_sf$hazard)), zero.policy = TRUE)

lisa_sf <- lisa_sf %>%
  dplyr::mutate(
    Ii     = bv_lisa[, i_col],
    p_sim  = bv_lisa[, p_col],
    quad   = dplyr::case_when(
      zx > 0 & lag_zy > 0 ~ "High vuln / High hazard",
      zx < 0 & lag_zy < 0 ~ "Low vuln / Low hazard",
      zx > 0 & lag_zy < 0 ~ "High vuln / Low hazard",
      zx < 0 & lag_zy > 0 ~ "Low vuln / High hazard",
      TRUE                ~ NA_character_
    ),
    quad_sig = factor(
      dplyr::if_else(!is.na(p_sim) & p_sim < 0.05, quad, "Not significant"),
      levels = c("High vuln / High hazard", "Low vuln / Low hazard",
                 "High vuln / Low hazard",  "Low vuln / High hazard",
                 "Not significant")
    )
  )

lisa_cols <- c(
  "High vuln / High hazard" = "#b2182b",   # both high - priority
  "Low vuln / Low hazard"   = "#2166ac",   # both low
  "High vuln / Low hazard"  = "#ef8a62",   # MISMATCH: latent risk
  "Low vuln / High hazard"  = "#67a9cf",   # MISMATCH: hazard, low vulnerability
  "Not significant"         = "grey90"
)

cat("\n=== bivariate LISA quadrants (p < 0.05) ===\n")
lisa_tab <- lisa_sf %>%
  st_drop_geometry() %>%
  dplyr::group_by(quad_sig) %>%
  dplyr::summarise(n_county = dplyr::n(),
                   pop_mil  = sum(pop, na.rm = TRUE) / 1e6,
                   .groups  = "drop")
print(as.data.frame(lisa_tab), row.names = FALSE)
write_csv(lisa_tab, file.path(season_dir, "lisa_quadrant_summary.csv"))

write_csv(
  lisa_sf %>% st_drop_geometry() %>%
    dplyr::select(GEOID, vuln, hazard, vuln_iqr, pop, Ii, p_sim, quad, quad_sig),
  file.path(season_dir, "lisa_bivariate_counties.csv")
)

p_lisa <- ggplot() +
  # base layer so counties without data read as light grey rather than as holes
  geom_sf(data = counties_sf, fill = "grey95", color = NA) +
  geom_sf(data = lisa_sf, aes(fill = quad_sig), color = NA) +
  scale_fill_manual(values = lisa_cols, name = NULL, drop = FALSE) +
  labs(
    title    = "Spatial mismatch between heat vulnerability and heatwave hazard",
    subtitle = paste0("Bivariate local Moran's I: county vulnerability vs. ",
                      "spatially lagged hazard.\nQueen contiguity, 999 permutations, ",
                      "p < 0.05. ", format(nrow(lisa_sf), big.mark = ","), " counties."),
    caption  = paste0("Orange = vulnerable populations with comparatively low surrounding hazard. ",
                      "Blue = high surrounding hazard with comparatively low vulnerability.\n",
                      "These two classes are the mismatch; red and dark blue are the aligned cases.")
  ) +
  theme_void(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        plot.title      = element_text(face = "bold", size = 13),
        plot.subtitle   = element_text(size = 9,   color = "grey30"),
        plot.caption    = element_text(size = 7.6, color = "grey40", hjust = 0),
        legend.position = "right")

if (interactive()) print(p_lisa)

ggsave(file.path(fig_dir, "lisa_bivariate_mismatch_map.png"),
       p_lisa, width = 11, height = 7.2, dpi = 600, device = png_dev)

#===========================================================================
# 12. MIXED-MODEL TEST FOR THE SECTION 6 HEATMAP
#===========================================================================
# The heatmap reports eta^2 and no p-values, because a tract-level test would
# treat 69k tracts as independent when the exposure cluster is assigned at the
# COUNTY level and is constant within a county. That inflates significance.
#
# Correct specification: cluster is a level-2 (between-county) predictor, so it
# must be tested against between-county variance with county as a random
# intercept:
#
#     score ~ factor(cluster) + (1 | StCoFIPS)
#
# lmerTest is not installed here, so the omnibus test is a likelihood-ratio
# test of the cluster term against the intercept-only model, both refit with
# ML (REML likelihoods are not comparable across different fixed effects).
# With 3,106 level-2 units the LRT's small-sample anticonservatism is moot.
#
# Also reported:
#   ICC   - share of tract variance sitting between counties (null model)
#   R2m   - marginal R2, variance explained by the cluster fixed effect
#           (Nakagawa & Schielzeth), the effect size to quote
#   p_naive - tract-level one-way ANOVA p, shown only to expose the inflation

suppressMessages(library(lme4))

mm_results <- lapply(pcs, function(p) {

  message("  fitting mixed model: ", p)

  dm <- d %>%
    dplyr::select(score = dplyr::all_of(p), cluster, StCoFIPS) %>%
    dplyr::filter(!is.na(score)) %>%
    dplyr::mutate(cluster = factor(cluster))

  m0 <- lmer(score ~ 1 + (1 | StCoFIPS), data = dm, REML = FALSE)
  m1 <- lmer(score ~ cluster + (1 | StCoFIPS), data = dm, REML = FALSE)

  lrt <- anova(m0, m1)

  vc0    <- as.data.frame(VarCorr(m0))
  tau0   <- vc0$vcov[vc0$grp == "StCoFIPS"]
  sig0   <- vc0$vcov[vc0$grp == "Residual"]

  vc1    <- as.data.frame(VarCorr(m1))
  tau1   <- vc1$vcov[vc1$grp == "StCoFIPS"]
  sig1   <- vc1$vcov[vc1$grp == "Residual"]
  var_f  <- var(as.vector(model.matrix(m1) %*% fixef(m1)))

  # naive tract-level ANOVA, for contrast only
  p_naive <- anova(aov(score ~ cluster, data = dm))[["Pr(>F)"]][1]

  data.frame(
    domain    = p,
    label     = unname(pc_labels[p]),
    n_tract   = nrow(dm),
    n_county  = dplyr::n_distinct(dm$StCoFIPS),
    ICC       = tau0 / (tau0 + sig0),
    eta2      = eta_sq(dm$score, dm$cluster),
    R2m       = var_f / (var_f + tau1 + sig1),
    chisq     = lrt$Chisq[2],
    df        = lrt$Df[2],
    p_mixed   = lrt$`Pr(>Chisq)`[2],
    p_naive   = p_naive
  )
}) %>% bind_rows() %>%
  dplyr::arrange(desc(R2m))

fmtp <- function(x) ifelse(x < 2e-16, "<2e-16", format.pval(x, digits = 3))

cat("\n=== mixed model: score ~ cluster + (1 | county) ===\n")
print(
  mm_results %>%
    dplyr::transmute(
      domain, label,
      ICC     = sprintf("%.3f", ICC),
      eta2    = sprintf("%.1f%%", 100 * eta2),
      R2m     = sprintf("%.1f%%", 100 * R2m),
      chisq   = sprintf("%.1f", chisq),
      df,
      p_mixed = fmtp(p_mixed),
      p_naive = fmtp(p_naive)
    ),
  row.names = FALSE
)
cat("\nn tracts:", mm_results$n_tract[1],
    "| n counties:", mm_results$n_county[1], "\n")
cat("ICC = between-county share of tract variance (null model).\n")
cat("R2m = marginal R2, the cluster fixed effect. Quote this as the effect size.\n")
cat("p_naive is the tract-level ANOVA and is NOT valid here - shown to expose\n")
cat("how far a test that ignores county clustering overstates the evidence.\n")

write_csv(mm_results, file.path(season_dir, "domain_exposure_mixedmodel.csv"))


#===========================================================================
# 13. TRACT-LEVEL BIVARIATE MAP
#===========================================================================
# Same construction and palette as section 10, at tract resolution. County
# hazard is assigned down to its tracts, so the hazard axis has no within-county
# variation - what this map adds over section 10 is the vulnerability axis,
# which is ~63% within-county and therefore invisible at county resolution.
#
# DESCRIPTIVE ONLY. No LISA, no permutation test, no p-values at this scale:
# 69k tracts carry only ~3k independent hazard values, so a permutation null
# would be far too narrow and would manufacture significance almost everywhere.
# Inference stays at county level in section 11.

# Simplified tract geometry, cached. At this extent one pixel is roughly 700 m,
# so a 300 m tolerance is sub-pixel while removing most vertices.
geom_cache <- file.path(season_dir, "tract_geom_simplified.rds")

if (!file.exists(geom_cache)) {
  message("building simplified tract geometry from tract_master_2020.rds ...")
  tm <- readRDS(file.path(data_dir, "tract_master_2020.rds"))
  tract_geom <- tm %>%
    dplyr::select(GEOID) %>%
    st_transform(5070) %>%
    st_simplify(dTolerance = 300, preserveTopology = TRUE)
  saveRDS(tract_geom, geom_cache)
  rm(tm); invisible(gc())
} else {
  tract_geom <- readRDS(geom_cache)
}

# Terciles are computed on the TRACT distribution, so each panel of the 3x3 is
# balanced in tracts. The hazard cut points therefore differ slightly from
# section 10, where counties were the unit - a county with many tracts carries
# proportionally more weight here.
tract_biv <- vuln_tract %>%
  dplyr::select(GEOID, StCoFIPS, hvi, pop) %>%
  left_join(county %>% dplyr::select(StCoFIPS, hazard), by = "StCoFIPS") %>%
  dplyr::filter(is.finite(hazard), is.finite(hvi)) %>%
  dplyr::mutate(
    haz_t  = terc(hazard),
    vuln_t = terc(hvi),
    biv    = paste(haz_t, vuln_t, sep = "-")
  )

message("tract bivariate map: ", format(nrow(tract_biv), big.mark = ","), " tracts")

tract_biv_summary <- tract_biv %>%
  dplyr::group_by(haz_t, vuln_t, biv) %>%
  dplyr::summarise(n_tract = dplyr::n(),
                   pop_mil = sum(pop, na.rm = TRUE) / 1e6,
                   .groups = "drop") %>%
  dplyr::arrange(desc(pop_mil))

cat("\n=== tract population by hazard x vulnerability tercile ===\n")
print(as.data.frame(tract_biv_summary), row.names = FALSE)
write_csv(tract_biv_summary, file.path(season_dir, "bivariate_tract_summary.csv"))

tract_biv_sf <- tract_geom %>%
  left_join(tract_biv, by = "GEOID")

p_biv_tract <- ggplot(tract_biv_sf) +
  geom_sf(aes(fill = biv), color = NA) +
  scale_fill_manual(values = biv_pal, na.value = no_data_col, guide = "none") +
  labs(
    title    = "Heat hazard and population vulnerability across CONUS census tracts",
    # subtitle = paste0("Terciles of county heatwave hazard crossed with tract-level ",
    #                   "heat vulnerability index.\n",
    #                   format(nrow(tract_biv), big.mark = ","),
    #                   " tracts scored; ",
    #                   format(sum(is.na(tract_biv_sf$biv)), big.mark = ","),
    #                   " shown grey (no data, mostly missing AC coverage).\n",
    #                   "Hazard is county-level; the vulnerability axis is what ",
    #                   "this adds over the county map.")
  ) +
  theme_void(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        plot.title      = element_text(face = "bold", size = 13),
        plot.subtitle   = element_text(size = 9, color = "grey30"))

# Reuses the 3x3 legend built for the county map, so both figures read alike.
p_biv_tract_full <- p_biv_tract +
  inset_element(p_leg, left = 0.02, bottom = 0.02, right = 0.20, top = 0.30)

if (interactive()) print(p_biv_tract_full)

ggsave(file.path(fig_dir, "bivariate_tract_hazard_vulnerability_map.png"),
       p_biv_tract_full, width = 11, height = 7.2, dpi = 600, device = png_dev)

message("done - figures written to: ", fig_dir)


#===========================================================================
# 14. METRO-AREA BIVARIATE MAPS - 8 metros across different climate types
#===========================================================================
# Same bivariate classification as section 13 (county hazard assigned down to
# tracts x tract-level vulnerability index, tract_biv_sf) - zoomed into 8
# metro areas chosen to span different US climate types, NOT the project's
# own 7-cluster exposure typology: hot-humid, hot-dry, warm-marine,
# mixed-humid, humid-continental (cold), humid-continental (very cold), and
# marine-cool. 8 panels, 4 columns x 2 rows.

suppressMessages(library(tigris))
options(tigris_use_cache = TRUE)

cbsa_all <- suppressMessages(
  core_based_statistical_areas(cb = TRUE, year = 2020, class = "sf")
) %>% st_transform(5070)

# Exact NAME match (not grepl) against the live tigris CBSA file, verified
# first - a partial match on "Miami" alone also catches an unrelated
# "Miami, OK" micropolitan area.
# ---------------------------------------------------------------- EDIT THIS --
metro_picks <- c(
  "Miami-Fort Lauderdale-Pompano Beach, FL" = "Miami (hot-humid)",
  "Phoenix-Mesa-Chandler, AZ"               = "Phoenix (hot-dry)",
  "Houston-The Woodlands-Sugar Land, TX"    = "Houston (hot-humid)",
  "Los Angeles-Long Beach-Anaheim, CA"      = "Los Angeles (warm-marine)",
  "Atlanta-Sandy Springs-Alpharetta, GA"    = "Atlanta (mixed-humid)",
  "Chicago-Naperville-Elgin, IL-IN-WI"      = "Chicago (cold)",
  "Minneapolis-St. Paul-Bloomington, MN-WI" = "Minneapolis (very cold)",
  "Seattle-Tacoma-Bellevue, WA"             = "Seattle (marine)"
)
# -----------------------------------------------------------------------------
stopifnot(all(names(metro_picks) %in% cbsa_all$NAME))

metro_sf <- cbsa_all %>%
  dplyr::filter(NAME %in% names(metro_picks)) %>%
  dplyr::mutate(metro_lab = factor(metro_picks[NAME], levels = unname(metro_picks)))

cat("\nmetro areas matched:\n")
print(as.data.frame(sf::st_drop_geometry(metro_sf["NAME"])))

# Tract -> metro membership via point-in-polygon on tract CENTROIDS (a tract
# straddling a metro boundary is assigned to whichever side its centroid
# falls on, not duplicated into both) - the membership lookup is then joined
# back onto the FULL tract polygons (tract_biv_sf, already carrying the
# section 13 bivariate class) so the plotted shapes are not clipped to points.
metro_membership <- tract_biv_sf %>%
  sf::st_centroid() %>%
  sf::st_join(metro_sf %>% dplyr::select(metro_lab), join = st_within, left = FALSE) %>%
  sf::st_drop_geometry() %>%
  dplyr::select(GEOID, metro_lab)

metro_tract_sf <- tract_biv_sf %>%
  dplyr::inner_join(metro_membership, by = "GEOID")

cat("\ntracts per metro:\n")
print(table(metro_tract_sf$metro_lab))

# facet_wrap(scales = "free") cannot be combined with coord_sf() - ggplot2
# raises "facet_wrap() can't use free scales with coord_sf()" (confirmed by
# running this). Independent per-panel extents therefore need 8 SEPARATE
# ggplot objects (one per metro, each its own geom_sf(), cropped to that
# metro's own bounding box via coord_sf(xlim=, ylim=)), combined with
# patchwork - the same pattern already used in section 7f for independent
# per-panel scales that facet_wrap couldn't give directly.
metro_plot <- function(lab) {
  dat_m <- metro_tract_sf %>% dplyr::filter(metro_lab == lab)
  bb    <- sf::st_bbox(dat_m)
  ggplot(dat_m) +
    geom_sf(aes(fill = biv), color = NA) +
    scale_fill_manual(values = biv_pal, na.value = no_data_col, guide = "none") +
    coord_sf(xlim = bb[c("xmin", "xmax")], ylim = bb[c("ymin", "ymax")],
            expand = FALSE) +
    labs(title = lab) +
    theme_void(base_size = 9) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.title       = element_text(face = "bold", size = 8.5, hjust = 0.5)
    )
}

metro_plots <- lapply(levels(metro_tract_sf$metro_lab), metro_plot)
p_metro_grid <- patchwork::wrap_plots(metro_plots, ncol = 4, nrow = 2)

# Reuses the same 3x3 legend as sections 10 and 13, so all three figures read
# alike. inset_element() (used there) is designed for inserting into a SINGLE
# plot's panel area - tried directly on this 8-panel patchwork first, and it
# landed on top of the last panel (Seattle) instead of floating independently
# (confirmed by rendering it). Stacking the legend as its own thin row below
# the whole 4x2 grid avoids that overlap entirely, regardless of patchwork's
# inset behavior on multi-panel objects.
p_metro_biv_full <- (p_metro_grid / p_leg) +
  patchwork::plot_layout(heights = c(6, 1)) +
  patchwork::plot_annotation(
    title = "Heat hazard and population vulnerability across 8 U.S. metro areas",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

if (interactive()) print(p_metro_biv_full)

ggsave(file.path(fig_dir, "bivariate_metro_hazard_vulnerability_map.png"),
       p_metro_biv_full, width = 14, height = 6, dpi = 600, device = png_dev)

message("done - metro figure written to: ", fig_dir)

#---- Script Metadata #----
# Title: Tract-level bivariate heat hazard x heat vulnerability map
# Author: Hang Li
# Purpose: Companion to sections 8-11 of domain_exposure_heatmap.R. That script
#          aggregates vulnerability UP to the county so hazard and vulnerability
#          share one scale, which is what the bivariate LISA requires. The cost
#          is that most of the vulnerability signal is thrown away: the index is
#          mostly a WITHIN-county quantity (see the decomposition printed in
#          section 4 below), so population-weighting to the county averages away
#          the contrast between a hot inner-city tract and the suburb next to it.
#
#          This script does the reverse assignment - county hazard pushed DOWN to
#          its tracts - purely so that within-county vulnerability texture is
#          visible on the map. It is DESCRIPTIVE ONLY. See the warning in
#          section 8 before adding any inferential statistic here.
#
# Inputs :  <season dir>/tract_vulnerability_index.csv   (built by domain_exposure_heatmap.R sec. 8)
#           <season dir>/county_clustered.csv
#           <data dir>/tract_master_2020.rds             (tract geometry)
# Outputs:  <fig dir>/bivariate_tract_hazard_vulnerability_map.png
#           <fig dir>/bivariate_tract_metro_insets.png
#           <fig dir>/bivariate_tract_national_with_metros.png
#           <season dir>/tract_bivariate_summary.csv
#---------------------------------------------------------------------------

#---- 1. Setup #----
rm(list = ls())

pkgs <- c("dplyr", "readr", "tidyr", "ggplot2", "sf", "patchwork", "ragg")
invisible(lapply(pkgs, library, character.only = TRUE))

# The default macOS png device drops UTF-8 glyphs; ragg::agg_png renders them.
png_dev <- ragg::agg_png

# Built via code points rather than typed literally: Rscript launched from a
# shell with LANG unset parses the source as US-ASCII and would mangle literal
# UTF-8 bytes. This keeps the source pure ASCII and is locale-independent.
arrow_sym <- intToUtf8(0x2192)              # rightwards arrow
times_sym <- intToUtf8(0x00D7)              # multiplication sign
dash_sym  <- intToUtf8(0x2014)              # em dash

data_dir <- "/Users/hangli/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Aune, Kyle T's files - 260204 - US Heat Clusters/Data"
fig_dir  <- "/Users/hangli/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Aune, Kyle T's files - 260204 - US Heat Clusters/Output/Figures"

season_dir <- file.path(data_dir, "processed", "season cluster")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Contiguous US, matching the county figure.
crs_conus <- 5070


#---- 2. Data #----
# The vulnerability index is read from cache only. It is built (and its variable
# signature validated) in section 8 of domain_exposure_heatmap.R; duplicating
# that logic here would risk the two figures silently diverging.
vuln_cache <- file.path(season_dir, "tract_vulnerability_index.csv")
if (!file.exists(vuln_cache)) {
  stop("missing ", vuln_cache,
       "\nRun section 8 of domain_exposure_heatmap.R first to build it.")
}

vuln_tract <- read_csv(vuln_cache,
                       col_types = cols(GEOID = col_character(),
                                        .default = col_double()))

county <- read_csv(
  file.path(season_dir, "county_clustered.csv"),
  col_types = cols(StCoFIPS = col_character(),
                   cluster  = col_integer(),
                   .default = col_double())
)

# Hazard composite, identical to domain_exposure_heatmap.R: each component
# z-scored ACROSS COUNTIES, then summed. Standardising at county level (not
# tract level) is what keeps this figure's hazard axis on the same footing as
# the county figure - tracts are not the population the z-scores refer to.
county <- county %>%
  mutate(
    hazard = as.numeric(scale(frequency_pct)) +
             as.numeric(scale(mean_Tmax_intensity)) +
             as.numeric(scale(longest_duration))
  )

message("vulnerability index: ", format(nrow(vuln_tract), big.mark = ","), " tracts")
message("county hazard:       ", format(nrow(county),     big.mark = ","), " counties")


#---- 3. Assign county hazard down to tracts #----
tract_dat <- vuln_tract %>%
  mutate(StCoFIPS = substr(GEOID, 1, 5)) %>%
  inner_join(
    county %>% dplyr::select(StCoFIPS, hazard, frequency_pct,
                             mean_Tmax_intensity, longest_duration, cluster),
    by = "StCoFIPS"
  ) %>%
  dplyr::filter(is.finite(hazard), is.finite(hvi))

message(sprintf("tracts with both hazard and vulnerability: %s across %s counties",
                format(nrow(tract_dat), big.mark = ","),
                format(n_distinct(tract_dat$StCoFIPS), big.mark = ",")))
message(sprintf("  dropped %s tracts whose county has no hazard record",
                format(nrow(vuln_tract) - nrow(tract_dat), big.mark = ",")))


#---- 4. Why this figure exists: variance decomposition of the index #----
# One-way ANOVA of hvi on county. The between-county share is the most that
# the population-weighted county index in domain_exposure_heatmap.R can carry;
# everything else is within-county contrast that only a tract map can show.
var_decomp <- function(x, g) {
  ok <- is.finite(x); x <- x[ok]; g <- g[ok]
  gm <- mean(x)
  ss_between <- sum(tapply(x, g, function(v) length(v) * (mean(v) - gm)^2))
  ss_total   <- sum((x - gm)^2)
  c(between = ss_between / ss_total, within = 1 - ss_between / ss_total)
}

vd_hvi <- var_decomp(tract_dat$hvi, tract_dat$StCoFIPS)
cat("\n=== variance decomposition of the tract vulnerability index ===\n")
cat(sprintf("  between counties: %5.1f%%\n", 100 * vd_hvi["between"]))
cat(sprintf("  within  counties: %5.1f%%   <- invisible on the county map\n",
            100 * vd_hvi["within"]))

# Same decomposition per component, so the caption can name which domains are
# the within-county ones rather than asserting it.
cat("\n  by component (within-county share):\n")
vcols <- grep("^V[0-9]+$", names(tract_dat), value = TRUE)
for (v in vcols) {
  cat(sprintf("    %-4s %5.1f%%\n", v,
              100 * var_decomp(tract_dat[[v]], tract_dat$StCoFIPS)["within"]))
}


#---- 5. Terciles #----
# Each variable is cut on the distribution of the units it was actually measured
# on: hazard on the 3,106 COUNTIES, vulnerability on the tracts.
#
# Cutting hazard on the tract distribution instead is tempting (a break would
# then mean "a third of tracts") but it silently reweights the hazard axis by
# tract count. High-hazard counties are largely rural and tract-poor while the
# low-hazard Northeast and upper Midwest are tract-rich, so the tract-weighted
# breaks sit far below the county ones (printed below) and the top hazard class
# swells. That would also break comparability with the county figure, which is
# the whole point of this map. Set haz_break_unit <- "tract" to see that version.
haz_break_unit <- "county"

terc_breaks <- function(x) quantile(x, c(0, 1/3, 2/3, 1), na.rm = TRUE)
cut_terc <- function(x, br) cut(x, breaks = br, labels = 1:3,
                                include.lowest = TRUE)

br_haz_county <- terc_breaks(county$hazard)
br_haz_tract  <- terc_breaks(tract_dat$hazard)
br_haz <- if (haz_break_unit == "county") br_haz_county else br_haz_tract

# Counties outside the tract-level hazard range would fall outside the county
# breaks; they cannot here (tract hazard is a subset of county hazard), but
# clamp the outer edges so the cut is total regardless.
br_haz[1] <- min(br_haz[1], min(tract_dat$hazard))
br_haz[4] <- max(br_haz[4], max(tract_dat$hazard))

tract_dat <- tract_dat %>%
  mutate(
    haz_t  = cut_terc(hazard, br_haz),
    vuln_t = cut_terc(hvi, terc_breaks(hvi)),
    biv    = paste(haz_t, vuln_t, sep = "-")
  )

stopifnot(!any(is.na(tract_dat$haz_t)), !any(is.na(tract_dat$vuln_t)))

cat("\n=== hazard tercile breaks ===\n")
cat(sprintf("  across counties: %+.3f, %+.3f%s\n",
            br_haz_county[2], br_haz_county[3],
            if (haz_break_unit == "county") "   <- used here" else ""))
cat(sprintf("  across tracts:   %+.3f, %+.3f%s\n",
            br_haz_tract[2], br_haz_tract[3],
            if (haz_break_unit == "tract") "   <- used here" else ""))
cat("\n  share of tracts per hazard tercile under each rule:\n")
cat(sprintf("    county breaks: %s\n",
            paste(sprintf("%.0f%%", 100 * prop.table(table(
              cut_terc(tract_dat$hazard, br_haz_county)))), collapse = " / ")))
cat(sprintf("    tract breaks:  %s\n",
            paste(sprintf("%.0f%%", 100 * prop.table(table(
              cut_terc(tract_dat$hazard, br_haz_tract)))), collapse = " / ")))

# Stevens-style bivariate palette, keyed "hazard-vulnerability". Identical to
# domain_exposure_heatmap.R section 10 so the two maps read as one series.
biv_pal <- c(
  "1-1" = "#e8e8e8", "2-1" = "#e4acac", "3-1" = "#c85a5a",
  "1-2" = "#b0d5df", "2-2" = "#ad9ea5", "3-2" = "#985356",
  "1-3" = "#64acbe", "2-3" = "#627f8c", "3-3" = "#574249"
)

# The county figure uses grey92 for missing counties, but 16% of tracts here
# have no vulnerability score and grey92 (#EBEBEB) is visually indistinguishable
# from the 1-1 class (#E8E8E8) - at tract scale that would read as a large
# low-hazard/low-vulnerability area rather than as absent data. Step it down far
# enough to be unambiguous while staying obviously outside the palette.
na_col <- "#c9c9c9"


#---- 6. Population per bivariate class #----
tot_pop <- sum(tract_dat$pop, na.rm = TRUE)

biv_summary <- tract_dat %>%
  group_by(haz_t, vuln_t, biv) %>%
  summarise(n_tract = n(),
            pop     = sum(pop, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(pop_mil = pop / 1e6,
         pop_pct = 100 * pop / tot_pop) %>%
  arrange(desc(pop))

cat("\n=== population by hazard ", times_sym, " vulnerability tercile (tracts) ===\n", sep = "")
print(as.data.frame(biv_summary %>%
                      dplyr::select(biv, haz_t, vuln_t, n_tract, pop_mil, pop_pct)),
      row.names = FALSE, digits = 4)

# Headline: high hazard AND high vulnerability, and the two mismatch corners.
hh <- biv_summary$pop_mil[biv_summary$biv == "3-3"]
lh <- biv_summary$pop_mil[biv_summary$biv == "1-3"]   # low hazard, high vuln
hl <- biv_summary$pop_mil[biv_summary$biv == "3-1"]   # high hazard, low vuln
cat(sprintf("\n  both high (3-3): %.1f M | latent risk, low hazard + high vuln (1-3): %.1f M | high hazard + low vuln (3-1): %.1f M\n",
            hh, lh, hl))

write_csv(biv_summary, file.path(season_dir, "tract_bivariate_summary.csv"))


#---- 7. Geometry #----
# Tract boundaries come from tract_master_2020.rds, which already carries an
# sf geometry column for CONUS only. Do NOT pull tracts from tigris - it is a
# ~50-state download for data already on disk.
message("reading tract geometry from tract_master_2020.rds ...")
tm <- readRDS(file.path(data_dir, "tract_master_2020.rds"))

# !valid_ct is dropped rather than drawn grey. Every land_km2 == 0 tract carries
# that flag, so these are the water polygons - Lake Michigan inside Cook County,
# New York harbour, the Channel Islands - plus other unpopulated tracts. Drawing
# them as "no data" would put a grey lake in the middle of the Chicago panel and
# imply a measurement gap where there is simply no population. This is the same
# filter the vulnerability index itself was built under.
tract_sf <- tm %>%
  dplyr::filter(valid_ct) %>%
  dplyr::select(GEOID, STATEFP) %>%
  st_transform(crs_conus)

n_invalid <- sum(!tm$valid_ct)
rm(tm); invisible(gc())

biv_sf <- tract_sf %>%
  left_join(tract_dat %>% dplyr::select(GEOID, StCoFIPS, hazard, hvi,
                                        haz_t, vuln_t, biv, pop),
            by = "GEOID")

# Account for every grey tract, so the caption states the reason rather than
# guessing at it. These two are mutually exclusive.
unfilled <- biv_sf %>%
  st_drop_geometry() %>%
  dplyr::filter(is.na(biv)) %>%
  mutate(reason = if_else(
    !substr(GEOID, 1, 5) %in% county$StCoFIPS,
    "county has no hazard record",
    "incomplete vulnerability block"
  ))

message(sprintf("tracts drawn: %s | classified: %s | grey: %s | excluded as unpopulated/water: %s",
                format(nrow(biv_sf),            big.mark = ","),
                format(sum(!is.na(biv_sf$biv)), big.mark = ","),
                format(nrow(unfilled),          big.mark = ","),
                format(n_invalid,               big.mark = ",")))
cat("\n=== unmapped (grey) tracts by reason ===\n")
print(as.data.frame(unfilled %>% count(reason, name = "n_tract") %>%
                      arrange(desc(n_tract))), row.names = FALSE)

# State outlines for the national panel. County outlines are deliberately left
# off: hazard is constant within a county, so drawing 3,100 county borders on
# top would just re-trace the hazard field and swamp the tract texture.
options(tigris_use_cache = TRUE)
suppressMessages(library(tigris))
states_sf <- suppressMessages(
  states(cb = TRUE, year = 2020, class = "sf", progress_bar = FALSE)
) %>%
  dplyr::filter(!STATEFP %in% c("02", "15", "60", "66", "69", "72", "78")) %>%
  st_transform(crs_conus)


#---- 8. WARNING - no inferential statistics at this scale #----
#===========================================================================
# Do NOT run a LISA, a permutation test, or any other inference on this table.
# Hazard is CONSTANT within each county: 69k tract rows carry only ~2.9k
# independent hazard values. A conditional-permutation null would be built by
# reshuffling replicated values, so its reference distribution is far too
# narrow and pseudo-p values collapse toward zero - the map would come back
# "significant mismatch" almost everywhere, as an artefact of the replication
# rather than a finding.
#
# All inferential claims about spatial mismatch stay with the county-level
# bivariate LISA in section 11 of domain_exposure_heatmap.R, where hazard and
# vulnerability are measured on the same units. This script is descriptive.
#===========================================================================


#---- 9. Legend #----
leg_df <- expand_grid(haz_t = 1:3, vuln_t = 1:3) %>%
  mutate(biv = paste(haz_t, vuln_t, sep = "-"))

make_legend <- function(base_size = 8) {
  ggplot(leg_df, aes(x = haz_t, y = vuln_t, fill = biv)) +
    geom_tile() +
    scale_fill_manual(values = biv_pal, guide = "none") +
    labs(x = paste0("Hazard ", arrow_sym),
         y = paste0("Vulnerability ", arrow_sym)) +
    coord_fixed() +
    theme_minimal(base_size = base_size) +
    theme(axis.text        = element_blank(),
          axis.title       = element_text(size = base_size - 0.5),
          panel.grid       = element_blank(),
          plot.background  = element_rect(fill = "white", color = NA))
}

p_leg <- make_legend()


#---- 10. National panel #----
sub_nat <- paste0(
  "Terciles of county heatwave hazard (frequency + intensity + duration), assigned to every tract in the county,\n",
  "crossed with terciles of the tract heat vulnerability index. ",
  format(nrow(tract_dat), big.mark = ","), " tracts in ",
  format(n_distinct(tract_dat$StCoFIPS), big.mark = ","), " counties.\n",
  "Hazard is cut on the county distribution (same breaks as the county figure), vulnerability on the tract distribution."
)

cap_nat <- paste0(
  "Vulnerability varies mostly WITHIN counties (", sprintf("%.0f%%", 100 * vd_hvi["within"]),
  " of index variance), which the population-weighted county map cannot show; hazard is available only at county\n",
  "level, so it is uniform inside each county here. Descriptive figure ", dash_sym,
  " no spatial inference is run at tract scale, because replicated hazard values would invalidate the permutation null.\n",
  "Grey = tract lacking a vulnerability score, ", sprintf("%s tracts", format(nrow(unfilled), big.mark = ",")),
  ", overwhelmingly because air-conditioning coverage is unavailable there. Water and unpopulated tracts are not drawn."
)

p_nat <- ggplot() +
  geom_sf(data = biv_sf, aes(fill = biv), color = NA) +
  geom_sf(data = states_sf, fill = NA, color = "white", linewidth = 0.18) +
  scale_fill_manual(values = biv_pal, na.value = na_col, guide = "none") +
  labs(
    title    = "Heat hazard and heat vulnerability across CONUS census tracts",
    subtitle = sub_nat,
    caption  = cap_nat
  ) +
  theme_void(base_size = 11) +
  # theme_void() leaves the canvas transparent, which renders as black in most
  # viewers and hides the black title/caption text. Force an opaque white canvas.
  theme(plot.background = element_rect(fill = "white", color = NA),
        plot.title      = element_text(face = "bold", size = 13),
        plot.subtitle   = element_text(size = 9,   color = "grey30"),
        plot.caption    = element_text(size = 7.4, color = "grey40", hjust = 0,
                                       margin = margin(t = 8)))

p_nat_full <- p_nat + inset_element(p_leg, left = 0.02, bottom = 0.04,
                                    right = 0.19, top = 0.32)

if (interactive()) print(p_nat_full)

ggsave(file.path(fig_dir, "bivariate_tract_hazard_vulnerability_map.png"),
       p_nat_full, width = 11, height = 7.0, dpi = 600, device = png_dev)
message("wrote national panel")


#---- 11. Metro insets #----
# Reid et al. (2009) mapped their vulnerability factors for Chicago, New York,
# Los Angeles and Dallas/Fort Worth; using the same four keeps this figure
# directly comparable to that paper. Metros are defined by county FIPS rather
# than by bounding box so every tract shown is one that carries a hazard value.
metros <- list(
  "Chicago" = c("17031", "17043", "17089", "17093", "17097", "17111", "17197",
                "18089", "18127"),
  "New York City" = c("36005", "36047", "36061", "36081", "36085", "36059",
                      "36119", "34003", "34013", "34017", "34031", "34039"),
  "Los Angeles" = c("06037", "06059"),
  "Dallas-Fort Worth" = c("48113", "48439", "48085", "48121", "48397", "48257",
                          "48139", "48251", "48367")
)

counties_sf <- suppressMessages(
  counties(cb = TRUE, year = 2020, class = "sf", progress_bar = FALSE)
) %>%
  st_transform(crs_conus)

frame_drop <- c("06037599000", "06037599100")   # see metro_panel() below

# Number of distinct hazard terciles present tells the reader whether the panel
# has any hazard contrast at all - several metros sit entirely in one tercile,
# which is exactly the point: inside them the map is reading vulnerability only.
metro_panel <- function(nm, fips) {
  # StCoFIPS is NA on unmapped tracts (it arrives via the join), so select the
  # panel's tracts on GEOID instead - otherwise the grey no-data tracts would
  # be dropped and the panels would show holes rather than missing data.
  sub_sf <- biv_sf %>% dplyr::filter(substr(GEOID, 1, 5) %in% fips)
  cty    <- counties_sf %>% dplyr::filter(GEOID %in% fips)

  d      <- sub_sf %>% st_drop_geometry() %>% dplyr::filter(!is.na(biv))
  n_grey <- nrow(sub_sf) - nrow(d)
  pop_m  <- sum(d$pop, na.rm = TRUE) / 1e6
  haz_lab <- if (n_distinct(d$haz_t) == 1) {
    paste0("hazard tercile ", as.character(d$haz_t[1]), " throughout")
  } else {
    paste0("hazard terciles ",
           paste(sort(unique(as.integer(d$haz_t))), collapse = "/"))
  }

  # Frame on the TRACTS, not on the county outlines. The cartographic-boundary
  # county polygons run out to the legal water limit, so Los Angeles County
  # alone would drag the extent far offshore and shrink the metro to a smudge.
  #
  # frame_drop is excluded from the viewport only - these tracts are still
  # counted, still summarised, and still drawn if they fall inside the frame.
  # It holds the two Channel Islands tracts (Avalon/Santa Catalina and San
  # Clemente), which sit ~70 km offshore and between them would cost a third of
  # the Los Angeles panel's height to show 3,875 of its 12.9 M residents.
  bb <- st_bbox(sub_sf %>% dplyr::filter(!GEOID %in% frame_drop))
  pad <- 0.02 * max(bb["xmax"] - bb["xmin"], bb["ymax"] - bb["ymin"])

  ggplot() +
    geom_sf(data = sub_sf, aes(fill = biv), color = NA) +
    geom_sf(data = cty, fill = NA, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = biv_pal, na.value = na_col, guide = "none") +
    coord_sf(xlim = c(bb["xmin"] - pad, bb["xmax"] + pad),
             ylim = c(bb["ymin"] - pad, bb["ymax"] + pad),
             expand = FALSE) +
    labs(title = nm,
         # Split across two short lines. Panels are only as wide as their own
         # map extent, so a single long line overruns the cell and collides
         # with the neighbouring panel's heading.
         subtitle = sprintf("%s tracts, %.1f M people\n%s; %s no data",
                            format(nrow(d), big.mark = ","), pop_m, haz_lab,
                            format(n_grey, big.mark = ","))) +
    theme_void(base_size = 10) +
    theme(plot.background = element_rect(fill = "white", color = NA),
          # Align titles to the whole cell rather than to the map's bounding
          # box, which differs per metro and leaves the headings ragged.
          plot.title.position = "plot",
          plot.title      = element_text(face = "bold", size = 11,
                                         margin = margin(b = 1)),
          plot.subtitle   = element_text(size = 7.4, color = "grey35",
                                         lineheight = 1.15,
                                         margin = margin(b = 4)),
          plot.margin     = margin(4, 6, 4, 6))
}

metro_plots <- Map(metro_panel, names(metros), metros)

# The legend gets its own short row rather than an inset_element. Adding an
# inset to an existing patchwork attaches it to the LAST panel, not the
# composition, so it lands inside the Dallas cell and is scaled to that cell -
# a real trap, and it silently produces a legible-looking but unreadable key.
legend_row <- make_legend(8) + plot_spacer() + plot_layout(widths = c(1, 3.2))

p_metro <- (metro_plots[[1]] | metro_plots[[2]]) /
  (metro_plots[[3]] | metro_plots[[4]]) /
  legend_row +
  plot_layout(heights = c(1, 1, 0.30)) +
  plot_annotation(
    title    = paste0("Within-county vulnerability texture ", dash_sym,
                      " four metropolitan areas"),
    subtitle = paste0(
      "Same bivariate classes as the national map. Hazard is uniform within each county (white outlines), so all\n",
      "variation inside a county is vulnerability. Metros follow Reid et al. (2009) for comparability.\n",
      "Grey tracts lack a vulnerability score, almost entirely because air-conditioning coverage is unavailable there."),
    theme = theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.title      = element_text(face = "bold", size = 13),
      plot.subtitle   = element_text(size = 9, color = "grey30",
                                     margin = margin(b = 6)))
  )

if (interactive()) print(p_metro)

ggsave(file.path(fig_dir, "bivariate_tract_metro_insets.png"),
       p_metro, width = 10.5, height = 8.6, dpi = 600, device = png_dev)
message("wrote metro insets")


#---- 12. Combined national + metro figure #----
# wrap_elements() renders the metro block to a single grob first. Nesting a
# patchwork directly discards its plot_annotation, so without this the metro
# half of the combined figure loses its heading and the note about grey tracts.
p_combined <- (p_nat_full / wrap_elements(p_metro)) +
  plot_layout(heights = c(1, 1.25)) &
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(fig_dir, "bivariate_tract_national_with_metros.png"),
       p_combined, width = 11, height = 15.5, dpi = 450, device = png_dev)
message("wrote combined figure")

message("done - figures written to: ", fig_dir)

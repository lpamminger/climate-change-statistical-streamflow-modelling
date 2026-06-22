# Analysis the wheat belt in WA to see if it is okay

# Load packages ----------------------------------------------------------------
pacman::p_load(tidyverse, ozmaps, sf, ggmagnify, trend, patchwork)


# Import functions -------------------------------------------------------------
source("./Functions/utility.R")


# Load data --------------------------------------------------------------------
gauge_information <- readr::read_csv(
  "./Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE
)

WA_gauges <- gauge_information |> 
  filter(state == "WA") |> 
  pull(gauge)

CAMELSAUS_boundary <- st_read(
  dsn = "./Data/Maps/02_location_boundary_area/shp/",
  as_tibble = TRUE
) |> 
  rename(gauge = CatchID) |> 
  filter(gauge %in% WA_gauges)


evidence_ratio <- read_csv(
  "./Modelling/Results/CMAES/evidence_ratio_results.csv",
  show_col_types = FALSE
) |> 
  # factor impact of CO2 term 
  mutate(
    impact_of_CO2_term = factor(
      impact_of_CO2_term,
      levels = c("No CO2 Term", "Negative-Intercept", "Positive-Intercept", "Negative-Slope", "Positive-Slope")
    )
  ) |>
  # factor qualitative evidence ratio values
  mutate(
    binned_evidence_ratio = factor(
      binned_evidence_ratio,
      levels = c("Weak", "Moderate", "Moderately Strong", "Strong", "Very Strong", "Extremely Strong")
    )
  )



WA_wheat_belt <- st_read(
  dsn = "Data/Maps/Wheatbelt_of_WA"
)

aus_map <- generate_aus_map_sf()

aus_map_crs <- st_crs(aus_map)

working_crs <- aus_map_crs[[1]]



# How many catchments are located in the WA wheat belt? ------------------------
WA_wheat_belt <- WA_wheat_belt |> 
  st_transform(crs = working_crs) 

CAMELSAUS_boundary <- CAMELSAUS_boundary |> 
  st_transform(crs = working_crs)

## Intersection
sf_use_s2(FALSE)
interect_wheat_belt <- st_intersection(
  CAMELSAUS_boundary,
  WA_wheat_belt
) 

gauges_that_interect_wheat_belt <- interect_wheat_belt |> 
  pull(gauge)

# 27/72 of the WA catchments intersect the wheat belt
# only one has a strong evidence ratio --> 701002
evidence_ratio |> 
  filter(gauge %in% gauges_that_interect_wheat_belt) |> 
  arrange(desc(evidence_ratio))

effected_catchments <- CAMELSAUS_boundary |> 
  filter(gauge %in% gauges_that_interect_wheat_belt) |> 
  add_column(
    type = "In Wheatbelt"
  )

uneffected_catchments <- CAMELSAUS_boundary |> 
  filter(!gauge %in% gauges_that_interect_wheat_belt) |> 
  add_column(
    type = "Not in Wheatbelt"
  )


  
## Compare evidence ratio of effected and unaffected catchments 
effected_gauges <- effected_catchments |> pull(gauge)
uneffected_gauges <- uneffected_catchments |> pull(gauge)

WA_evidence_ratio <- evidence_ratio |> 
  filter(gauge %in% c(effected_gauges, uneffected_gauges))

wheatbelt_evi_comparison <- WA_evidence_ratio |> 
  mutate(
    in_wheatbelt = case_when(
      gauge %in% effected_gauges ~ TRUE,
      gauge %in% uneffected_gauges ~ FALSE,
      .default = NA
    )
  ) |> 
  filter(state == "WA")

### Nice breaks for boxplot
evi_ratio_range <- evidence_ratio |> pull(evidence_ratio) |> range()
y_axis_evi_min <- round_any(evi_ratio_range[1], accuracy = 10, f = floor)
y_axis_evi_max <- round_any(evi_ratio_range[2], accuracy = 1E16, f = ceiling)
nice_breaks_log10 <- seq(from = log10(abs(y_axis_evi_min)), to = log10(y_axis_evi_max), by = 3)
nice_breaks_log10[1] <- -1
nice_breaks <- 10^nice_breaks_log10
y_axis_scale_transform <- scale_y_continuous(
  transform = scales::pseudo_log_trans(base = 10),
  breaks = nice_breaks
)


### Add label
single_label <- function(x_pos, y_pos, label_name) { # for adding a, b, c labels
  tribble(
    ~x_pos, ~y_pos, ~label_name,
    x_pos,  y_pos,  label_name
  )
}


wheatbelt_evi_boxplot <- wheatbelt_evi_comparison |> 
  ggplot(aes(x = in_wheatbelt, y = evidence_ratio)) +
  geom_boxplot(
    staplewidth = 0.5,
    outlier.alpha = 0.7,
    outlier.colour = "black",
    outlier.fill = "grey",
    outlier.shape = 21,
    whisker.linewidth = 0.2,
    box.linewidth = 0.2,
    median.linewidth = 0.4,
    staple.linewidth = 0.2,
    fill = "grey90"
  ) +
  geom_text(
    data = single_label(x_pos = FALSE, y_pos = 1E16, label_name = "a"),
    mapping = aes(x = x_pos, y = y_pos, label = label_name),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 12,
    size.unit = "pt",
    hjust = 7,
    vjust = -0.7
  ) +
  y_axis_scale_transform +
  theme_bw() +
  labs(
    x = "In Wheatbelt",
    y = "Evidence Ratio"
  ) +
  theme(
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8)
  )


# Make make of wheat belt ------------------------------------------------------
# make WA_wheat_belt have same layout as effected catchments
WA_altered_wheat_belt <- WA_wheat_belt |> 
  select(geometry) |> 
  add_column(
    gauge = "gauge",
    type = "Wheatbelt"
  ) |> 
  relocate(
    gauge, geometry, type
  )

# add geom_sf single then alter fill
# for illustrative purposes remove nested catchments
# the evidence ratio dots shows if there is a nested catchment
grouped_catchments <- effected_catchments |> 
  rbind(uneffected_catchments) |> 
  rbind(WA_altered_wheat_belt)

remove_nested <- gauge_information |> 
  filter(state == "WA") |> 
  filter(status == "sub_catchment") |> 
  pull(gauge)

grouped_catchments <- grouped_catchments |> 
  filter(!gauge %in% remove_nested)

custom_palette <- function(x) {
  rev(c("yellow", "blue", "red", "#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#f7f7f7"))
}


wheat_belt_map <- aus_map |> 
  filter(state == "WA") |> 
  ggplot() +
  geom_sf() +
  geom_sf(
    data = grouped_catchments,
    aes(fill = type),
    colour = "black",
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  geom_point(
    aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term),
    data = WA_evidence_ratio,
    colour = "black",
    size = 2,
    inherit.aes = FALSE,
    stroke = 0.1
  ) +
  geom_text(
    data = single_label(x_pos = 115, y_pos = -27.5, label_name = "b"),
    mapping = aes(x = x_pos, y = y_pos, label = label_name),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 12,
    size.unit = "pt",
    hjust = 4
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    # large space is to centre the title
    fill = "                    Evidence Ratio                        Wheatbelt Catchment Status",
    shape = bquote("Impact of" ~ CO[2] ~ "Term")
  ) +
  # zoom it down
  coord_sf(xlim = c(114, 120), ylim = c(-35, -27.5)) +
  theme_bw() +
  scale_fill_manual(
    values = custom_palette(),
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c(21, 22, 23, 25, 24)#,
    #drop = FALSE
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 5, shape = 21), nrow = 3), # Wrap legend with nrow
    shape = guide_legend(override.aes = list(size = 5, fill = "grey50"), nrow = 2)
  ) +
  theme(
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    legend.title = element_text(hjust = 0.5)
  )


wheatbelt_analysis <- wheatbelt_evi_boxplot + wheat_belt_map +
  theme(
    legend.title.position = "top",
    legend.position = "bottom",
    legend.text = element_text(size = 8, hjust = 0.5),
    legend.title = element_text(size = 9),
    legend.box = "horizontal",
    legend.justification = c(1,1)
    )

ggsave(
  filename = "./Figures/Other/wheat_belt_map.pdf",
  plot = wheatbelt_analysis,
  device = "pdf",
  width = 180,
  height = 170,
  units = "mm"
)


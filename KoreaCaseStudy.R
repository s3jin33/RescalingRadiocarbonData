# This R script generates Figures 9 and 10.
# To produce Figure 9, load the file "Yeongsanriver.csv".
# To produce Figure 10, load the file "Geumriver.csv".
# For spatio-temporal KDE, make sure to use a different spatial window for each dataset,
# as the appropriate window size differs between the Yeongsan and Geum River Basins.

# Load required libraries
library(dplyr)
library(rcarbon)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(readr)
library(patchwork)
library(zoo)
library(ggrepel)
library(ggsci)
library(ggrepel)
library(RColorBrewer)
library(purrr)

# set directory
setwd("/Users/daljaepark/Library/CloudStorage/OneDrive-개인/논문투고/2025 Bootstrap resampling 리부트/마한백제데이터")


#########################################################################
##################  Figure 9 and 10 (a): SPD comparison 
#########################################################################

#------------------------------------------------------------------------
#------------  step 1. prepare datasets
#------------------------------------------------------------------------

# Load published radiocarbon dates.
# Load one of the datasets: either "Yeongsanriver.csv" or "Geumriver.csv",
# depending on the target river basin.
sampled_data <- read.csv("C14_data_Yeongsanriver.csv", , fileEncoding="cp949") 
sampled_data <- read.csv("C14_data_Geumriver.csv", , fileEncoding="cp949") 


house_counts <- sampled_data %>%
  distinct(SiteCode, total_houses)
print(house_counts)


rad_count <- sampled_data %>%
  distinct(SiteCode, count)
print(rad_count)



# calculate weights for each settlement
weights <- sampled_data %>%
  group_by(SiteCode) %>%
  summarise(total_houses = unique(total_houses), count = n()) %>%
  mutate(weight = total_houses / count)

# prepare bootstrap dataset
# define function for resampling
resample_data <- function(sampled_data, n_resamples, seed_rs) {
  set.seed(seed_rs) 
  resampled_data <- data.frame()
  
  for (resample_iter in 1:n_resamples) {
    print(paste("Starting resampling iteration:", resample_iter))
    
    # resampling by site(SiteCode)
    resampled_iteration <- sampled_data %>%
      group_by(SiteCode) %>%
      group_modify(~ {
        n_houses <- unique(.x$`total_houses`)  
        
        if (length(n_houses) == 0 || is.na(n_houses) || n_houses <= 0) {
          return(NULL)
        }
        
        resampled_subset <- .x[sample(nrow(.x), n_houses, replace = TRUE), ]
        
        resampled_subset <- resampled_subset %>%
          mutate(
            Resample_Iteration = resample_iter,
            Sampling_Num = seq_len(n())
          )
        
        return(resampled_subset)
      }) %>%
      ungroup()
   
    resampled_data <- bind_rows(resampled_data, resampled_iteration)
    print(paste("Completed resampling iteration:", resample_iter))
  }
  return(resampled_data)
}


# resample 30 times, set seed as 1234
resampled_output <- resample_data(sampled_data, n_resamples = 30, seed_rs = 1234)



#------------------------------------------------------------------------
#------------  step 2. calibrate dates and calculate SPD 
#------------------------------------------------------------------------


# Function to calibrate dates with errors
calibrate_dates <- function(settlement_data) {
  if (nrow(settlement_data) > 0) {
    calibrated <- calibrate(settlement_data$BP, errors = settlement_data$Error, timeRange = c(2600, 1200), normalised = TRUE)
    return(calibrated)
  }
  return(NULL)
}

# Function to extract SPD info while retaining `SiteCode`
extract_spd_info_with_site <- function(cal_date, site) {
  if (inherits(cal_date, "CalDates")) {
    spd_result <- spd(cal_date, timeRange = c(2600, 1200), spdnormalised = FALSE)
    spd_data <- data.frame(
      calBP = spd_result$grid$calBP,
      PrDens = spd_result$grid$PrDens,
      SiteCode = site  
    )
    return(spd_data)
  }
  return(NULL)
}

# Compute SPD for Sampled Data 
sampled_spd_results <- sampled_data %>%
  group_by(SiteCode) %>%
  summarise(Calibrated_Dates = list(calibrate_dates(pick(everything()))), .groups = 'drop') %>%
  mutate(SPD = map2(Calibrated_Dates, SiteCode, extract_spd_info_with_site))

# Compute SPD for Resampled Data 
resampled_spd_results <- resampled_output %>%
  group_by(SiteCode) %>%
  summarise(Calibrated_Dates = list(calibrate_dates(pick(everything()))), .groups = 'drop') %>%
  mutate(SPD = map2(Calibrated_Dates, SiteCode, extract_spd_info_with_site))

# Compute Weighted SPD using Sampled Data 
weighted_spd_results <- sampled_spd_results %>%
  left_join(weights, by = "SiteCode") %>%
  mutate(SPD = pmap(list(SPD, weight, SiteCode), function(spd, w, site) {
    if (!is.null(spd)) {
      spd$PrDens <- spd$PrDens * w  
      spd$SiteCode <- site
    }
    return(spd)
  }))


#drawing SPD. by settlement with 25-year rolling mean
combine_spd_per_settlement <- function(spd_results) {
  bind_rows(spd_results$SPD) %>%
    group_by(calBP, SiteCode) %>%
    summarise(PrDens = sum(PrDens, na.rm = TRUE), .groups = 'drop') %>%
    arrange(SiteCode, desc(calBP)) %>%
    group_by(SiteCode) %>%
    mutate(PrDens = zoo::rollmean(PrDens, k = 25, fill = "extend", align = "center")) %>%
    ungroup()
}



# ✅ Aggregate SPDs for Sampled, Resampled, and Weighted
sampled_spd_total_by_settlement <- combine_spd_per_settlement(sampled_spd_results)
resampled_spd_total_by_settlement <- combine_spd_per_settlement(resampled_spd_results)
weighted_spd_total_by_settlement <- combine_spd_per_settlement(weighted_spd_results)


# ✅ Function to normalize SPDs 
normalize_spd <- function(spd_df) {
  total_density <- sum(spd_df$PrDens, na.rm = TRUE)
  spd_df <- spd_df %>%
    mutate(PrDens = PrDens / total_density)  # ✅ Normalize SPD
  return(spd_df)
}

# ✅ Normalize the SPDs
sampled_spd_normalized <- normalize_spd(sampled_spd_total_by_settlement %>%
                                          group_by(calBP) %>%
                                          summarise(PrDens = sum(PrDens), .groups = 'drop') %>%
                                          mutate(Dataset = "Unrescaled original"))

resampled_spd_normalized <- normalize_spd(resampled_spd_total_by_settlement %>%
                                            group_by(calBP) %>%
                                            summarise(PrDens = sum(PrDens), .groups = 'drop') %>%
                                            mutate(Dataset = "Bootstrapped"))

weighted_spd_normalized <- normalize_spd(weighted_spd_total_by_settlement %>%
                                           group_by(calBP) %>%
                                           summarise(PrDens = sum(PrDens), .groups = 'drop') %>%
                                           mutate(Dataset = "Weighted"))



#------------------------------------------------------------------------
#------------  step 3. draw total SPD to compare three datasets 
#------------------------------------------------------------------------

total_spd_normalized <- bind_rows(sampled_spd_normalized, resampled_spd_normalized, weighted_spd_normalized)
spd_colors <- c("Unrescaled original" = "#ffbb6f", "Bootstrapped" = "#999999", "Weighted" = "#5e4c5f")
total_spd_normalized$Dataset <- factor(total_spd_normalized$Dataset, levels = c("Unrescaled original", "Bootstrapped", "Weighted"))

plot_total_spd_normalized <- ggplot(total_spd_normalized, aes(x = calBP, y = PrDens, color = Dataset, linetype = Dataset)) +
  geom_line(linewidth = 1) +  
  scale_color_manual(values = spd_colors) +
  scale_linetype_manual(values = c("Unrescaled original" = "solid", "Bootstrapped" = "dashed", "Weighted" = "dotdash")) +
  scale_x_reverse(limits = c(2200, 1400)) +  
  labs(title = "Total SPD Comparison",
       x = "cal BP", y = "Normalized SPD") +
  theme_minimal(base_family = "noto") +
  theme(
    text = element_text(family = "noto"),
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA),  
    legend.position = "top",  # ✅ Move legend to the top
    legend.text = element_text(size = 12, face = "bold"),  
    legend.title = element_blank(), 
    legend.key.width = unit(1.5, "cm"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),  
    axis.text = element_text(size = 11, face = "bold"), 
    axis.title = element_text(size = 11, face = "bold")  
  ) 

print(plot_total_spd_normalized)
ggsave("Figure9a.png", plot_total_spd_normalized, width = 12, height = 4, dpi = 300)



#########################################################################
##################  Figure 9 and 10 (b): settlement-SPD comparison 
#########################################################################

#------------------------------------------------------------------------
#------------  step 1. normalize settlement-level SPD by dataset
#------------------------------------------------------------------------

normalize_preserve_shape <- function(all_spd_df) {
  # Compute total SPD contribution per settlement
  settlement_totals <- all_spd_df %>%
    group_by(Dataset, SiteCode) %>%
    summarise(Total_PrDens = sum(PrDens), .groups = 'drop') %>%
    group_by(Dataset) %>%
    mutate(
      Dataset_Total = sum(Total_PrDens),
      Target_Proportion = Total_PrDens / Dataset_Total  
    )
  
  # Merge with original data and scale accordingly
  normalized_spd <- all_spd_df %>%
    left_join(settlement_totals %>% dplyr::select(Dataset, SiteCode, Target_Proportion), 
              by = c("Dataset", "SiteCode")) %>%
    group_by(Dataset, SiteCode) %>%
    mutate(
      Adjusted_PrDens = PrDens * (Target_Proportion / sum(PrDens))  
    ) %>%
    ungroup()
  
  return(normalized_spd)
}


# Combine all SPD datasets
all_spd_df <- bind_rows(
  sampled_spd_total_by_settlement %>% mutate(Dataset = "Unrescaled original"),
  resampled_spd_total_by_settlement %>% mutate(Dataset = "Bootstrapped"),
  weighted_spd_total_by_settlement %>% mutate(Dataset = "Weighted")
)


# ✅ Normalize the SPDs
normalized_spd_df <- normalize_preserve_shape(all_spd_df)
normalized_spd_df$SiteCode <- factor(normalized_spd_df$SiteCode)  


#------------------------------------------------------------------------
#------------  step 2. plot SPD with labels
#------------------------------------------------------------------------

# set color palette
n_settlements <- length(unique(normalized_spd_df$SiteCode))
distinct_colors <- pal_d3("category20")(n_settlements)


# Choose ONE of the two plotting functions below (SiteCode vs. SiteID labels)

#----- 1) plot with label (SiteCode version)
plot_normalized_spds <- function(dataset_name) {
  df_filtered <- normalized_spd_df %>%
    filter(Dataset == dataset_name)
  
  # Find peak positions for each SiteCode for labeling
  label_positions <- df_filtered %>%
    group_by(SiteCode) %>%
    summarise(calBP = calBP[which.max(Adjusted_PrDens)],  
              Adjusted_PrDens = max(Adjusted_PrDens)) 
  
  ggplot(df_filtered, aes(x = calBP, y = Adjusted_PrDens, fill = factor(SiteCode), color = factor(SiteCode))) +
    geom_area(alpha = 0.4, position = "identity") +
    geom_line(linewidth = 0.5) +
    scale_fill_manual(name = "Settlement", values = distinct_colors) +
    scale_color_manual(name = "Settlement", values = distinct_colors) +
    labs(title = paste(dataset_name, " dataset (rescale to 1)"),
         x = "cal BP", y = "Scaled Probability Density") +
    scale_x_reverse(limits = c(2200, 1400)) +
    theme_minimal() +
    theme(text = element_text(family = "noto"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          #legend.spacing.x = unit(0.6, "cm"),
          #legend.key.size = unit(0.6, "cm"),
          #legend.key.width = unit(0.6, "cm"),
          legend.text = element_text(size = 30, face = "bold"),
          legend.title = element_text(size = 30, face = "bold"),
          plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
          axis.text = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size = 25, face = "bold", margin = margin(r = 0)),
          axis.title.x = element_text(size = 25, face = "bold")
    ) +
    geom_text_repel(
      data = label_positions,
      aes(x = calBP, y = Adjusted_PrDens, label = SiteCode),
      size = 11, family = "noto", fontface = "bold",
      box.padding = 0.5, point.padding = 0.5,
      segment.color = "black",           
      segment.size = 1,                  
      max.overlaps = Inf,
      force = 20
    )
}


#----- 2) plot with label (SiteID version)
# Assign IDs 1, 2, 3... in alphanumeric order
site_lookup <- normalized_spd_df %>%
  distinct(SiteCode) %>%
  arrange(SiteCode) %>%            
  mutate(SiteID = row_number())

normalized_spd_df <- normalized_spd_df %>%
  left_join(site_lookup, by = "SiteCode") %>%
  mutate(SiteID_f = factor(SiteID, levels = sort(unique(SiteID))))


site_lookup %>% arrange(SiteID) %>% print(n = Inf)

plot_normalized_spds <- function(dataset_name) {
  df_filtered <- normalized_spd_df %>%
    filter(Dataset == dataset_name)
  
  label_positions <- df_filtered %>%
    group_by(SiteID_f) %>%
    summarise(
      calBP = calBP[which.max(Adjusted_PrDens)],
      Adjusted_PrDens = max(Adjusted_PrDens),
      .groups = "drop"
    )
  
  ggplot(df_filtered, aes(x = calBP, y = Adjusted_PrDens,
                          fill = SiteID_f, color = SiteID_f)) +
    geom_area(alpha = 0.4, position = "identity") +
    geom_line(linewidth = 1.5) +
    scale_fill_manual(name = "Settlement", values = distinct_colors) +
    scale_color_manual(name = "Settlement", values = distinct_colors) +
    labs(title = paste(dataset_name),
         x = "cal BP", y = "Probability Density") +
    scale_x_reverse(limits = c(2200, 1400)) +
    scale_y_continuous(breaks = c(0, 0.001), limits = c(0, NA))+ 
    theme_minimal() +
    theme(text = element_text(family = "noto"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 30, face = "bold"),
          legend.title = element_text(size = 30, face = "bold"),
          plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
          axis.text = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size = 25, face = "bold", margin = margin(r = 0)),
          axis.title.x = element_text(size = 25, face = "bold")
    ) 
  +
    geom_text_repel(
      data = label_positions,
      aes(x = calBP, y = Adjusted_PrDens, label = SiteID_f),
      size = 11, family = "noto", fontface = "bold",
      box.padding = 0.5, point.padding = 0.5,
      segment.color = "black",
      segment.size = 1,
      max.overlaps = Inf,
      force = 20
    )
}




# ✅ Create plots
plot_sampled <- plot_normalized_spds("Unrescaled original")
plot_resampled <- plot_normalized_spds("Bootstrapped")
plot_weighted <- plot_normalized_spds("Weighted")


# Combine plots horizontally and share legend from one plot only
combined_normalized_plot <- plot_sampled + plot_resampled + plot_weighted +
  plot_layout(ncol = 3, guides = "collect") & 
  theme(legend.position = "bottom",
        legend.text = element_text(size = 30, face = "bold"),
        legend.title = element_text(size = 30, face = "bold"),
        legend.key.size = unit(1.5, "cm"),
        legend.key.width = unit(1.5, "cm")) &
  guides(fill = guide_legend(nrow = 2), color = guide_legend(nrow = 2))

print(combined_normalized_plot)
ggsave("Figure10b.png", combined_normalized_plot, width = 30, height = 12, dpi = 300, bg = "white")
#The figures in the manuscript were edited using Adobe Illustrator.


#########################################################################
##################  Figure 9 and 10 (c): Spatio-temporal KDE 
#########################################################################

library(rworldmap)
library(spatstat)
library(sf)
library(raster)
library(gridExtra)
library(rworldmap)
library(rnaturalearth)
library(terra)
# install.packages("spatstat.geom")  # Required for `as.owin()`
# install.packages("maptools")  # Helps with conversion
library(sp)
library(spatstat.geom)
library(gridExtra)

#------------------------------------------------------------------------
#------------  step 1. prepare background map
#------------------------------------------------------------------------

# Set Korea Coordinate System (EPSG:5179)
bng <- CRS("+init=EPSG:5179")
skorea <- getMap(resolution = "high") 
skorea <- skorea[skorea$GEOUNIT == "South Korea", ]  
skorea_sf <- st_as_sf(skorea) 
skorea_sf <- st_transform(skorea_sf, 5179) 
korea_owin <- as.owin(st_geometry(skorea_sf))


# Download and load shapefile for administrative boundaries of Korea
# you can download shapefile here: http://www.gisdeveloper.co.kr/download/admin_shp/sig_20230729.zip
map_korea_total <- st_read("sig_20230729/sig.shp", options = "ENCODING=CP949") 
st_crs(map_korea_total) <- st_crs(5179)  # Korea's Unified Coordinate System
map_korea_total$SIG_KOR_NM <- enc2utf8(map_korea_total$SIG_KOR_NM)



###### filtering region: Choose Yeongsan river basin or Geum river basin

# 1) Yeongsan river basin
map_korea_filtered <- subset(map_korea_total, substr(SIG_CD, 1, 2) == "29" | SIG_CD == "46710")

# 2) Geum river basin 
sig_cd_values <- c("43110", "36110", "43111", "43112", "43113","43114",  "30140", "30170", "30200","30230", "43710") 
map_korea_filtered <- map_korea_total %>%
  filter(SIG_CD %in% sig_cd_values)



# Merge multiple geometries into a single polygon. Convert sf to SpatioPolygons and then to owin
map_korea_filtered <- st_union(map_korea_filtered)
map_korea_sp <- as(map_korea_filtered, "Spatial")
map_korea_spatstat <- as.owin(map_korea_filtered) 


#------------------------------------------------------------------------
#------------ step 2. prepare radiocarbon dates for spatial analysis
#------------------------------------------------------------------------
# Convert sampled_data to an sf object with initial CRS as WGS 84 (EPSG:4326)
sampled_data_sf <- st_as_sf(sampled_data, coords = c("x", "y"), crs = 4326)

# Transform to the same CRS as map_korea_filtered (EPSG:5179)
sampled_data_sf <- st_transform(sampled_data_sf, crs = 5179)
unique_sampled_data <- sampled_data %>% distinct(SiteCode, .keep_all = TRUE)

# Convert the unique data to an sf object
unique_sampled_data_sf <- st_as_sf(
  unique_sampled_data, 
  coords = c("x", "y"), 
  crs = 4326  # WGS 84
) %>%
  st_transform(crs = 5179)

# Plot with non-overlapping text labels for unique points
ggplot() +
  geom_sf(data = map_korea_filtered, fill = "lightgray", color = "black", alpha = 0.7) +
  geom_sf(data = unique_sampled_data_sf, color = "red", size = 3, alpha = 0.8) +
  geom_text_repel(data = unique_sampled_data_sf, aes(label = SiteCode, geometry = geometry), 
                  stat = "sf_coordinates", size = 2, family = "nanumgothic") +
  labs(title = "Sampled Data Points with Labels on Map", x = "Longitude", y = "Latitude") +
  theme_minimal(base_family = "nanumgothic") +
  coord_sf(expand = FALSE) +
  theme(panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        panel.grid.minor = element_line(color = "gray90", linetype = "dashed"))


resampled_output_sf <- st_as_sf(
  resampled_output,
  coords = c("x", "y"),
  crs = 4326  # WGS 84
)

# Transform to the same CRS as map_korea_filtered (EPSG:5179)
resampled_output_sf <- st_transform(resampled_output_sf, crs = 5179)


#------------------------------------------------------------------------
#------------ step 3. draw spatio-temporal KDE with radiocarbon dates
#------------------------------------------------------------------------

# 1) stKDE for published data

# Calibrate resampled data and perform stkde analysis for published data
calibrated_sampled <- calibrate(x = sampled_data$BP,  errors = sampled_data$Error,  
  timeRange = c(2600, 1200), normalised = FALSE)
stkde_sampled <- stkde( x = calibrated_sampled, coords = st_coordinates(sampled_data_sf),  
  win = map_korea_spatstat, sbw=5000, cellres = 50, focalyears = seq(2100, 1400, -100),
  tbw = 50, timeRange = c(2400, 1400), dateRange = range(filtered_sampled_data$BP, na.rm = TRUE))


# Define the file path for the output image
output_file <- "Figure 10c_unrescaled original data 2000-1500.png"
png(filename = output_file, width = 6000, height = 1000, res = 500)
par(mfrow = c(1, 6), mar = c(0.5, 0, 2, 0.5), oma = c(0, 0, 3, 0), font = 2)
#plot(stkde_sampled, 2100, type = "focal", cex.main = 1.8, cex.axis = 1.2, cex.lab = 1.2, main = "Year 2100")
plot(stkde_sampled, 2000, type = "focal", cex.main = 1.8, cex.axis = 1.2, cex.lab = 1.2, main = "Year 2000")
plot(stkde_sampled, 1900, type = "focal", cex.main = 1.8, cex.axis = 1.2, cex.lab = 1.2, main = "Year 1900")
plot(stkde_sampled, 1800, type = "focal", cex.main = 1.8, cex.axis = 1.2, cex.lab = 1.2, main = "Year 1800")
plot(stkde_sampled, 1700, type = "focal", cex.main = 1.8, cex.axis = 1.2, cex.lab = 1.2, main = "Year 1700")
plot(stkde_sampled, 1600, type = "focal", cex.main = 1.8, cex.axis = 1.2, cex.lab = 1.2, main = "Year 1600")
plot(stkde_sampled, 1500, type = "focal", cex.main = 1.8, cex.axis = 1.2, cex.lab = 1.2, main = "Year 1500")
mtext("Unrescaled original", outer = TRUE, cex = 1.2, font = 2, line = 1)
dev.off()


# 2) stKED for bootstrap resampled data

# Calibrate resampled data and perform stkde analysis for resampled data
calibrated_resampled <- calibrate(x = resampled_output$BP, errors = resampled_output$Error, 
  timeRange = c(2400, 1400), normalised = FALSE)
stkde_resampled <- stkde( x = calibrated_resampled, coords = st_coordinates(resampled_output_sf), 
  win = map_korea_spatstat, sbw = 5000, cellres = 50, focalyears = seq(2100, 1400, -100),
  tbw = 50, timeRange = c(2400, 1400), dateRange = range(resampled_output$BP, na.rm = TRUE))



# Plot the density for each focal year
output_file_resampled <- "Figure 10c_Bootstrapped data 2100-1600.png"
png(filename = output_file_resampled, width = 6000, height = 1000, res = 500)
par(mfrow = c(1, 6), mar = c(0.5, 0, 2, 0.5), oma = c(0, 0, 3, 0), font = 2)
plot(stkde_resampled, 2100, type = "focal", cex.main = 1.8, cex.axis = 1.2, cex.lab = 1.2, main = "Year 2100")
plot(stkde_resampled, 2000, type = "focal", cex.main = 1.8, cex.axis = 1.2, cex.lab = 1.2, main = "Year 2000")
plot(stkde_resampled, 1900, type = "focal", cex.main = 1.8, cex.axis = 1.2, cex.lab = 1.2, main = "Year 1900")
plot(stkde_resampled, 1800, type = "focal", cex.main = 1.8, cex.axis = 1.2, cex.lab = 1.2, main = "Year 1800")
plot(stkde_resampled, 1700, type = "focal", cex.main = 1.8, cex.axis = 1.2, cex.lab = 1.2, main = "Year 1700")
plot(stkde_resampled, 1600, type = "focal", cex.main = 1.8, cex.axis = 1.2, cex.lab = 1.2, main = "Year 1600")
#plot(stkde_resampled, 1500, type = "focal", cex.main = 1.8, cex.axis = 1.2, cex.lab = 1.2, main = "Year 1500")
mtext("Bootstrapped", outer = TRUE, cex = 1.3, font = 2, line = 1)
dev.off()



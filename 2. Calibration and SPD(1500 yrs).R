
# ------------------------------------------------
# Load required packages (install if needed)
# ------------------------------------------------
library(tidyverse)
library(rcarbon)
library(truncnorm)
library(ggplot2)
library(readr)
library(gridExtra)
library(corrplot)
library(RColorBrewer)
library(viridis)
library(reshape2)
library(ggpubr)
library(readr)


getwd()

# ---------------------------------------------
# Step 1: Load the dataset (by percentage. e.g., data for maximum sampling fraction 20, 30, 40)
# ---------------------------------------------
original_data <- read.csv("")
sampled_data <- read.csv("")
resampled_data <- read.csvread.csv("")
resampled_30_data <- read.csvread.csv("")

# Add a fixed error column (in this case, 40)
original_data$Error <- 40
sampled_data$Error <- 40
resampled_data$Error <- 40
resampled_30_data$Error <- 40



base_path<-("")


# Create new directories
dir.create(file.path(base_path, "CalData"), recursive = TRUE)

dir.create(file.path(base_path, "CalPlots"), recursive = TRUE)



# Function to calibrate dates
calibrate_dates <- function(settlement_data) {
  calibrated <- calibrate(settlement_data$Date_BP, errors = rep(40, nrow(settlement_data)), normalised = TRUE)
  return(calibrated)
}


extract_spd_info <- function(cal_date) {
  if(inherits(cal_date, "CalDates")) {
    spd_result <- spd(cal_date, 
                      timeRange = c(5500, 3000), 
                      spdnormalised = FALSE)  # Set spdnormalised to FALSE
    return(data.frame(
      calBP = spd_result$grid$calBP,
      PrDens = spd_result$grid$PrDen
    ))
  }
  return(NULL)
}


# Function to create plotting dataframe for each dataset
create_plot_df <- function(spd_results, dataset_name) {
  do.call(rbind, lapply(1:nrow(spd_results), function(i) {
    spd_info <- spd_results$SPD[[i]]
    if(!is.null(spd_info)) {
      spd_info$Settlement <- spd_results$Settlement[i]
      spd_info$Dataset <- dataset_name
      return(spd_info)
    }
    return(NULL)
  }))
}

# Function to normalize SPDs while preserving shapes (total = 1)
normalize_preserve_shape <- function(all_plot_df) {
  # First calculate total contribution of each settlement
  settlement_totals <- all_plot_df %>%
    group_by(Dataset, Settlement) %>%
    summarise(Total_PrDens = sum(PrDens), .groups = 'drop') %>%
    group_by(Dataset) %>%
    mutate(
      Dataset_Total = sum(Total_PrDens),
      Target_Proportion = round(Total_PrDens/Dataset_Total,3) # Now proportion will be between 0 and 1
    )
  
  # Join with original data and scale
  normalized_spd <- all_plot_df %>%
    left_join(
      settlement_totals %>% select(Dataset, Settlement, Target_Proportion),
      by = c("Dataset", "Settlement")
    ) %>%
    group_by(Dataset, Settlement) %>%
    mutate(
      # Scale each settlement's SPD while preserving its shape
      Adjusted_PrDens = PrDens * (Target_Proportion / sum(PrDens))
    ) %>%
    ungroup()
  
  return(normalized_spd)
}



# Function to process one combination
process_combination <- function(seed_s, seed_rs, settlement_dist, bp_dist) {
  

  # Create combination name for file naming
  combo_name <- sprintf("S%d_RS%d_%s_%s", 
                        seed_s, seed_rs, settlement_dist, bp_dist)
  
  # Filter data for this combination
  original_filtered <- original_data %>%
    filter(Settlement_Dist == settlement_dist, BP_Dist == bp_dist)
  
  sampled_filtered <- sampled_data %>%
    filter(SeedS == seed_s, 
           Settlement_Dist == settlement_dist, 
           BP_Dist == bp_dist)
  
  resampled_filtered <- resampled_data %>%
    filter(SeedS == seed_s, 
           SeedRS == seed_rs, 
           Settlement_Dist == settlement_dist, 
           BP_Dist == bp_dist)
  
  resampled_30_filtered <- resampled_30_data %>%
    filter(SeedS == seed_s, 
           SeedRS == seed_rs, 
           Settlement_Dist == settlement_dist, 
           BP_Dist == bp_dist)
  
  # Calculate SPDs
  original_spd_results <- original_filtered %>%
    group_by(Settlement) %>%
    summarise(Calibrated_Dates = list(calibrate_dates(cur_data())), .groups = 'drop') %>%
    mutate(SPD = lapply(Calibrated_Dates, function(x) extract_spd_info(x)))
  
  sampled_spd_results <- sampled_filtered %>%
    group_by(SeedS, Settlement) %>%
    summarise(Calibrated_Dates = list(calibrate_dates(cur_data())), .groups = 'drop') %>%
    mutate(SPD = lapply(Calibrated_Dates, function(x) extract_spd_info(x)))
  
  weighted_sampled_spd_results <- sampled_filtered %>%
    group_by(SeedS, Settlement) %>%
    summarise(
      Calibrated_Dates = list(calibrate_dates(cur_data())),
      Total_Houses = first(Total_Houses),
      Sampled_Houses = n(),
      Weight = Total_Houses / Sampled_Houses,
      .groups = 'drop'
    ) %>%
    mutate(
      SPD = mapply(function(x, w) {
        spd_info <- extract_spd_info(x)
        if(!is.null(spd_info)) {
          spd_info$PrDens <- spd_info$PrDens * w
        }
        return(spd_info)
      }, Calibrated_Dates, Weight, SIMPLIFY = FALSE)
    )
  
  
  resampled_30_spd_results <- resampled_30_filtered %>%
    group_by(SeedS, SeedRS, Settlement) %>%
    summarise(Calibrated_Dates = list(calibrate_dates(cur_data())), .groups = 'drop') %>%
    mutate(SPD = lapply(Calibrated_Dates, function(x) extract_spd_info(x)))
  
  # Create plotting dataframes
  original_plot_df <- create_plot_df(original_spd_results, "Original")
  sampled_plot_df <- create_plot_df(sampled_spd_results, "Sampled")
  weighted_sampled_plot_df <- create_plot_df(weighted_sampled_spd_results, "Weighted Sampled")
  resampled_30plot_df <- create_plot_df(resampled_30_spd_results, "Resampled30")

  # Combine all plotting dataframes
  all_plot_df <- rbind(
    original_plot_df,
    sampled_plot_df,
    weighted_sampled_plot_df,
    resampled_30plot_df
  )
  
  
  # Create and save normalized plots
  normalized_results <- normalize_preserve_shape(all_plot_df)
 
  # Save all results
  write.csv(all_plot_df, 
            file.path(base_path, "CalData",
                      paste0("all_plot_df_", combo_name, ".csv")), 
            row.names = FALSE)
  
  write.csv(normalized_results,
            file.path(base_path, "CalData",
                      paste0("normalized_", combo_name, ".csv")),
            row.names = FALSE)


  return(list(
    all_plot_df = all_plot_df,
    normalized_results = normalized_results
  ))
}

# Main function to process all combinations
process_all_combinations <- function(SeedS, SeedRS, Settlement_Dist, BP_Dist) {
  results_list <- list()
  
  for(s in SeedS) {
    for(rs in SeedRS) {
      for(sd in Settlement_Dist) {
        for(bp in BP_Dist) {
          cat(sprintf("\nProcessing: SeedS=%d, SeedRS=%d, %s, %s\n", 
                      s, rs, sd, bp))
          
          result <- tryCatch({
            process_combination(s, rs, sd, bp)
          }, error = function(e) {
            warning(sprintf("Error in combination S%d_RS%d_%s_%s: %s", 
                            s, rs, sd, bp, e$message))
            return(NULL)
          })
          
          if(!is.null(result)) {
            results_list[[length(results_list) + 1]] <- result
            cat("Successfully processed combination\n")
          }
        }
      }
    }
  }
  
  return(results_list)
}

# Process all combinations
results <- process_all_combinations(
  SeedS = 31:50,
  SeedRS = c(1),#,2),
  Settlement_Dist = c("normal","power_law"),   
  BP_Dist = c("uniform", "normal","skewed")) 







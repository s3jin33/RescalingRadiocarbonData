library(tidyverse)
library(rcarbon)
library(truncnorm)
library(ggplot2)
library(readr)
library(zoo)


#Define Global Seed
set.seed(42)

#Define Original Seed
seeds_og <- c(42)  
seeds_sp <- c(1:50)
seeds_resample <- c(1)
n_settlements <- 10
total_houses <- 1000



# Define all distributions
settlement_distributions <- c("normal", "power_law")
bp_distributions <- c("uniform", "normal", "skewed")


# Create output directory
if (!dir.exists("Plots")) dir.create("Plots")
if (!dir.exists("Data")) dir.create("Data")

############Step 1: Functions to Generate Settlement Size Distributions##############################

###Setting minmum number of houses. In this case, 5
adjust_distribution <- function(dist, total_houses) {
  # Ensure each settlement has at least one house
  dist[dist < 5] <- 5
  
  # Distribute remaining houses to match total_houses 
  difference <- total_houses - sum(dist)
  
  while (difference != 0) {
    for (i in seq_along(dist)) {
      if (difference == 0) break
      if (difference > 0) {
        dist[i] <- dist[i] + 1  # Add houses if total is short
        difference <- difference - 1
      } else if (difference < 0 && dist[i] > 5) {
        dist[i] <- dist[i] - 1  # Subtract if over
        difference <- difference + 1
      }
    }
  }
  
  return(dist)
}


generate_normal <- function(n_settlements, total_houses) {
  dist <- round(rnorm(n_settlements, mean = total_houses / n_settlements, sd = total_houses / (n_settlements * 2.5)))
  adjust_distribution(dist, total_houses)
}


generate_power_law <- function(n_settlements, total_houses) {
  # Generate random weights based on a power-law distribution
  random_weights <- (1:n_settlements)^(-1)  # Base power-law weights
  random_weights <- random_weights * runif(n_settlements, min = 0.8, max = 1.2)  # Add randomness
  
  # Scale weights to sum to total_houses
  dist <- round(random_weights / sum(random_weights) * total_houses)
  
  # Ensure the total matches total_houses exactly
  difference <- total_houses - sum(dist)
  while (difference != 0) {
    for (i in seq_along(dist)) {
      if (difference == 0) break
      if (difference > 0) {
        dist[i] <- dist[i] + 1
        difference <- difference - 1
      } else if (difference < 0 && dist[i] > 1) {
        dist[i] <- dist[i] - 1
        difference <- difference + 1
      }
    }
  }
  
  return(dist)
}



####### Step 2: Functions to Generate BP Date Distributions#########################

# Generate valid start and end dates for each settlement
generate_life_cycle_dates <- function() {
  while (TRUE) {
    # Generate the start date (older, larger BP)
    start_date <- sample(4500:3100, 1)
    
    # Generate the end date (younger, smaller BP)
    end_date_min <- max(3000, start_date - 500)  # At most 500 years younger
    end_date_max <- min(4500, start_date - 100)  # At least 50 years younger
    
    # Ensure a valid range for the end date
    if (end_date_min <= end_date_max) {
      if (end_date_min == end_date_max) {
        end_date <- end_date_min
      } else {
        end_date <- sample(end_date_min:end_date_max, 1)
      }
      return(c(start_date, end_date))
    }
  }
}


# Function to generate BP dates within a given start and end date
generate_uniform_dates <- function(n, start_date, end_date) {
  # Ensure the range is strictly between start_date (larger BP) and end_date (smaller BP)
  sample(end_date:start_date, n, replace = TRUE)
}

generate_normal_dates <- function(n, start_date, end_date) {
  mean_date <- (start_date + end_date) / 2
  sd_date <- (start_date - end_date) / 4
  round(pmin(pmax(rtruncnorm(n, a = end_date, b = start_date, mean = mean_date, sd = sd_date), end_date), start_date))
}

generate_skewed_dates <- function(n, start_date, end_date) {
  round(pmin(pmax(end_date + rbeta(n, 2, 5) * (start_date - end_date), end_date), start_date))
}


validate_settlement_data <- function(data) {
  # Validate start and end dates
  invalid_start_end <- data %>%
    filter(Start_Date <= End_Date)  # Start date (larger BP) must be greater than End date (smaller BP)
  
  if (nrow(invalid_start_end) > 0) {
    warning("❌ Invalid Start and End dates detected.")
    print(invalid_start_end)
  } else {
    print("✅ All Start and End dates are valid.")
    invalid_start_end <- NULL  # Explicitly set to NULL if valid
  }
  
  # Validate life spans
  invalid_life_span <- data %>%
    group_by(Settlement) %>%
    summarise(Life_Span = max(Start_Date) - min(End_Date)) %>%
    filter(Life_Span < 100 | Life_Span > 500)
  
  if (nrow(invalid_life_span) > 0) {
    warning("❌ Invalid life spans detected.")
    print(invalid_life_span)
  } else {
    print("✅ All life spans are valid.")
    invalid_life_span <- NULL  # Explicitly set to NULL if valid
  }
  
  # Return both invalid cases for review
  list(
    invalid_start_end = invalid_start_end,
    invalid_life_span = invalid_life_span
  )
}



assign_bp_dates_varying_life_cycle <- function(settlement_dist, bp_dist_name) {
  settlement_dist <- settlement_dist[settlement_dist > 0]
  settlement_dist[settlement_dist < 5] <- 5
  
  total_houses <- sum(settlement_dist)
  settlement_start_end_dates <- t(replicate(length(settlement_dist), generate_life_cycle_dates()))
  
  settlement_details <- data.frame(Settlement = integer(),
                                   Total_Houses = integer(),
                                   House_ID = integer(),
                                   Start_Date = integer(),
                                   End_Date = integer(),
                                   Life_Span = integer(),
                                   Date_BP = integer())
  
  for (sett in 1:length(settlement_dist)) {
    start_date <- settlement_start_end_dates[sett, 1]
    end_date <- settlement_start_end_dates[sett, 2]
    
    # Generate BP dates and ensure bounds are included
    settlement_bpd <- switch(bp_dist_name,
                             "uniform" = generate_uniform_dates(settlement_dist[sett] - 2, start_date, end_date),
                             "normal" = generate_normal_dates(settlement_dist[sett] - 2, start_date, end_date),
                             "skewed" = generate_skewed_dates(settlement_dist[sett] - 2, start_date, end_date))
    
    # Force inclusion of start_date and end_date
    settlement_bpd <- c(start_date, end_date, settlement_bpd)
    
    # Shuffle to distribute start and end dates randomly within the settlement
    settlement_bpd <- sample(settlement_bpd)
    
    settlement_details <- rbind(settlement_details,
                                data.frame(Settlement = rep(sett, settlement_dist[sett]),
                                           Total_Houses = rep(settlement_dist[sett], settlement_dist[sett]),
                                           House_ID = 1:settlement_dist[sett],
                                           Start_Date = rep(start_date, settlement_dist[sett]),
                                           End_Date = rep(end_date, settlement_dist[sett]),
                                           Life_Span = rep(start_date - end_date, settlement_dist[sett]),
                                           Date_BP = settlement_bpd))
  }
  
  return(settlement_details)
}


##########Step 3: Sampling and Resampling Functions #################################

sample_houses <- function(data, sampling_pct, settlement_dist_name, bp_dist_name, seed_s) {
  # Create a unique seed for this sampling step
  base_seed <- sum(utf8ToInt(paste0(settlement_dist_name, bp_dist_name)))
  combined_seed <- base_seed + seed_s + sampling_pct  # Incorporate sampling percentage in the seed
  set.seed(combined_seed)
  
  # Convert percentage to a decimal (e.g., 20 -> 0.20)
  sampling_proportion <- sampling_pct / 100
  
  # Perform stratified sampling by Settlement
  sampled_data <- data %>%
    group_by(Settlement) %>%
    group_modify(~ {
      n_total <- nrow(.x)
      
      # Calculate maximum allowed sample size based on the percentage.
      # Using floor() here so that for 100 houses and 20% sampling, max_sample is 20.
      max_sample <- floor(n_total * sampling_proportion)
      
      # Ensure that at least 5 houses are sampled.
      # If max_sample is less than 5 (unlikely if each settlement has at least 5 houses), use 5.
      if (max_sample < 5) {
        final_sample_size <- 5
      } else {
        final_sample_size <- sample(5:max_sample, 1)
      }
      
      # Randomly sample the rows using the chosen sample size
      .x %>%
        mutate(rank = runif(n_total)) %>%
        arrange(rank) %>%
        slice_head(n = final_sample_size)
    }) %>%
    ungroup()
  
  sampled_counts <- sampled_data %>%
    group_by(Settlement) %>%
    summarise(Sampled_Houses = n(), .groups = "drop")
  
  sampled_data <- sampled_data %>%
    left_join(sampled_counts, by = "Settlement")
  
  return(sampled_data)
}


resample_data <- function(sampled_data, original_data, n_resamples, settlement_dist_name, bp_dist_name, seed_rs) {
  # Create a unique seed for each distribution combination and resample seed
  seed_value <- sum(utf8ToInt(paste0(settlement_dist_name, bp_dist_name))) + seed_rs
  set.seed(seed_value)

  print("=== Debugging Resample Data ===")
  print(paste("Total rows in sampled_data:", nrow(sampled_data)))
  print(paste("Total rows in original_data:", nrow(original_data)))
  
  resampled_data <- data.frame()
  
  for (resample_iter in 1:n_resamples) {
    print(paste("Starting resampling iteration:", resample_iter))
    
    for (sett in unique(original_data$Settlement)) {
      n_houses <- nrow(subset(original_data, Settlement == sett))
      sampled_subset <- subset(sampled_data, Settlement == sett)
      
      print(paste("Settlement:", sett, "Original Houses:", n_houses, "Sampled Houses:", nrow(sampled_subset)))
      
      if (nrow(sampled_subset) == 0) {
        warning(paste("No sampled data for Settlement:", sett))
        next
      }
      
      resampled_indices <- sample(seq_len(nrow(sampled_subset)), n_houses, replace = TRUE)
      resampled_subset <- sampled_subset[resampled_indices, ]
      
      resampled_subset <- resampled_subset %>%
        mutate(
          Resample_Iteration = resample_iter,
          Sampling_Num = seq_len(nrow(resampled_subset))
        )
      
      resampled_data <- rbind(resampled_data, resampled_subset)
    }
    
    print(paste("Completed resampling iteration:", resample_iter))
  }
  
  print("Final resampled data distribution:")
  print(resampled_data %>% group_by(Settlement) %>% summarise(Resampled_Houses = n()))
  
  return(resampled_data)
}



###################Step 4: Main Loop with Debugging and Plotting#######################

# Initialize data frames for original data
combined_original_data <- data.frame()

# First loop: Generate original datasets
for (seeds_original in seeds_og) {
  set.seed(seeds_original)  # Set the seed for generating original data
  for (settlement_dist_name in settlement_distributions) {
    set.seed(42)
    print(paste("Processing Settlement Distribution:", settlement_dist_name))
    
    # Step 1: Generate original settlement data
    settlement_dist <- get(paste0("generate_", settlement_dist_name))(n_settlements, total_houses)
    
    # Verify total houses
    current_total <- sum(settlement_dist)
    print(paste("Settlement distribution total houses:", current_total))
    stopifnot(current_total == total_houses)
    
    for (bp_dist_name in bp_distributions) {
      print(paste("Processing BP Distribution:", bp_dist_name))
      
      original_data <- assign_bp_dates_varying_life_cycle(settlement_dist, bp_dist_name)
      original_data <- original_data %>%
        mutate(SeedOG = seeds_original,
               Settlement_Dist = settlement_dist_name,
               BP_Dist=bp_dist_name)
      
      combined_original_data <- bind_rows(combined_original_data, original_data)
      
      # Verify total houses in original data
      original_total <- nrow(original_data)
      print(paste("Original data total houses:", original_total))
      stopifnot(original_total == total_houses)
      
      # Save original and sampled data for reference
      output_dir <- paste0("Data/Original")
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      write_csv(original_data, paste0(output_dir, "/Original_", settlement_dist_name, "_",bp_dist_name,"_","OG", seeds_original, ".csv"))


      # Check if validation passed
      validation_results <- validate_settlement_data(original_data)
      if ((!is.null(validation_results$invalid_start_end) && nrow(validation_results$invalid_start_end) > 0) || 
          (!is.null(validation_results$invalid_life_span) && nrow(validation_results$invalid_life_span) > 0)) {
        print("❌ Validation errors detected. Skipping further processing.")
        next
      }
    }
  }
}


saveRDS(combined_original_data, file.path(output_dir, "combined_original_data.rds"))

write.csv(combined_original_data, file.path(output_dir,"combined_original_data.csv"))






#############Part 2: Sampling, resampling, and plotting#######

#Retrieving Original Data
combined_original_data <- readRDS(file.path(output_dir,"combined_original_data.rds"))


#Change maximum sampling fractions. e.g., 20, 30, 40
sampling_pct <- c(40)

####Reset dataframes
combined_sampled_data <- data.frame()
combined_resampled_data <- data.frame()
combined_resampled30_data <- data.frame()
combined_spd_comparison_norm <- data.frame()
summary_rows <- list()


######## Second loop: Sampling, resampling, and plotting#####
for (seeds_original in seeds_og) {
  for (settlement_dist_name in settlement_distributions) {
    for (bp_dist_name in bp_distributions) {
      # Get the corresponding original data
      original_data <- combined_original_data %>%
        filter(SeedOG == seeds_original,
               Settlement_Dist == settlement_dist_name,
               BP_Dist == bp_dist_name)

      # Plot 1: Frequency of Settlements by House Size
      plots_dir <- paste0("Plots/OG Seed ",seeds_original, "_sampling ", sampling_pct)
      if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

      freq_original <- original_data %>%
        group_by(Settlement) %>%
        summarise(Houses = n(), .groups = "drop")

      png(filename = paste0(plots_dir, "/Frequency_House_Size_", settlement_dist_name, "_", bp_dist_name, "_OG", seeds_original, ".png"),
          width = 800, height = 600)
      print(
        ggplot(freq_original, aes(x = Houses)) +
          geom_histogram(binwidth = 10, fill = "steelblue", color = "black", boundary = 0) +
          labs(title = paste("Frequency of Settlements by House Size\n", settlement_dist_name, " Distribution (OG: ", seeds_original, ")", sep = ""),
               x = "Number of Houses per Settlement", y = "Frequency") +
          theme_minimal() +
          theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
      )
      dev.off()

      # Step 2: Sample houses from the original data
      for (seeds_sample in seeds_sp) {
        set.seed(seeds_sample)  # Set the seed for sampling
        print(paste("Processing Seed Sample:", seeds_sample))
        sampled_data <- sample_houses(original_data, 
                                      sampling_pct = sampling_pct, 
                                      settlement_dist_name = settlement_dist_name,
                                      bp_dist_name = bp_dist_name,
                                      seed_s = seeds_sample)
        
        sampled_data <- sampled_data %>%
          mutate(SeedS = seeds_sample)  # Add sampled seed for tracking
        
        # Combine data into accumulators
        
        combined_sampled_data <- bind_rows(combined_sampled_data, sampled_data)
        

        sample_output_dir <- paste0("Data/OG Seed ",seeds_original,"_sampling ", sampling_pct,"/S Seed_",seeds_sample)
        if (!dir.exists(sample_output_dir)) dir.create(sample_output_dir, recursive = TRUE)
        write_csv(sampled_data, paste0(sample_output_dir, "/Sampled_", settlement_dist_name, "_",bp_dist_name,"_OG", seeds_original, "_S", seeds_sample,"_P", sampling_pct,".csv"))


        # Step 3:Perform resampling with different seeds
        for (i in seq_along(seeds_resample)) {
          set.seed(seeds_resample[i])
          print(paste("Processing Seed Resample:", seeds_resample))
          
          # Step 3-1: Resample once
          resampled_1 <- resample_data(sampled_data, original_data, 
                                       n_resamples = 1,
                                       settlement_dist_name = settlement_dist_name,
                                       bp_dist_name = bp_dist_name,
                                       seed_rs = seeds_resample[i])
          
                    # Combine resampled 10 data across seeds
          resampled_1 <- resampled_1 %>%
            mutate(SeedRS = seeds_resample[i])  # Add resample seed information
          combined_resampled_data <- bind_rows(combined_resampled_data, resampled_1)
          
          # Step 3-2: Resample 30. 
          resampled_30 <- resample_data(sampled_data, original_data, 
                                        n_resamples = 30,
                                        settlement_dist_name = settlement_dist_name,
                                        bp_dist_name = bp_dist_name,
                                        seed_rs = seeds_resample[i])

          # Combine resampled 30 data across seeds
          resampled_30 <- resampled_30 %>%
            mutate(SeedRS = seeds_resample[i])  # Add resample seed information
          combined_resampled30_data <- bind_rows(combined_resampled30_data, resampled_30)
          
          #Save resampled outputs
         resample_output_dir <-paste0("Data/OG Seed ",seeds_original,"_sampling ", sampling_pct,"/S Seed_",seeds_sample)
         if (!dir.exists(resample_output_dir)) dir.create(resample_output_dir, recursive = TRUE)
         write_csv(resampled_1, paste0(resample_output_dir, "/Resampled_1_", settlement_dist_name, "_", bp_dist_name,"_OG", seeds_original, "_S", seeds_sample, "_RS",seeds_resample[i],"_P", sampling_pct,".csv"))
         write_csv(resampled_30, paste0(resample_output_dir, "/Resampled_30_",settlement_dist_name,"_", bp_dist_name,"_OG", seeds_original, "_S", seeds_sample, "_RS",seeds_resample[i],"_P", sampling_pct,".csv"))


          # Plot 2: Overlay Histogram: original, sampled, resampled 1.
          png(filename = paste0(plots_dir, "/Overlay_BP_Histogram_", settlement_dist_name, "_", bp_dist_name,"_OG", seeds_original,"_S", seeds_sample, "_RS",seeds_resample[i],"_P", sampling_pct,".png"),
              width = 800, height = 600)
          print(
            ggplot() +
              geom_histogram(data = original_data, aes(x = Date_BP, fill = "Original"), binwidth = 25, alpha = 0.5) +
              geom_histogram(data = sampled_data, aes(x = Date_BP, fill = "Sampled"), binwidth = 25, alpha = 0.7) +
              geom_histogram(data = resampled_1, aes(x = Date_BP, fill = "Resampled 1"), binwidth = 25, alpha = 0.6) +
              scale_fill_manual(name = "Data Type", values = c("Original" = "steelblue", "Sampled" = "forestgreen", "Resampled 1" = "orange")) +
              scale_x_reverse(expand = expansion(mult = c(0, 0.05))) +
              labs(title = paste("Overlay Histogram (", settlement_dist_name, "-", bp_dist_name, "-", "OG: ", seeds_original,"_S: ", seeds_sample, "_RS: ", seeds_resample[i], "_PCT: ", sampling_pct,"%)", sep = ""),
                   x = "Uncalibrated BP Date", y = "Count") +
              theme_minimal() +
              theme(legend.position = "top", plot.title = element_text(size = 18, face = "bold", hjust = 0.5))

          )
          dev.off()


          # Plot 3: Side-by-Side Histogram: original, sampled, resampled, resampled 30.
          png(filename = paste0(plots_dir, "/Side_by_Side_BP_Comparison_", settlement_dist_name,"_", bp_dist_name,"_OG", seeds_original,"_S", seeds_sample,  "_RS",seeds_resample[i], "_P", sampling_pct,".png"),width = 1800, height = 600)

          # Combine data and explicitly order the Type levels
          combined_bp_data <- bind_rows(
            original_data %>% mutate(Type = "Original"),
            sampled_data %>% mutate(Type = "Sampled"),
            resampled_1 %>% mutate(Type = "Resampled 1")
          ) %>%
            mutate(Type = factor(Type, levels = c("Original", "Sampled", "Resampled 1"))) # Explicit ordering

          print(
            ggplot(combined_bp_data, aes(x = Date_BP, fill = Type)) +
              geom_histogram(binwidth = 25, color = "black", alpha = 0.7, position = "identity") +
              scale_x_reverse(expand = expansion(mult = c(0, 0.05))) +
              facet_grid(Type ~ Settlement) +
              scale_fill_manual(values = c("Original" = "steelblue", "Sampled" = "forestgreen", "Resampled 1" = "orange")) +
              labs(title = paste("Side-by-Side BP Date Comparison by Settlement (OG: ", seeds_original,"_S: ", seeds_sample,"_RS: ",seeds_resample[i],"_PCT: ", sampling_pct,"%)", sep = ""),
                   x = "Uncalibrated BP Date", y = "Number of Houses") +
              theme_minimal() +
              theme(legend.position = "none", plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
          )
          dev.off()

          # Plot 3-2: Side-by-Side Histogram(with resampel10)
          png(filename = paste0(plots_dir, "/Side_by_Side_BP_Comparison2_", settlement_dist_name, "_", bp_dist_name, "_OG", seeds_original, "_S", seeds_sample,"_RS", seeds_resample[i],"_P", sampling_pct,".png"),width = 1800, height = 600)

          # Combine data and explicitly order the Type levels
          combined_bp_data <- bind_rows(
            original_data %>% mutate(Type = "Original"),
            sampled_data %>% mutate(Type = "Sampled"),
            resampled_1 %>% mutate(Type = "Resampled 1"),
            resampled_30 %>% mutate(Type = "Resampled 30")
          ) %>%
            mutate(Type = factor(Type, levels = c("Original", "Sampled", "Resampled 1", "Resampled 30"))) # Explicit ordering

          print(
            ggplot(combined_bp_data, aes(x = Date_BP, fill = Type)) +
              geom_histogram(binwidth = 25, color = "black", alpha = 0.7, position = "identity") +
              scale_x_reverse(expand = expansion(mult = c(0, 0.05))) +
              facet_grid(Type ~ Settlement, scales = "free_y") +  # Use free y scales for better control of large datasets
              scale_fill_manual(values = c("Original" = "steelblue", "Sampled" = "forestgreen", "Resampled 1" = "orange", "Resampled 30" = "purple")) +
              labs(title = paste("Side-by-Side BP Date Comparison by Settlement (OG: ", seeds_original,"_S: ", seeds_sample, "_RS: ", seeds_resample[i], "_PCT: ", sampling_pct,"%)", sep = ""),
                   x = "Uncalibrated BP Date", y = "Number of Houses") +
              theme_minimal() +
              theme(legend.position = "none", plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
          )
          dev.off()



          # Step 4: Record summary information
          for (settlement_id in unique(original_data$Settlement)) {
            # Filter subsets
            original_subset <- original_data[original_data$Settlement == settlement_id, ]
            sampled_subset <- sampled_data[sampled_data$Settlement == settlement_id, ]
            resampled_1_subset <- resampled_1[resampled_1$Settlement == settlement_id, ]
            resampled_30_subset <- resampled_30[resampled_30$Settlement == settlement_id, ]

            # Debugging: Print filtered subsets for inspection
            cat("Settlement ID:", settlement_id, "\n")
            cat("Original subset rows:", nrow(original_subset), "\n")
            cat("Sampled subset rows:", nrow(sampled_subset), "\n")
            cat("Resampled 1 subset rows:", nrow(resampled_1_subset), "\n")
            cat("Resampled 30 subset rows:", nrow(resampled_30_subset), "\n\n")

            # Check if Date_BP exists in subsets
            if (!"Date_BP" %in% colnames(original_subset)) print("Missing Date_BP in original_subset")
            if (!"Date_BP" %in% colnames(sampled_subset)) print("Missing Date_BP in sampled_subset")
            if (!"Date_BP" %in% colnames(resampled_1_subset)) print("Missing Date_BP in resampled_1_subset")
            if (!"Date_BP" %in% colnames(resampled_30_subset)) print("Missing Date_BP in resampled_30_subset")

            # Append to summary
            summary_rows[[length(summary_rows) + 1]] <- data.frame(
              SeedOG = seeds_original,
              SeedS = seeds_sample,
              SeedRS = seeds_resample,
              Settlement_Dist = settlement_dist_name,
              BP_Dist = bp_dist_name,
              Settlement = settlement_id,
              Total_Houses = nrow(original_subset),
              Start_Date_Original = if (nrow(original_subset) > 0) max(original_subset$Date_BP) else NA,
              End_Date_Original = if (nrow(original_subset) > 0) min(original_subset$Date_BP) else NA,
              Life_Span_Original = if (nrow(original_subset) > 0) max(original_subset$Date_BP) - min(original_subset$Date_BP) else NA,
              Sampled_Houses = nrow(sampled_subset),
              Start_Date_Sampled = if (nrow(sampled_subset) > 0) max(sampled_subset$Date_BP) else NA,
              End_Date_Sampled = if (nrow(sampled_subset) > 0) min(sampled_subset$Date_BP) else NA,
              Life_Span_Sampled = if (nrow(sampled_subset) > 0) max(sampled_subset$Date_BP) - min(sampled_subset$Date_BP) else NA,
              Start_Date_Resampled_1 = if (nrow(resampled_1_subset) > 0) max(resampled_1_subset$Date_BP) else NA,
              End_Date_Resampled_1 = if (nrow(resampled_1_subset) > 0) min(resampled_1_subset$Date_BP) else NA,
              Life_Span_Resampled_1 = if (nrow(resampled_1_subset) > 0) max(resampled_1_subset$Date_BP) - min(resampled_1_subset$Date_BP) else NA
            )


        }
      }
    }
  }
}
}


# Save combined original and sampled data

summary_output_dir <- paste0("Data/Summary")
if (!dir.exists(summary_output_dir)) dir.create(summary_output_dir, recursive = TRUE)

summary_data <- bind_rows(summary_rows)
write_csv(summary_data, file.path(summary_output_dir, paste0("Summary_Settlements_OG", seeds_original,"_P", sampling_pct,".csv")))

write_csv(combined_sampled_data, file.path(summary_output_dir,  paste0("Combined_Sampled_Data_OG", seeds_original,"_P", sampling_pct,".csv")))

# Export combined resampled data
write_csv(combined_resampled_data, paste0(summary_output_dir, "/Combined_Resampled_1_Data_OG", seeds_original,"_P", sampling_pct,".csv"))

write_csv(combined_resampled30_data, paste0(summary_output_dir, "/Combined_Resampled_30_Data_OG", seeds_original,"_P", sampling_pct,".csv"))



################ Verification with standard data################

# Function to verify and compare datasets
verify_multiple_datasets <- function(file_paths) {
  # Validate input paths
  existing_paths <- file_paths[file.exists(file_paths)]
  if (length(existing_paths) == 0) {
    stop("No valid file paths found")
  }
  
  # Print files being compared
  cat("\n=== Files Being Compared ===\n")
  for (path in existing_paths) {
    cat("- ", basename(path), "\n")
  }
  
  # Read datasets
  datasets <- list()
  for (path in existing_paths) {
    tryCatch({
      datasets[[basename(path)]] <- read_csv(path)
      cat("✓ Successfully read:", basename(path), "\n")
    }, error = function(e) {
      warning("❌ Failed to read ", basename(path), ": ", e$message)
    })
  }
  
  # Validate number of datasets
  n_datasets <- length(datasets)
  if (n_datasets < 2) {
    stop("Need at least 2 valid datasets to compare. Only found ", n_datasets)
  }
  
  # Initialize results dataframe
  comparison_results <- data.frame(
    Dataset1 = character(),
    Dataset2 = character(),
    Identical = logical(),
    Differences = character(),
    stringsAsFactors = FALSE
  )
  
  # Compare datasets
  dataset_names <- names(datasets)
  cat("\n=== Comparing Datasets ===\n")
  for (i in 1:(length(dataset_names)-1)) {
    for (j in (i+1):length(dataset_names)) {
      tryCatch({
        name1 <- dataset_names[i]
        name2 <- dataset_names[j]
        
        cat("\nComparing", basename(name1), "with", basename(name2), "\n")
        
        comparison <- all.equal(datasets[[name1]], 
                                datasets[[name2]])
        
        comparison_results <- rbind(comparison_results, data.frame(
          Dataset1 = name1,
          Dataset2 = name2,
          Identical = isTRUE(comparison),
          Differences = ifelse(isTRUE(comparison), 
                               "None", 
                               paste(comparison, collapse = "; ")),
          stringsAsFactors = FALSE
        ))
      }, error = function(e) {
        warning("❌ Error comparing ", name1, " and ", name2, ": ", e$message)
      })
    }
  }
  
  # Print dataset summaries
  cat("\n=== Dataset Summaries ===\n")
  for (name in names(datasets)) {
    cat("\n", basename(name), ":\n")
    cat("Rows:", nrow(datasets[[name]]), "\n")
    cat("Columns:", ncol(datasets[[name]]), "\n")
    cat("Column names:", paste(colnames(datasets[[name]]), collapse = ", "), "\n")
  }
  
  return(list(
    comparison_results = comparison_results,
    datasets = datasets
  ))
}

# Function to perform verification for specific comparison
perform_verification <- function(description, path1, path2) {
  cat("\n=== Verifying", description, "===\n")
  
  file_paths <- c(path1, path2)
  
  results <- tryCatch({
    verify_multiple_datasets(file_paths)
  }, error = function(e) {
    cat("❌ Error in verification:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(results)) {
    cat("\n=== Results for", description, "===\n")
    print(results$comparison_results)
    
    cat("\n=== Sample Data for", description, "===\n")
    for (name in names(results$datasets)) {
      cat("\nFirst few rows of", basename(name), ":\n")
      print(head(results$datasets[[name]]))
    }
  }
  
  return(results)
}

# Base paths
base_path_1 <- ""
base_path_2 <- ""
base_path_3 <- ""
base_path_4 <- ""



# Perform verifications
# 1. Original Data
original_results <- perform_verification(
  "Original Data",
  file.path(),
  file.path()
)

# 2. Sampled Data
sampled_results_20 <- perform_verification(
  "Sampled Data",
  file.path(),
  file.path()
)

sampled_results_30 <- perform_verification(
  "Sampled Data",
  file.path(),
  file.path()
)

sampled_results_40 <- perform_verification(
  "Sampled Data",
  file.path(),
  file.path()
)
# 3. Resampled Data
resampled_results_20 <- perform_verification(
  "Resampled Data",
  file.path(),
  file.path()
)

resampled_results_30 <- perform_verification(
  "Resampled Data",
  file.path(),
  file.path()
)

resampled_results_40 <- perform_verification(
  "Resampled Data",
  file.path(),
  file.path()
)

# 4. Summary Data
summary_results_20 <- perform_verification(
  "Summary Data",
  file.path(),
  file.path()
)

summary_results_30 <- perform_verification(
  "Summary Data",
  file.path(),
  file.path()
)

summary_results_40 <- perform_verification(
  "Summary Data",
  file.path(),
  file.path()
)


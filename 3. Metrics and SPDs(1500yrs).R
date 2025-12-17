library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(zoo)
library(patchwork)
library(tidyr)
library(randomcoloR)
library(ggsci)

########################STEP 1: Merging all calibrated data (only normalized Probability density) ########################################
parse_filename <- function(filename) {
  # Extract seed numbers in order from highest to lowest to avoid overlapping
  seed_pattern <- paste0(
    "S50_RS1|",                    # 50 first
    "S4[0-9]_RS1|",               # 40-49
    "S3[0-9]_RS1|",               # 30-39
    "S2[0-9]_RS1|",               # 20-29
    "S1[0-9]_RS1|",               # 10-19
    "S[1-9]_RS1"                  # 1-9 last
  )
  
  seed_match <- str_extract(filename, seed_pattern)
  seed_num <- str_extract(seed_match, "S\\d+")
  
  # Extract distributions
  settlement_dist <- ifelse(grepl("power_law", filename), "power_law", "normal")
  
  bp_dist <- case_when(
    grepl("_normal.csv$", filename) ~ "normal",
    grepl("_skewed.csv$", filename) ~ "skewed",
    grepl("_uniform.csv$", filename) ~ "uniform",
    TRUE ~ "unknown"
  )
  
  # Add error checking
  if(is.na(seed_num)) {
    warning(paste("Could not extract seed number from filename:", filename))
    seed_num <- "Unknown"
  }
  
  result <- list(
    seed = seed_num,
    settlement_dist = settlement_dist,
    bp_dist = bp_dist
  )
  
  return(result)
}

# Test the function with various cases
test_filenames <- c(
  "all_plot_df_S1_RS1_normal_normal.csv",
  "all_plot_df_S9_RS1_normal_normal.csv",
  "all_plot_df_S10_RS1_normal_normal.csv",
  "all_plot_df_S19_RS1_normal_normal.csv",
  "all_plot_df_S20_RS1_normal_normal.csv",
  "all_plot_df_S29_RS1_normal_normal.csv",
  "all_plot_df_S30_RS1_normal_normal.csv",
  "all_plot_df_S39_RS1_normal_normal.csv",
  "all_plot_df_S40_RS1_normal_normal.csv",
  "all_plot_df_S49_RS1_normal_normal.csv",
  "all_plot_df_S50_RS1_normal_normal.csv"
)

# Run tests
print("Testing seed extraction:")
for(test_file in test_filenames) {
  result <- parse_filename(test_file)
  print(paste(test_file, "->", result$seed))
}

# Define file paths for each sampling percentage
file_paths <- list(
  "20" = "",
  "30" = "",
  "40" =""
)

# Read normalized data only
normalized_data <- list()

for(pct in names(file_paths)) {
  print(paste("\nProcessing sampling percentage:", pct, "%"))
  path <- file_paths[[pct]]
  
  print(paste("Reading from directory:", path))
  
  files <- list.files(path = path,
                      pattern = "normalized_.*csv$",
                      full.names = TRUE)
  
  print(paste("Number of files found:", length(files)))
  
  data_list <- map_df(files, function(file) {
    tryCatch({
      print(paste("\nProcessing file:", basename(file)))
      
      data <- read.csv(file)
      print(paste("  Rows in raw data:", nrow(data)))
      
      filename <- basename(file)
      info <- parse_filename(filename)
      print(paste("  Extracted seed:", info$seed))
      print(paste("  Settlement distribution:", info$settlement_dist))
      print(paste("  BP distribution:", info$bp_dist))
      
      processed_data <- data %>%
        mutate(
          Settlement_Dist = info$settlement_dist,
          BP_Dist = info$bp_dist,
          Sample_Seed = as.numeric(str_extract(info$seed, "\\d+")),
          Sampling_Pct = as.numeric(pct)
        )
      
      print(paste("  Rows in processed data:", nrow(processed_data)))
      print(paste("  Unique Sample Seeds:", paste(unique(processed_data$Sample_Seed), collapse = ", ")))
      
      return(processed_data)
      
    }, error = function(e) {
      print(paste("ERROR processing file:", basename(file)))
      print(paste("Error message:", e$message))
      return(NULL)
    })
  })
  
  print(paste("\nSummary for", pct, "% sampling:"))
  print(paste("Total rows in combined data:", nrow(data_list)))
  print("Sample Seed distribution:")
  print(table(data_list$Sample_Seed))
  
  normalized_data[[pct]] <- data_list
}

# Combine all normalized data
print("\nCombining all data...")
final_data_1 <- bind_rows(normalized_data)

# Final verification steps
print("\nFinal Data Summary:")
print("Total number of rows:")
print(nrow(final_data_1))

print("\nSample_Seed distribution:")
print(table(final_data_1$Sample_Seed))

print("\nCombinations of Sample_Seed and Sampling_Pct:")
print(table(final_data_1$Sample_Seed, final_data_1$Sampling_Pct))

print("\nUnique combinations of Settlement_Dist and BP_Dist:")
print(table(final_data_1$Settlement_Dist, final_data_1$BP_Dist))

# Save the combined data
print("\nSaving combined data...")
write.csv(final_data_1, "combined_normalized_samples_1to50_1500yrs.csv", row.names = FALSE)
print("Data saved successfully!")

##In case when the final data is already saved, load it in
final_data_1<-read.csv("")


########################STEP 2-1: RMSE calculation by settlement########################################
# Function to calculate RMSE norm
calculate_RMSE_norm_per_year <- function(original, comparison) {
  if (length(original) != length(comparison)) {
    stop("Vectors must be of the same length")
  }
  n <- length(original)  
  sqrt(sum((original - comparison)^2) / n)
}

### Calculate RMSE by settlement with rollingmean
## in this case rolling mean by 50 years

calculate_metrics_by_settlement_rm <- function(data, k = 50, keep_area = TRUE) {
  print("\nStarting metrics calculation...")
  
  result <- data %>% 
    {
      print("\nAfter filtering for Normalized data:")
      print(paste("Number of rows:", nrow(.)))
      print("\nUnique combinations:")
      print("Settlement_Dist:")
      print(unique(.$Settlement_Dist))
      print("BP_Dist:")
      print(unique(.$BP_Dist))
      print("Sampling_Pct:")
      print(unique(.$Sampling_Pct))
      print("Sample_Seeds:")
      print(unique(.$Sample_Seed))
      .
    } %>%
    #──────────────────────── ① Aggregate Probability Density by settlement ────────────────────────
    group_by(calBP, Dataset, Sample_Seed, Settlement,
             Settlement_Dist, BP_Dist, Sampling_Pct) %>% 
    summarise(Aggregated_PrDens = sum(Adjusted_PrDens), .groups = "drop") %>% 
    {
      print("\nAfter initial aggregation:")
      print("Number of unique combinations:")
      print(paste("Sample Seeds:", n_distinct(.$Sample_Seed)))
      print(paste("Settlements:", n_distinct(.$Settlement)))
      print(paste("Settlement_Dist:", n_distinct(.$Settlement_Dist)))
      print(paste("BP_Dist:", n_distinct(.$BP_Dist)))
      print(paste("Sampling_Pct:", n_distinct(.$Sampling_Pct)))
      .
    } %>%
    #──────────────────────── ② rolling mean  ───────────────────
    group_by(Dataset, Sample_Seed, Settlement,
             Settlement_Dist, BP_Dist, Sampling_Pct) %>% 
    arrange(calBP, .by_group = TRUE) %>% 
    mutate(
      smoothed = zoo::rollmean(Aggregated_PrDens, k = k,
                               fill = NA, align = "center"),
      smoothed = if (keep_area)
        smoothed / sum(smoothed, na.rm = TRUE) *
        sum(Aggregated_PrDens, na.rm = TRUE)
      else smoothed
    ) %>% 
    ungroup() %>% 
    
    #──────────────────────── ③ RMSE ───────────────────────────
    group_by(Sample_Seed, Settlement,
             Settlement_Dist, BP_Dist, Sampling_Pct) %>% 
    summarise(
      Original_Values  = list(smoothed[Dataset == "Original"]),
      Sampled_Values   = list(smoothed[Dataset == "Sampled"]),
      Weighted_Values  = list(smoothed[Dataset == "Weighted Sampled"]),
      Resampled30_Values = list(smoothed[Dataset == "Resampled30"]),
      .groups = "drop"
    ) %>% 
    mutate(
      RMSE_Sampled  = map2_dbl(Original_Values, Sampled_Values,  ~sqrt(mean((.x - .y)^2, na.rm = TRUE))),
      RMSE_Weighted = map2_dbl(Original_Values, Weighted_Values, ~sqrt(mean((.x - .y)^2, na.rm = TRUE))),
      RMSE_Resampled30 = map2_dbl(Original_Values, Resampled30_Values, ~sqrt(mean((.x - .y)^2, na.rm = TRUE)))
    ) %>% 
    select(-ends_with("_Values"))
  
  return(result)
}


# Process each sampling percentage separately
sampling_pcts <- unique(final_data_1$Sampling_Pct)
result_list <- list()
for(pct in sampling_pcts) {
  print(paste("\nProcessing sampling percentage:", pct))
  
  pct_data <- final_data_1 %>% filter(Sampling_Pct == pct)
  result_list[[as.character(pct)]] <- calculate_metrics_by_settlement_rm(pct_data)
  
  # Clean up
  rm(pct_data)
  gc()
}

# Combine results
metrics_by_settlement_1 <- bind_rows(result_list)

# Save raw results (all KL and RMSE values for each settlement/combination)
write.csv(metrics_by_settlement_1, "raw_metrics_by_settlement_1500yrs_RM50.csv", row.names = FALSE)


# Check if SPDs sum to 1
check_spd_sums <- function(data) {
  spd_sums <- data %>%
    group_by(Dataset, Sample_Seed, Settlement_Dist, BP_Dist, Sampling_Pct) %>%
    summarise(
      SPD_Sum = sum(Adjusted_PrDens),
      .groups = "drop"
    )
  
  # Print any SPDs that do not sum to 1
  incorrect_sums <- spd_sums %>%
    filter(abs(SPD_Sum - 1) > 1e-10)  # Allowing for a small numerical tolerance
  
  if (nrow(incorrect_sums) > 0) {
    print("SPDs that do not sum to 1:")
    print(incorrect_sums)
  } else {
    print("All SPDs sum to 1.")
  }
}

# Call the function with your data
check_spd_sums(final_data_1)

# Check the number of rows for each sample seed
metrics_by_settlement_1 %>%
  group_by(Sample_Seed) %>%
  summarise(n = n()) %>%
  print(n = Inf)

# Check the complete structure
metrics_by_settlement_1 %>%
  group_by(Sample_Seed, Settlement,Settlement_Dist,BP_Dist) %>%
  summarise(n = n()) %>%
  print(n = Inf)


metrics_by_settlement_1 <- metrics_by_settlement_1 %>%
  rename(
    RMSE_WeightedSampled = RMSE_Weighted,
    RMSE_Resampled30 = RMSE_Resampled30
  )



# Function to perform both Wilcoxon and t-tests
perform_statistical_tests <- function(data) {
  # Function to calculate rank details for Wilcoxon test
  calculate_rank_details <- function(x, y) {
    differences <- x - y
    abs_diff <- abs(differences)
    ranks <- rank(abs_diff)
    
    # Split ranks by sign of difference
    positive_ranks <- sum(ranks[differences > 0])
    negative_ranks <- sum(ranks[differences < 0])
    
    # Count positive and negative differences
    n_positive <- sum(differences > 0)
    n_negative <- sum(differences < 0)
    n_ties <- sum(differences == 0)
    
    return(list(
      positive_ranks_sum = positive_ranks,
      negative_ranks_sum = negative_ranks,
      n_positive = n_positive,
      n_negative = n_negative,
      n_ties = n_ties
    ))
  }
  
  # ============== WILCOXON TESTS ==============
  # Perform Wilcoxon tests and calculate ranks
  sampled_vs_weighted_wilcox <- wilcox.test(
    data$RMSE_Sampled,
    data$RMSE_WeightedSampled,
    paired = TRUE
  )
  weighted_ranks <- calculate_rank_details(data$RMSE_Sampled, data$RMSE_WeightedSampled)
  
  sampled_vs_resampled_wilcox <- wilcox.test(
    data$RMSE_Sampled,
    data$RMSE_Resampled30,
    paired = TRUE
  )
  resampled_ranks <- calculate_rank_details(data$RMSE_Sampled, data$RMSE_Resampled30)
  
  # Create Wilcoxon results data frame
  wilcox_results <- data.frame(
    Comparison = c("Sampled vs WeightedSampled", "Sampled vs Resampled30"),
    Test = "Wilcoxon",
    p_value = c(sampled_vs_weighted_wilcox$p.value, sampled_vs_resampled_wilcox$p.value),
    statistic = c(sampled_vs_weighted_wilcox$statistic, sampled_vs_resampled_wilcox$statistic),
    mean_first = c(mean(data$RMSE_Sampled), mean(data$RMSE_Sampled)),
    mean_second = c(mean(data$RMSE_WeightedSampled), mean(data$RMSE_Resampled30)),
    sd_first = c(sd(data$RMSE_Sampled), sd(data$RMSE_Sampled)),
    sd_second = c(sd(data$RMSE_WeightedSampled), sd(data$RMSE_Resampled30)),
    positive_ranks_sum = c(weighted_ranks$positive_ranks_sum, resampled_ranks$positive_ranks_sum),
    negative_ranks_sum = c(weighted_ranks$negative_ranks_sum, resampled_ranks$negative_ranks_sum),
    n_positive = c(weighted_ranks$n_positive, resampled_ranks$n_positive),
    n_negative = c(weighted_ranks$n_negative, resampled_ranks$n_negative),
    n_ties = c(weighted_ranks$n_ties, resampled_ranks$n_ties)
  ) %>%
    mutate(
      significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      difference = mean_second - mean_first
    )
  
  # ============== T-TESTS ==============
  # Perform paired t-tests
  sampled_vs_weighted_t <- t.test(
    data$RMSE_Sampled,
    data$RMSE_WeightedSampled,
    paired = TRUE
  )
  
  sampled_vs_resampled_t <- t.test(
    data$RMSE_Sampled,
    data$RMSE_Resampled30,
    paired = TRUE
  )
  
  # Create t-test results data frame
  t_results <- data.frame(
    Comparison = c("Sampled vs WeightedSampled", "Sampled vs Resampled30"),
    Test = "t-test",
    p_value = c(sampled_vs_weighted_t$p.value, sampled_vs_resampled_t$p.value),
    statistic = c(sampled_vs_weighted_t$statistic, sampled_vs_resampled_t$statistic),
    df = c(sampled_vs_weighted_t$parameter, sampled_vs_resampled_t$parameter),
    conf_int_lower = c(sampled_vs_weighted_t$conf.int[1], sampled_vs_resampled_t$conf.int[1]),
    conf_int_upper = c(sampled_vs_weighted_t$conf.int[2], sampled_vs_resampled_t$conf.int[2]),
    mean_first = c(mean(data$RMSE_Sampled), mean(data$RMSE_Sampled)),
    mean_second = c(mean(data$RMSE_WeightedSampled), mean(data$RMSE_Resampled30)),
    sd_first = c(sd(data$RMSE_Sampled), sd(data$RMSE_Sampled)),
    sd_second = c(sd(data$RMSE_WeightedSampled), sd(data$RMSE_Resampled30))
  ) %>%
    mutate(
      significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      difference = mean_second - mean_first
    )
  
  # Print detailed results
  cat("\nWILCOXON SIGNED-RANK TEST RESULTS:\n")
  for(i in 1:nrow(wilcox_results)) {
    cat(sprintf("\n%s:\n", wilcox_results$Comparison[i]))
    cat(sprintf("p-value: %.4f %s\n", 
                wilcox_results$p_value[i],
                wilcox_results$significance[i]))
    cat(sprintf("Mean difference: %.4f\n", 
                wilcox_results$difference[i]))
    cat(sprintf("First method mean (SD): %.4f (%.4f)\n", 
                wilcox_results$mean_first[i],
                wilcox_results$sd_first[i]))
    cat(sprintf("Second method mean (SD): %.4f (%.4f)\n", 
                wilcox_results$mean_second[i],
                wilcox_results$sd_second[i]))
    cat("\nRank Details:\n")
    cat(sprintf("Positive ranks sum: %.1f (n = %d)\n",
                wilcox_results$positive_ranks_sum[i],
                wilcox_results$n_positive[i]))
    cat(sprintf("Negative ranks sum: %.1f (n = %d)\n",
                wilcox_results$negative_ranks_sum[i],
                wilcox_results$n_negative[i]))
    cat(sprintf("Number of ties: %d\n",
                wilcox_results$n_ties[i]))
  }
  
  cat("\nPAIRED T-TEST RESULTS:\n")
  for(i in 1:nrow(t_results)) {
    cat(sprintf("\n%s:\n", t_results$Comparison[i]))
    cat(sprintf("p-value: %.4f %s\n", 
                t_results$p_value[i],
                t_results$significance[i]))
    cat(sprintf("t-statistic: %.4f (df = %.1f)\n", 
                t_results$statistic[i],
                t_results$df[i]))
    cat(sprintf("Mean difference: %.4f\n", 
                t_results$difference[i]))
    cat(sprintf("95%% Confidence Interval: [%.4f, %.4f]\n",
                t_results$conf_int_lower[i],
                t_results$conf_int_upper[i]))
    cat(sprintf("First method mean (SD): %.4f (%.4f)\n", 
                t_results$mean_first[i],
                t_results$sd_first[i]))
    cat(sprintf("Second method mean (SD): %.4f (%.4f)\n", 
                t_results$mean_second[i],
                t_results$sd_second[i]))
  }
  
  return(list(wilcoxon = wilcox_results, t_test = t_results))
}

# Usage:
statistical_results <- perform_statistical_tests(metrics_by_settlement_1)


# Save results
write.csv(statistical_results, "wilcoxon_ttest_results_by_settlement_RM50(1500yrs).csv", row.names = FALSE)


########################STEP 2-2: Plots for RMSE by settlement ########################################

metrics_by_settlement_1<-read.csv("")

# Create long format for RMSE metric
metrics_long_1 <- metrics_by_settlement_1 %>%
  pivot_longer(
    cols = starts_with("RMSE_"),
    names_to = "Dataset_Type",
    names_pattern = "RMSE_(.*)",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = "RMSE",  
    Dataset_Type = factor(Dataset_Type,
                          levels = c("Sampled", "Weighted", "Resampled30"),
                          labels = c("Random Sample", "Weighted", "Bootstrapped"))
  )

# Get unique combinations of Settlement_Dist, BP_Dist, and Sampling_Pct
dist_combinations <- metrics_by_settlement_1 %>%
  select(Settlement_Dist, BP_Dist, Sampling_Pct) %>%
  distinct()

# Create plots for each combination
for(i in 1:nrow(dist_combinations)) {
  settlement_dist <- dist_combinations$Settlement_Dist[i]
  bp_dist <- dist_combinations$BP_Dist[i]
  samp_pct <- dist_combinations$Sampling_Pct[i]
  
  # RMSE norm plot
  p_RMSE <- metrics_long_1 %>%
    filter(Settlement_Dist == settlement_dist,
           BP_Dist == bp_dist,
           Sampling_Pct == samp_pct,
           Metric == "RMSE") %>%
    ggplot(aes(x = factor(Settlement), y = Value, fill = Dataset_Type)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    theme_bw() +
    labs(title = paste("Boxplot of RMSE by Settlement\n",
                       settlement_dist, "-", bp_dist, "at", samp_pct, "%"),
         x = "Settlement",
         y = "RMSE",
         fill = "Dataset Type") +
    theme(
      axis.text.x = element_text(size = 18, angle = 0,face = "bold"),
      axis.text.y = element_text(size = 18, angle = 0, face = "bold"),
      axis.title.x = element_text(size = 15,margin = margin(t = 15, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(size = 15),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      strip.text = element_text(size = 15),
      strip.background = element_rect(fill = "lightgray"),
      legend.position = "bottom",
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 12),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA)
    ) +
    scale_fill_bmj() +
    scale_x_discrete(limits = as.character(1:10))
  
  
  ggsave(paste0("rmse_norm_boxplots_", 
                settlement_dist, "_", 
                bp_dist, "_",
                samp_pct, "percent.png"),
         p_RMSE, width = 10, height = 6)
}

# Create summary statistics
summary_stats <- metrics_long_1 %>%
  group_by(Settlement_Dist, BP_Dist, Sampling_Pct, Settlement, Metric, Dataset_Type) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    Median = median(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    Q1 = quantile(Value, 0.25, na.rm = TRUE),
    Q3 = quantile(Value, 0.75, na.rm = TRUE),
    Min = min(Value, na.rm = TRUE),
    Max = max(Value, na.rm = TRUE),
    N = n()  # Number of seeds
  ) %>%
  ungroup() %>%
  arrange(Metric, Settlement_Dist, BP_Dist, Sampling_Pct, Settlement, Dataset_Type)

# Save summary statistics
write.csv(summary_stats, "metrics_by_settlement_summary_1500yrs.csv", row.names = FALSE)


#########################STEP 2-3 (Figure 7b) : Compute RMSE by settlement and Boxplot visualization with facet_wrap#########################
metrics_long_1 <- metrics_long_1 %>%
  mutate(
    Settlement_Dist = ifelse(Settlement_Dist == "power_law", "power law", Settlement_Dist),
    Combined_Dist = paste(Settlement_Dist, BP_Dist, sep = "-")
  )


# RMSE norm plot
p_RMSE <- metrics_long_1 %>%
  filter(Metric =="RMSE") %>%  # Adjust percentage as needed
  ggplot(aes(x = factor(Settlement), y = Value, fill = Dataset_Type)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  facet_wrap(~ Combined_Dist, ncol = 3) +  # Facet by Combined_Dist
  theme_bw() +
  labs(title = "RMSE by Settlement (Maximum Sample Fraction 30%)",
       x = "Settlement",
       y = "RMSE",
       fill = "Dataset Type") +
  theme(
    axis.text.x = element_text(size = 30, angle = 0, face = "bold"),
    axis.text.y = element_text(size = 30, angle = 0, face = "bold"),
    axis.title.x = element_text(size = 30,  face = "bold", margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 30, face="bold"),
    plot.title = element_text(size = 35, face = "bold", hjust=0.5),
    strip.text = element_text(size = 40, face = "bold"),
    strip.background = element_rect(fill = "lightgray"),
    legend.position = "bottom",
    legend.title = element_text(size = 30,  face = "bold"),
    legend.text = element_text(size = 30, face="bold"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA)
  ) +
  scale_fill_bmj() +
  scale_x_discrete(limits = as.character(1:10))

# Save plots
ggsave("rmse_boxplots_facet_1500yrs_RM50.png", p_RMSE, width = 20, height = 12)


########################STEP 3-1: Combine SPDs by settlement to create an overall SPD, then compute RMSE########################################

#rolling mean, in this case rolling mean by 50 years
calculate_combined_metrics_rm <- function(data, k = 50, keep_area = TRUE) {
  result <- data %>%
    group_by(calBP, Dataset, Sample_Seed, Settlement_Dist, BP_Dist, Sampling_Pct) %>%
    summarise(Aggregated_PrDens = sum(Adjusted_PrDens), .groups = "drop") %>%
    group_by(Dataset, Sample_Seed, Settlement_Dist, BP_Dist, Sampling_Pct) %>%
    arrange(calBP, .by_group = TRUE) %>%
    mutate(
      smoothed = zoo::rollmean(Aggregated_PrDens, k = k, fill = NA, align = "center"),
      smoothed = if (keep_area)
        smoothed / sum(smoothed, na.rm = TRUE) * sum(Aggregated_PrDens, na.rm = TRUE)
      else smoothed
    ) %>%
    ungroup() %>%
    group_by(Sample_Seed, Settlement_Dist, BP_Dist, Sampling_Pct) %>%
    summarise(
      Original_Values = list(smoothed[Dataset == "Original"]),
      Sampled_Values = list(smoothed[Dataset == "Sampled"]),
      WeightedSampled_Values = list(smoothed[Dataset == "Weighted Sampled"]),
      Resampled30_Values = list(smoothed[Dataset == "Resampled30"]),
      .groups = "drop"
    ) %>%
    mutate(
      RMSE_Sampled = map2_dbl(Original_Values, Sampled_Values, ~sqrt(mean((.x - .y)^2, na.rm = TRUE))),
      RMSE_WeightedSampled = map2_dbl(Original_Values, WeightedSampled_Values, ~sqrt(mean((.x - .y)^2, na.rm = TRUE))),
      RMSE_Resampled30 = map2_dbl(Original_Values, Resampled30_Values, ~sqrt(mean((.x - .y)^2, na.rm = TRUE)))
    ) %>%
    select(-ends_with("_Values"))
  
  return(result)
}


# Calculate metrics
metrics_combined_1 <- calculate_combined_metrics_rm(final_data_1)

# Save raw results
write.csv(metrics_combined_1, "raw_metrics_combined__settlements_1500yrs_overallSPD_RM50.csv", row.names = FALSE)

# Create summary statistics for both KL and RMSE metrics
summary_stats <- metrics_combined_1 %>%
  group_by(Settlement_Dist, BP_Dist, Sampling_Pct) %>%
  summarise(
    # RMSE Sampled
    RMSE_Sampled_Mean = mean(RMSE_Sampled),
    RMSE_Sampled_SD = sd(RMSE_Sampled),
    RMSE_Sampled_Median = median(RMSE_Sampled),
    RMSE_Sampled_Q1 = quantile(RMSE_Sampled, 0.25),
    RMSE_Sampled_Q3 = quantile(RMSE_Sampled, 0.75),
    
    # RMSE WeightedSampled
    RMSE_WeightedSampled_Mean = mean(RMSE_WeightedSampled),
    RMSE_WeightedSampled_SD = sd(RMSE_WeightedSampled),
    RMSE_WeightedSampled_Median = median(RMSE_WeightedSampled),
    RMSE_WeightedSampled_Q1 = quantile(RMSE_WeightedSampled, 0.25),
    RMSE_WeightedSampled_Q3 = quantile(RMSE_WeightedSampled, 0.75),
    
    # RMSE Resampled30
    RMSE_Resampled30_Mean = mean(RMSE_Resampled30),
    RMSE_Resampled30_SD = sd(RMSE_Resampled30),
    RMSE_Resampled30_Median = median(RMSE_Resampled30),
    RMSE_Resampled30_Q1 = quantile(RMSE_Resampled30, 0.25),
    RMSE_Resampled30_Q3 = quantile(RMSE_Resampled30, 0.75),
    
    # Number of samples
    n_samples = n()
  ) %>%
  ungroup()

# Save summary statistics
write.csv(summary_stats, "summary_metrics_combined_settlements_1500yrs_RM50.csv", row.names = FALSE)



########################STEP 3-2: Merge settlement SPDs into an overall SPD and plot RMSE as boxplots########################################
metrics_long_combined_1 <- metrics_combined_1 %>%
  pivot_longer(
    cols = starts_with("RMSE_"),
    names_to = "Dataset_Type",
    names_pattern = "RMSE_(.*)",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = "RMSE",  # Add this column since we're filtering by it in the plot
    Dataset_Type = factor(Dataset_Type,
                          levels = c("Sampled", "WeightedSampled", "Resampled30"),
                          labels = c("Random Sample", "Weighted", "Bootstrapped"))
  )


# Get unique combinations
dist_combinations <- metrics_combined_1 %>%
  select(Settlement_Dist, BP_Dist) %>%
  distinct()

# Create plots for each combination
for(i in 1:nrow(dist_combinations)) {
  settlement_dist <- dist_combinations$Settlement_Dist[i]
  bp_dist <- dist_combinations$BP_Dist[i]
  
  # Set square dimensions
  plot_size <- 8  # This will be used for both width and height
  
  # RMSE norm plot
  p_RMSE <- metrics_long_combined_1 %>%
    filter(Settlement_Dist == settlement_dist,
           BP_Dist == bp_dist) %>%
    ggplot(aes(x = factor(Sampling_Pct), y = Value, fill = Dataset_Type)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    theme_bw() +
    labs(title = paste("Boxplot of RMSE by Sampling Percentage\n",
                       settlement_dist, "-", bp_dist),
         x = "Sampling Percentage",
         y = "RMSE",
         fill = "Dataset Type") +
    theme(
      axis.text.x = element_text(size = 18, angle = 0, face = "bold"),
      axis.text.y = element_text(size = 18, angle = 0, face = "bold"),
      axis.title.x = element_text(size = 15,
                                  margin = margin(t = 15, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(size = 15),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      legend.position = "bottom",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      panel.grid.minor = element_blank(),
      # Add extra margin to accommodate the legend in square format
      plot.margin = margin(t = 20, r = 20, b = 40, l = 20)
    ) +
    scale_fill_bmj() +
    scale_x_discrete(limits = c("20", "30", "40"))
  
  ggsave(paste0("RMSE_norm_combined_", 
                settlement_dist, "_", 
                bp_dist, ".png"),
         p_RMSE, 
         width = plot_size, 
         height = plot_size)
}

######################## STEP 3-2 (Figure 4b): Combine settlement SPDs into overall SPD and draw RMSE boxplots with facet_wrap ########################################
create_faceted_comparison_plots <- function(metrics_long_combined_1) {
  # RMSE norm plot with facets
  p_RMSE_faceted <- metrics_long_combined_1 %>%
    filter(Metric == "RMSE") %>%
    ggplot(aes(x = factor(Sampling_Pct), y = Value, fill = Dataset_Type)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    facet_grid(Settlement_Dist ~ BP_Dist) +
    theme_bw() +
    labs(title = "Boxplot of RMSE by Maximum Sample Fraction",
         x = "Maximum Sample Fractions",
         y = "RMSE",
         fill = "Dataset Type") +
    theme(
      axis.text.x = element_text(size = 25, angle = 0, face = "bold"),
      axis.text.y = element_text(size = 25, angle = 0, face = "bold"),
      axis.title.x = element_text(size = 23,face = "bold",
                                  margin = margin(t = 15, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(size = 25,face = "bold", margin = margin(t = 0, r = 15, b = 0, l = 0)),
      plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_text(size = 25, face = "bold"),    # Increase from 15
      legend.text = element_text(size = 25),                    # Increase from 12
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 25, face = "bold"),
      strip.background = element_rect(fill = "lightgray"),
      plot.margin = margin(t = 20, r = 20, b = 40, l = 20)
    ) +
    scale_fill_bmj() +
    scale_x_discrete(limits = c("20", "30", "40"))
  
  plot_width <- 12
  plot_height <- 8
  
  ggsave("RMSE_norm_faceted_all_1500yrs.png", 
         p_RMSE_faceted, 
         width = plot_width, 
         height = plot_height)
  
  return(list(RMSE = p_RMSE_faceted))
}

faceted_plots <- create_faceted_comparison_plots(metrics_long_combined_1)


######################## STEP 3-3 (Figure 1b): Combine settlement SPDs into overall SPD and draw RMSE boxplots (original vs. sample only) #######################################

rm_colors <- pal_rickandmorty()(12)
modified_colors <- c(
  "20" = rm_colors[2],  # Blue
  "30" = rm_colors[3],  # Purple 
  "40" = rm_colors[7]   # Green
)

p_RMSE <- metrics_long_combined_1 %>%
  filter(Dataset_Type == "Random Sample",
         Metric == "RMSE") %>%
  ggplot(aes(x = factor(Sampling_Pct), y = Value, fill = factor(Sampling_Pct))) +
  geom_boxplot(alpha = 0.35) +
  geom_jitter(width = 0.1, aes(color = factor(Sampling_Pct))) + 
  facet_grid(Settlement_Dist ~ BP_Dist) +
  theme_bw() +
  labs(title = "Hypothetical Population vs. Random Sample Set",
       x = "Maximum Sample Fractions",
       y = "RMSE",
       fill = "Maximum Sample Fractions %",
       color ="Maximum Sample Fractions %") +
  theme(
    axis.text.x = element_text(size = 30, angle = 0, face = "bold"),
    axis.text.y = element_text(size = 30, angle = 0, face = "bold"),
    axis.title.x = element_text(size = 30, face = "bold",
                                margin = margin(t = 20, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 30, face = "bold",
                                margin = margin(t = 0, r = 20, b = 0, l = 0)),
    plot.title = element_text(size = 33, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 35, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size =30, face = "bold"),
    legend.text = element_text(size = 30),
    legend.key.size = unit(1.5, "cm"),
    legend.margin = margin(t = 10, r = 0, b = 0, l = 0),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 30, r = 30, b = 30, l = 30)
  ) + 
  scale_x_discrete(limits = c("20", "30", "40")) +
  scale_fill_manual(values = modified_colors) +
  scale_color_manual(values = modified_colors)

width <- 15
height <- width * 0.9

ggsave("RMSE_norm_sampled_all_combos_jitter_1500yrs_RM50.png", p_RMSE, width = width, height = height)



######################### STEP 4-1: Merge settlement SPDs into an overall SPD and visualize it for original, resampled, weighted, and sample datasets########################################

aggregate_spd_data <- function(data) {
  aggregated_data <- data %>%
    group_by(calBP, Dataset, Settlement_Dist, BP_Dist, Sample_Seed, Sampling_Pct) %>%
    summarise(Adjusted_PrDens = sum(Adjusted_PrDens), .groups = "drop")
  
  return(aggregated_data)
}


agg_data_1 <- aggregate_spd_data(final_data_1)
head(agg_data_1)
unique(agg_data_1$Dataset)


##### Plot separately for maximum sampling fractions: 20%, 30%, and 40% #####

create_spd_comparison_plots_1500 <- function(agg_data,
                                        smooth_k    = 50,
                                        percentages = c(20, 30, 40)) {
  combinations <- agg_data %>%
    select(Settlement_Dist, BP_Dist, Sample_Seed) %>%
    distinct()
  
  dir.create("SPD_plots", showWarnings = FALSE)
  
  for (i in seq_len(nrow(combinations))) {
    combo <- combinations[i, ]
    
    for (pct in percentages) {
      plot_data <- agg_data %>%
        filter(Settlement_Dist == combo$Settlement_Dist,
               BP_Dist == combo$BP_Dist,
               ((Dataset %in% c("Original", "Sampled", "Weighted Sampled", "Resampled30") &
                   Sample_Seed == combo$Sample_Seed)),
               Sampling_Pct == pct) %>%
        group_by(calBP, Dataset) %>%
        summarise(plot_density = sum(Adjusted_PrDens), .groups = "drop") %>%
        group_by(Dataset) %>%
        arrange(calBP, .by_group = TRUE) %>%
        mutate(
          plot_density_smooth =
            zoo::rollmean(plot_density, k = smooth_k, fill = 0, align = "center"),
          plot_density_smooth =
            plot_density_smooth / sum(plot_density_smooth) * sum(plot_density)
        ) %>%
        ungroup()
      
      
      if (nrow(plot_data) == 0) {
        message("Skipped (no data): ",
                combo$Settlement_Dist, "_",
                combo$BP_Dist, "_seed", combo$Sample_Seed, "_", pct, "pct")
        next
      }
      
      original_data <- plot_data %>% filter(Dataset == "Original")
      other_data    <- plot_data %>% filter(Dataset != "Original") %>%
        mutate(Dataset = factor(Dataset,
                                levels = c("Sampled", "Resampled30", "Weighted Sampled")))
      weighted_data    <- other_data %>% filter(Dataset == "Weighted Sampled")
      other_lines_data <- other_data %>% filter(Dataset != "Weighted Sampled")
      
      p <- ggplot() +
        geom_area(data = original_data,
                  aes(calBP, plot_density_smooth, fill = "Original"),
                  alpha = 0.8) +
        geom_line(data = other_lines_data,
                  aes(calBP, plot_density_smooth, colour = Dataset),
                  linewidth = 2) +
        geom_line(data = weighted_data,
                  aes(calBP, plot_density_smooth, colour = Dataset),
                  linewidth = 2, linetype = "dashed") +
        scale_x_reverse(limits = c(5500, 3000),
                        breaks = seq(5500, 3000, -200)) +
        scale_color_bmj(limits = c("Sampled","Resampled30","Weighted Sampled")) +
        scale_fill_manual(values = c("Original" = "azure4")) +
        guides(fill  = guide_legend(order = 1, override.aes = list(alpha = 1)),
               colour= guide_legend(order = 2, override.aes = list(
                 linetype = c("solid","solid","dashed"),
                 linewidth = c(6,6,3)))) +
        labs(title  = sprintf("SPD Comparison (%s - %s)\nSeed %s - %d%%",
                              combo$Settlement_Dist, combo$BP_Dist,
                              combo$Sample_Seed, pct),
             x = "Years BP",
             y = "Probability Density",
             fill = NULL, colour = NULL) +
        theme_bw() +
        theme(axis.text   = element_text(size = 18, face = "bold"),
              axis.title  = element_text(size = 20, face = "bold"),
              plot.title  = element_text(size = 20, face = "bold", hjust = .5),
              legend.position = "bottom",
              legend.text = element_text(size = 15),
              legend.title= element_text(size = 15, face = "bold"),
              panel.grid.minor = element_blank())
      
      out_name <- sprintf("SPD_plots/SPD_%s_%s_seed%s_%dpct.png",
                          combo$Settlement_Dist, combo$BP_Dist,
                          combo$Sample_Seed, pct)
      ggsave(out_name, p, width = 15, height = 8, dpi = 300)
    }
  }
}


# Run the functions
create_spd_comparison_plots_1500(agg_data_1)


create_spd_comparison_plots(agg_data_1, smooth_k = 50,
                            percentages = c(30))             # when testing for specific sampling fraction



# ##### Extract and plot just 4 specific combinations (Figure 3) #####

create_spd_comparison_plots_facets <- function(final_data_1) {

  # --- Fix labels early in both datasets to match properly
  agg_data_1 <- agg_data_1 %>%
    mutate(
      Settlement_Dist = ifelse(Settlement_Dist == "power_law", "power law", Settlement_Dist)  # Fix label in main data
    )
  
  specific_combos <- tibble(
    Settlement_Dist = c("power law", "power law", "power law", "normal"),
    BP_Dist = c( "uniform","normal","skewed", "uniform" ),
    Sample_Seed = c(45, 50, 4, 30),
    Sampling_Pct = 30,
    Facet_Order = 1:4  # Set order for facets
  ) %>%
    mutate(
      Alphabet = c("(a)", "(b)", "(c)", "(d)"),  
      Combined_Label = paste0(Alphabet, " ", Settlement_Dist, " - ", BP_Dist, " (Seed: ", Sample_Seed, ")")
    )
  
  dir.create("SPD_plots", showWarnings = FALSE)
  
  plot_data <- agg_data_1 %>%
    inner_join(specific_combos,
               by = c("Settlement_Dist", "BP_Dist", "Sample_Seed", "Sampling_Pct")) %>%
    filter(Dataset %in% c("Original", "Sampled", "Weighted Sampled", "Resampled30")) %>%
    group_by(Settlement_Dist, BP_Dist, Sample_Seed, Sampling_Pct, calBP, Dataset, Combined_Label, Facet_Order) %>%
    summarise(plot_density = sum(Adjusted_PrDens), .groups = "drop") %>%
    mutate(
      Combined_Label = factor(Combined_Label, levels = specific_combos$Combined_Label)  
    )
  
  plot_data <- plot_data %>%
    group_by(Dataset, Combined_Label) %>%    
    arrange(calBP, .by_group = TRUE) %>%
    mutate(
      plot_density_smooth =
        zoo::rollmean(plot_density, k = 50, fill = NA, align = "center"),
      plot_density_smooth =
        plot_density_smooth /
        sum(plot_density_smooth, na.rm = TRUE) * sum(plot_density)
    ) %>%
    ungroup()
  
  original_data <- plot_data %>% filter(Dataset == "Original")
  other_data <- plot_data %>% filter(Dataset != "Original")
  
  # --- Set legend order
  other_data <- other_data %>%
    mutate(Dataset = recode(Dataset,
                            "Sampled" = "Random Sample",
                            "Resampled30" = "Bootstrapped",
                            "Weighted Sampled" = "Weighted"))
  
  original_data <- original_data %>%
    mutate(Dataset = recode(Dataset, "Original" = "Hypothetical Population"))
  
  weighted_data <- other_data %>% filter(Dataset == "Weighted")
  other_lines_data <- other_data %>% filter(Dataset != "Weighted")
  
  
  p <- ggplot() +
    geom_area(data = original_data,
              aes(calBP, plot_density_smooth, fill = Dataset),  
              alpha = 0.8) +
    geom_line(data = other_lines_data,
              aes(calBP, plot_density_smooth, colour = Dataset),  
              linewidth = 3) +
    geom_line(data = weighted_data,
              aes(calBP, plot_density_smooth, colour = Dataset), 
              linewidth = 3, linetype = "dashed") +
    
    facet_wrap(~ Combined_Label, ncol = 2, scales = "free_y") +

        scale_x_reverse(limits = c(5500, 3000),
                    breaks = seq(5500, 3000, -200))+

    
    scale_color_bmj(limits = c("Random Sample", "Weighted",  "Bootstrapped")) +
    scale_fill_manual(values = c("Hypothetical Population" = "azure4"))+
    guides(
      fill = guide_legend(order = 1, override.aes = list(alpha = 1)),
      color = guide_legend(order = 2, override.aes = list(
        linetype = c("solid", "solid", "dashed"),
        linewidth = c(6, 6, 3)
      ))
    )+
    labs(
      title = "SPD Comparison (Maximum Sample Fraction 30%)",
      x = "Cal BP",
      y = "Probability Density",
      fill = NULL,
      color = NULL
    ) +
    theme_bw(base_size = 34) +  
    theme(
        axis.text.x = element_text(size = 45, angle = 45, hjust = 1, face = "bold",margin = margin(t = 10)),
        axis.text.y = element_text(size = 45, face = "bold",margin = margin(r = 10)),
        axis.title.x = element_text(size = 40, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 40, face = "bold", margin = margin(r = 15)),
        plot.title = element_text(size = 45, face = "bold", hjust = 0.5),
        strip.text = element_text(size = 50, face = "bold"), 
        legend.position = "bottom",
        legend.text = element_text(size = 40, face="bold"),
        legend.title = element_text(size = 32, face = "bold"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 33, r = 33, b = 33, l = 33),
        legend.key.width = unit(3.5, "cm"),
        legend.key.height = unit(2, "cm")
        
      )
  # --- Save plot
  ggsave("SPD_Comparison_Faceted_1500_RM50.png", p,
         width = 35, height = 25, dpi = 300)
  
  return(p)
}

getwd()
# Run the function
combined_plot <- create_spd_comparison_plots_facets(agg_data_1)

# STEP 4-2 (Supplementary): Merge settlement SPDs into overall SPD and plot SPD for all 50 seeds together ######

# For Original vs Sample
create_original_vs_sampled_plots <- function(agg_data_1) {
  # Get unique combinations of Settlement_Dist and BP_Dist
  combinations <- agg_data_1 %>%
    select(Settlement_Dist, BP_Dist) %>%
    distinct()
  
  
  # For each combination
  for(i in 1:nrow(combinations)) {
    combo <- combinations[i,]
    
    # Get all data for this combination
    plot_data <- agg_data_1 %>%
      filter(Settlement_Dist == combo$Settlement_Dist,
             BP_Dist == combo$BP_Dist,
             Dataset %in% c("Original", "Sampled")) %>%
      group_by(calBP, Dataset, Sample_Seed, Sampling_Pct) %>%
      summarise(plot_density = sum(Adjusted_PrDens), .groups = "drop")
    
    # Separate original and sampled data
    original_data <- plot_data %>% 
      filter(Dataset == "Original") %>%
      distinct(calBP, plot_density, Sampling_Pct)
    
    sampled_data <- plot_data %>% 
      filter(Dataset == "Sampled")
    
    # Create plot
    p <- ggplot() +
      # Add sampled lines with continuous color scale
      geom_line(data = sampled_data,
                aes(x = calBP, y = plot_density, 
                    group = Sample_Seed,
                    color = as.numeric(Sample_Seed)),
                linewidth = 0.3,
                alpha = 0.6) +
      # Add original line
      geom_line(data = original_data,
                aes(x = calBP, y = plot_density),
                color = "tomato",
                linewidth = 2.2) +
      facet_wrap(~Sampling_Pct, 
                 labeller = labeller(Sampling_Pct = function(x) paste0(x, "% Sampling")),
                 nrow = 1) +
      scale_x_reverse(limits = c(5500, 3000),
                      breaks = seq(5500, 3000, -200))+
      scale_color_viridis_c(name = "Sampling Seed") +
      labs(title = sprintf("Original - Sampled Comparison\n%s - %s",
                           combo$Settlement_Dist,
                           combo$BP_Dist),
           x = "Years BP",
           y = "Probability Density") +
      theme_bw() +
      theme(
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18,
                                    face = "bold",
                                    margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 18,
                                    face = "bold",
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_text(size = 15, face = "bold"),  # Increase title size
        legend.text = element_text(size = 13),  # Increase text size
        legend.key.width = unit(2, "cm"),
        legend.key.height = unit(0.5, "cm"),    # Optionally adjust height
        legend.margin = margin(t = 10, b = 10),  # Add some margin
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
        strip.text = element_text(size = 25, face = "bold"),
        strip.background = element_rect(fill = "white")
      )
    
    # Save plot
    filename <- sprintf("SPD_plots/Original_vs_Sampled_%s_%s_allseeds_faceted.png",
                        combo$Settlement_Dist,
                        combo$BP_Dist)
    
    ggsave(filename, p, width = 20, height = 8, dpi = 300)
  }
}


create_original_vs_sampled_plots_wide <- function(agg_data_1) {
  distinct_colors <- c(
    "#FF4136",
    "#FF851B",
    "#FFDC00",
    "#2ECC40",
    "#0074D9",
    "#B10DC9",
    "#85144b",
    "#3D9970",
    "#39CCCC",
    "#01FF70",
    "#F012BE",
    "#7FDBFF",
    "#870C25",
    "#FF6B6B",
    "#4ECDC4",
    "#45B7D1",
    "#FDBB5A",
    "#6A5ACD",
    "#20B2AA",
    "#FF69B4",
    "#8A2BE2",
    "#00CED1",
    "#FF4500",
    "#32CD32",
    "#DC143C",
    "#FF7F50",
    "#00FF7F",
    "#1E90FF",
    "#FF00FF",
    "#00FF00",
    "#FF1493",
    "#00BFFF",
    "#1E1E1E",
    "#8B4513",
    "#4682B4",
    "#D2691E",
    "#5F9EA0",
    "#7B68EE",
    "#00FA9A",
    "#FF4500",
    "#2F4F4F",
    "#DAA520",
    "#008B8B",
    "#9400D3",
    "#FF6347",
    "#40E0D0",
    "#EE82EE",
    "#F0E68C",
    "#D8BFD8",
    "#FF9999",
    "#66CDAA"
  )
  

  
  strip_colors <- c("#E6F2FF", "#CCE5FF", "#B3D9FF") 
  
  # Get unique combinations of Settlement_Dist and BP_Dist
  combinations <- agg_data_1 %>%
    select(Settlement_Dist, BP_Dist) %>%
    distinct()
  
  # Create a directory for plots if it doesn't exist
  dir.create("SPD_plots", showWarnings = FALSE)
  
  # For each combination
  for(i in 1:nrow(combinations)) {
    combo <- combinations[i,]
    
    # Get all data for this combination
    plot_data <- agg_data_1 %>%
      filter(Settlement_Dist == combo$Settlement_Dist,
             BP_Dist == combo$BP_Dist,
             Dataset %in% c("Original", "Sampled"),
             Sampling_Pct %in% c(20, 30, 40)) %>%
      group_by(calBP, Dataset, Sample_Seed, Sampling_Pct) %>%
      summarise(plot_density = sum(Adjusted_PrDens), .groups = "drop")
    
    # Separate original and sampled data
    original_data <- plot_data %>% 
      filter(Dataset == "Original") %>%
      mutate(line_type = "Original SPD") %>%
      distinct(calBP, plot_density, Sampling_Pct, line_type)
    
    sampled_data <- plot_data %>% 
      filter(Dataset == "Sampled")
    
    # Create plot
    p <- ggplot() +
      geom_line(data = sampled_data,
                aes(x = calBP, y = plot_density, 
                    group = interaction(Sample_Seed, Sampling_Pct),
                    color = factor(Sample_Seed)),
                linewidth = 1.25,
                alpha = 0.7,
                show.legend = FALSE) +
      geom_line(data = original_data,
                aes(x = calBP, y = plot_density, 
                    linetype = line_type),
                color = "black",
                linewidth = 3) +
      facet_grid(Sampling_Pct ~ ., 
                 labeller = labeller(Sampling_Pct = function(x) paste0(x, "% Sampling"))) +
      scale_x_reverse(limits = c(5500, 3000),
                      breaks = seq(5500, 3000, -200)) +
      scale_color_manual(values = distinct_colors) +
      scale_linetype_manual(values = c("Original SPD" = "solid")) +
      labs(title = sprintf("Original - Sampled Comparison\n%s - %s",
                           combo$Settlement_Dist,
                           combo$BP_Dist),
           x = "Years BP",
           y = "Probability Density") +
      theme_bw() +
      theme(
        axis.text.x = element_text(size =32, face = "bold"),
        axis.text.y = element_text(size = 32, face = "bold"),
        axis.title.x = element_text(size = 30,
                                    face = "bold",
                                    margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 30,
                                    face = "bold",
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_text(size = 30, face = "bold"),
        legend.text = element_text(size = 30),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 20, r = 30, b = 20, l = 30),
        strip.text = element_text(size = 35, face = "bold", color = "black"),
        strip.background = element_rect(fill = strip_colors)
      ) +
      guides(linetype = guide_legend(title = "Line Type"),
             color = "none")
    
    filename <- sprintf("SPD_plots/Original_vs_Sampled_%s_%s_vertically_stacked.png",
                        combo$Settlement_Dist,
                        combo$BP_Dist)
    
    ggsave(filename, p, width = 20, height = 20, dpi = 300)
  }
}

wide_plot <- create_original_vs_sampled_plots_wide(agg_data_1)



create_combined_original_vs_sampled_spd_plots <- function(agg_data_1) {
  strip_colors <- c("#E6F2FF", "#CCE5FF", "#B3D9FF") 
  
  agg_data_1 <- agg_data_1 %>%
    mutate(
      Settlement_Dist = ifelse(Settlement_Dist == "power_law", "power law", Settlement_Dist),
      Combined_Dist = paste(Settlement_Dist, BP_Dist, sep = "-")
    )
  
  # Prepare the data
  plot_data <- agg_data_1 %>%
    filter(Dataset %in% c("Original", "Sampled")) %>%
    group_by(Combined_Dist, calBP, Dataset, Sample_Seed, Sampling_Pct) %>%
    summarise(plot_density = sum(Adjusted_PrDens), .groups = "drop")
  
  # Separate original and sampled data
  original_data <- plot_data %>% 
    filter(Dataset == "Original") %>%
    group_by(Combined_Dist, calBP, Sampling_Pct) %>%
    summarise(plot_density = mean(plot_density), .groups = "drop") %>%
    mutate(line_type = "Original SPD")
  
  sampled_data <- plot_data %>% 
    filter(Dataset == "Sampled")
  
  # Generate 50 distinct, vibrant colors
  distinct_colors <- c(
    "#FF4136",
    "#FF851B",
    "#FFDC00",
    "#2ECC40",
    "#0074D9",
    "#B10DC9",
    "#85144b",
    "#3D9970",
    "#39CCCC",
    "#01FF70",
    "#F012BE",
    "#7FDBFF",
    "#870C25",
    "#FF6B6B",
    "#4ECDC4",
    "#45B7D1",
    "#FDBB5A",
    "#6A5ACD",
    "#20B2AA",
    "#FF69B4",
    "#8A2BE2",
    "#00CED1",
    "#FF4500",
    "#32CD32",
    "#DC143C",
    "#FF7F50",
    "#00FF7F",
    "#1E90FF",
    "#FF00FF",
    "#00FF00",
    "#FF1493",
    "#00BFFF",
    "#1E1E1E",
    "#8B4513",
    "#4682B4",
    "#D2691E",
    "#5F9EA0",
    "#7B68EE",
    "#00FA9A",
    "#FF4500",
    "#2F4F4F",
    "#DAA520",
    "#008B8B",
    "#9400D3",
    "#FF6347",
    "#40E0D0",
    "#EE82EE",
    "#F0E68C",
    "#D8BFD8",
    "#FF9999",
    "#66CDAA"
  )
  
  # Create combined dataset
  combined_data <- left_join(
    sampled_data, 
    original_data, 
    by = c("Combined_Dist", "calBP", "Sampling_Pct"),
    suffix = c("_sampled", "_original")
  )
  
  # Create plot
  p <- ggplot(combined_data) +
    geom_line(aes(x = calBP, y = plot_density_sampled, 
                  group = interaction(Sample_Seed, Sampling_Pct),
                  color = factor(Sample_Seed)),
              linewidth = 1.5,
              alpha = 0.9,
              show.legend = FALSE) +
    # Add original line with legend
    geom_line(data = original_data,
              aes(x = calBP, y = plot_density, 
                  linetype = line_type),
              color = "black",
              linewidth = 2.5) +
    facet_grid(Combined_Dist ~ Sampling_Pct, 
               labeller = labeller(
                 Sampling_Pct = function(x) paste0(x, "%")
               ),
               scales = "free_y") +
    scale_x_reverse(limits = c(5500, 3000),
                    breaks = seq(5500, 3000, -200)) +
    scale_color_manual(values = distinct_colors) +
    scale_linetype_manual(values = c("Original SPD" = "solid")) +
    labs(title = "SPD Comparison Across Distributions",
         x = "Years BP",
         y = "Probability Density") +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 40, face = "bold", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 40, face = "bold"),
      axis.title.x = element_text(size = 40,
                                  face = "bold",
                                  margin = margin(t = 15, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(size = 40,
                                  face = "bold",
                                  margin = margin(t = 0, r = 15, b = 0, l = 0)),
      plot.title = element_text(size = 45, face = "bold", hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_text(size = 40, face = "bold"),
      legend.text = element_text(size = 40, face="bold"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 30, r = 30, b = 30, l = 30),
      strip.text.x = element_text(size = 45, face = "bold", color = "black"), 
      strip.text.y = element_text(size = 40, face = "bold", color = "black"), 
      strip.background.x = element_rect(fill = strip_colors),
      strip.background.y = element_rect(fill = "white"), legend.key.width = unit(3, "cm"),
    ) +
    guides(
      linetype = guide_legend(title = NULL),  
      color = "none",
      linewidth= 12
    )
  
  ggsave("Combined_SPD_Plots.png", p, 
         width = 35, height = 45, dpi = 300)
  
  return(p)
}


combined_plot <- create_combined_original_vs_sampled_spd_plots(agg_data_1)



# For Original vs Resampled30
create_original_vs_resampled_plots <- function(agg_data_1) {
  # Get unique combinations of Settlement_Dist and BP_Dist
  combinations <- agg_data_1 %>%
    select(Settlement_Dist, BP_Dist) %>%
    distinct()
  
  # For each combination
  for(i in 1:nrow(combinations)) {
    combo <- combinations[i,]
    
    # Get all data for this combination
    plot_data <- agg_data_1 %>%
      filter(Settlement_Dist == combo$Settlement_Dist,
             BP_Dist == combo$BP_Dist,
             Dataset %in% c("Original", "Resampled30")) %>%
      group_by(calBP, Dataset, Sample_Seed, Sampling_Pct) %>%
      summarise(plot_density = sum(Adjusted_PrDens), .groups = "drop")
    
    # Create plot
    p <- ggplot() +
      # Add resampled lines with continuous color scale
      geom_line(data = plot_data %>% filter(Dataset == "Resampled30"),
                aes(x = calBP, y = plot_density, 
                    group = Sample_Seed,
                    color = as.numeric(Sample_Seed)),
                linewidth = 0.3,
                alpha = 0.5) +
      # Add original line
      geom_line(data = plot_data %>% filter(Dataset == "Original"),
                aes(x = calBP, y = plot_density),
                color = "#FF4B00",
                linewidth = 2.5) +
      facet_wrap(~Sampling_Pct, 
                 labeller = labeller(Sampling_Pct = function(x) paste0(x, "% Sampling")),
                 nrow = 1) +
      scale_x_reverse(limits = c(5500, 3000),
                      breaks = seq(5500, 3000, -200))+
      scale_color_viridis_c(
        name = "Sampling Seed",
        option = "mako",  # Blue scale
        guide = guide_colorbar(
          frame.colour = "black",
          ticks.colour = "black",
          barheight = 0.3,
          barwidth = 10
        )
      ) +
      labs(title = sprintf("Original - Resampled30 Comparison\n%s - %s",
                           combo$Settlement_Dist,
                           combo$BP_Dist),
           x = "Years BP",
           y = "Probability Density") +
      theme_bw() +
      theme(
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18,
                                    face = "bold",
                                    margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 18,
                                    face = "bold",
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 13),
        legend.key.width = unit(2, "cm"),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(t = 10, b = 10),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
        strip.text = element_text(size = 25, face = "bold"),
        strip.background = element_rect(fill = "white")
      )
    
    # Save plot
    filename <- sprintf("SPD_plots/Original_vs_Resampled30_%s_%s_allseeds_faceted.png",
                        combo$Settlement_Dist,
                        combo$BP_Dist)
    
    ggsave(filename, p, width = 20, height = 8, dpi = 300)
  }
}

# For Original vs Weighted
create_original_vs_weighted_plots <- function(agg_data_1) {
  # Get unique combinations of Settlement_Dist and BP_Dist
  combinations <- agg_data_1 %>%
    select(Settlement_Dist, BP_Dist) %>%
    distinct()
  
  # For each combination
  for(i in 1:nrow(combinations)) {
    combo <- combinations[i,]
    
    # Get all data for this combination
    plot_data <- agg_data_1 %>%
      filter(Settlement_Dist == combo$Settlement_Dist,
             BP_Dist == combo$BP_Dist,
             Dataset %in% c("Original", "Weighted Sampled")) %>%
      group_by(calBP, Dataset, Sample_Seed, Sampling_Pct) %>%
      summarise(plot_density = sum(Adjusted_PrDens), .groups = "drop")
    
    # Create plot
    p <- ggplot() +
      # Add weighted sampled lines with continuous color scale
      geom_line(data = plot_data %>% filter(Dataset == "Weighted Sampled"),
                aes(x = calBP, y = plot_density, 
                    group = Sample_Seed,
                    color = as.numeric(Sample_Seed)),
                linewidth = 0.3,
                alpha = 0.5) +
      # Add original line
      geom_line(data = plot_data %>% filter(Dataset == "Original"),
                aes(x = calBP, y = plot_density),
                color = "#FF4B00",
                linewidth = 2.5) +
      facet_wrap(~Sampling_Pct, 
                 labeller = labeller(Sampling_Pct = function(x) paste0(x, "% Sampling")),
                 nrow = 1) +
      scale_x_reverse(limits = c(5500, 3000),
                      breaks = seq(5500, 3000, -200))+
      scale_color_viridis_c(
        name = "Sampling Seed",
        option = "cividis", 
        guide = guide_colorbar(
          frame.colour = "black",
          ticks.colour = "black",
          barheight = 0.3,
          barwidth = 10
        )
      ) +
      labs(title = sprintf("Original - Weighted Sampled Comparison\n%s - %s",
                           combo$Settlement_Dist,
                           combo$BP_Dist),
           x = "Years BP",
           y = "Probability Density") +
      theme_bw() +
      theme(
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 18,
                                    face = "bold",
                                    margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 18,
                                    face = "bold",
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 13),
        legend.key.width = unit(2, "cm"),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(t = 10, b = 10),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
        strip.text = element_text(size = 25, face = "bold"),
        strip.background = element_rect(fill = "white")
      )
    
    # Save plot
    filename <- sprintf("SPD_plots/Original_vs_WeightedSampled_%s_%s_allseeds_faceted.png",
                        combo$Settlement_Dist,
                        combo$BP_Dist)
    
    ggsave(filename, p, width = 20, height = 8, dpi = 300)
  }
}



create_original_vs_resampled_plots_wide <- function(agg_data_1) {
  # 50 distinct, vibrant colors
  distinct_colors <- c(
    "#FF4136", "#FF851B", "#FFDC00", "#2ECC40", "#0074D9", "#B10DC9", "#85144b", 
    "#3D9970", "#39CCCC", "#01FF70", "#F012BE", "#7FDBFF", "#870C25", "#FF6B6B", 
    "#4ECDC4", "#45B7D1", "#FDBB5A", "#6A5ACD", "#20B2AA", "#FF69B4", "#8A2BE2", 
    "#00CED1", "#FF4500", "#32CD32", "#DC143C", "#FF7F50", "#00FF7F", "#1E90FF", 
    "#FF00FF", "#00FF00", "#FF1493", "#00BFFF", "#1E1E1E", "#8B4513", "#4682B4", 
    "#D2691E", "#5F9EA0", "#7B68EE", "#00FA9A", "#FF4500", "#2F4F4F", "#DAA520", 
    "#008B8B", "#9400D3", "#FF6347", "#40E0D0", "#EE82EE", "#F0E68C", "#D8BFD8", 
    "#FF9999", "#66CDAA"
  )
  
  strip_colors <- c("#E6F2FF", "#CCE5FF", "#B3D9FF")  
  
  # Get unique combinations of Settlement_Dist and BP_Dist
  combinations <- agg_data_1 %>%
    select(Settlement_Dist, BP_Dist) %>%
    distinct()
  
  # Create a directory for plots if it doesn't exist
  dir.create("SPD_plots", showWarnings = FALSE)
  
  # For each combination
  for(i in 1:nrow(combinations)) {
    combo <- combinations[i,]
    
    # Get all data for this combination
    plot_data <- agg_data_1 %>%
      filter(Settlement_Dist == combo$Settlement_Dist,
             BP_Dist == combo$BP_Dist,
             Dataset %in% c("Original", "Resampled30"),
             Sampling_Pct %in% c(20, 30, 40)) %>%
      group_by(calBP, Dataset, Sample_Seed, Sampling_Pct) %>%
      summarise(plot_density = sum(Adjusted_PrDens), .groups = "drop")
    
    # Separate original and resampled data
    original_data <- plot_data %>% 
      filter(Dataset == "Original") %>%
      mutate(line_type = "Original SPD") %>%
      distinct(calBP, plot_density, Sampling_Pct, line_type)
    
    resampled_data <- plot_data %>% 
      filter(Dataset == "Resampled30")
    
    # Create plot
    p <- ggplot() +
      # Add resampled lines with distinct color palette
      geom_line(data = resampled_data,
                aes(x = calBP, y = plot_density, 
                    group = interaction(Sample_Seed, Sampling_Pct),
                    color = factor(Sample_Seed)),
                linewidth = 1.25,
                alpha = 0.7,
                show.legend = FALSE) +
      # Add original line with legend
      geom_line(data = original_data,
                aes(x = calBP, y = plot_density, 
                    linetype = line_type),
                color = "black",
                linewidth = 3) +
      facet_grid(Sampling_Pct ~ ., 
                 labeller = labeller(Sampling_Pct = function(x) paste0(x, "% Sampling"))) +
      scale_x_reverse(limits = c(5500, 3000),
                      breaks = seq(5500, 3000, -200)) +
      # Use the distinct color palette for resampled lines
      scale_color_manual(values = distinct_colors) +
      # Customize the linetype for original SPD
      scale_linetype_manual(values = c("Original SPD" = "solid")) +
      labs(title = sprintf("Original - Resampled30 Comparison\n%s - %s",
                           combo$Settlement_Dist,
                           combo$BP_Dist),
           x = "Years BP",
           y = "Probability Density") +
      theme_bw() +
      theme(
        axis.text.x = element_text(size = 32, face = "bold"),
        axis.text.y = element_text(size = 32, face = "bold"),
        axis.title.x = element_text(size = 30,
                                    face = "bold",
                                    margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 30,
                                    face = "bold",
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_text(size = 30, face = "bold"),
        legend.text = element_text(size = 30),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 20, r = 30, b = 20, l = 30),
        strip.text = element_text(size = 35, face = "bold", color = "black"),
        strip.background = element_rect(fill = strip_colors)
      ) +
      # Only show legend for original SPD
      guides(linetype = guide_legend(title = "Line Type"),
             color = "none")
    
    # Save plot
    filename <- sprintf("SPD_plots/Original_vs_Resampled30_%s_%s_vertically_stacked.png",
                        combo$Settlement_Dist,
                        combo$BP_Dist)
    
    ggsave(filename, p, width = 20, height = 20, dpi = 300)
  }
}

create_original_vs_weighted_plots_wide <- function(agg_data_1) {
  # 50 distinct, vibrant colors
  distinct_colors <- c(
    "#FF4136", "#FF851B", "#FFDC00", "#2ECC40", "#0074D9", "#B10DC9", "#85144b", 
    "#3D9970", "#39CCCC", "#01FF70", "#F012BE", "#7FDBFF", "#870C25", "#FF6B6B", 
    "#4ECDC4", "#45B7D1", "#FDBB5A", "#6A5ACD", "#20B2AA", "#FF69B4", "#8A2BE2", 
    "#00CED1", "#FF4500", "#32CD32", "#DC143C", "#FF7F50", "#00FF7F", "#1E90FF", 
    "#FF00FF", "#00FF00", "#FF1493", "#00BFFF", "#1E1E1E", "#8B4513", "#4682B4", 
    "#D2691E", "#5F9EA0", "#7B68EE", "#00FA9A", "#FF4500", "#2F4F4F", "#DAA520", 
    "#008B8B", "#9400D3", "#FF6347", "#40E0D0", "#EE82EE", "#F0E68C", "#D8BFD8", 
    "#FF9999", "#66CDAA"
  )
  
  strip_colors <- c("#E6F2FF", "#CCE5FF", "#B3D9FF")  
  
  # Get unique combinations of Settlement_Dist and BP_Dist
  combinations <- agg_data_1 %>%
    select(Settlement_Dist, BP_Dist) %>%
    distinct()
  
  # Create a directory for plots if it doesn't exist
  dir.create("SPD_plots", showWarnings = FALSE)
  
  # For each combination
  for(i in 1:nrow(combinations)) {
    combo <- combinations[i,]
    
    # Get all data for this combination
    plot_data <- agg_data_1 %>%
      filter(Settlement_Dist == combo$Settlement_Dist,
             BP_Dist == combo$BP_Dist,
             Dataset %in% c("Original", "Weighted Sampled"),
             Sampling_Pct %in% c(20, 30, 40)) %>%
      group_by(calBP, Dataset, Sample_Seed, Sampling_Pct) %>%
      summarise(plot_density = sum(Adjusted_PrDens), .groups = "drop")
    
    # Separate original and weighted sampled data
    original_data <- plot_data %>% 
      filter(Dataset == "Original") %>%
      mutate(line_type = "Original SPD") %>%
      distinct(calBP, plot_density, Sampling_Pct, line_type)
    
    weighted_data <- plot_data %>% 
      filter(Dataset == "Weighted Sampled")
    
    # Create plot
    p <- ggplot() +
      # Add weighted sampled lines with distinct color palette
      geom_line(data = weighted_data,
                aes(x = calBP, y = plot_density, 
                    group = interaction(Sample_Seed, Sampling_Pct),
                    color = factor(Sample_Seed)),
                linewidth = 1.25,
                alpha = 0.7,
                show.legend = FALSE) +
      # Add original line with legend
      geom_line(data = original_data,
                aes(x = calBP, y = plot_density, 
                    linetype = line_type),
                color = "black",
                linewidth = 3) +
      facet_grid(Sampling_Pct ~ ., 
                 labeller = labeller(Sampling_Pct = function(x) paste0(x, "% Sampling"))) +
      scale_x_reverse(limits = c(5500, 3000),
                      breaks = seq(5500, 3000, -200)) +
      # Use the distinct color palette for weighted sampled lines
      scale_color_manual(values = distinct_colors) +
      # Customize the linetype for original SPD
      scale_linetype_manual(values = c("Original SPD" = "solid")) +
      labs(title = sprintf("Original - Weighted Sampled Comparison\n%s - %s",
                           combo$Settlement_Dist,
                           combo$BP_Dist),
           x = "Years BP",
           y = "Probability Density") +
      theme_bw() +
      theme(
        axis.text.x = element_text(size = 32, face = "bold"),
        axis.text.y = element_text(size = 32, face = "bold"),
        axis.title.x = element_text(size = 30,
                                    face = "bold",
                                    margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 30,
                                    face = "bold",
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_text(size = 30, face = "bold"),
        legend.text = element_text(size = 30),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 20, r = 30, b = 20, l = 30),
        strip.text = element_text(size = 35, face = "bold", color = "black"),
        strip.background = element_rect(fill = strip_colors)
      ) +
      # Only show legend for original SPD
      guides(linetype = guide_legend(title = "Line Type"),
             color = "none")
    
    # Save plot
    filename <- sprintf("SPD_plots/Original_vs_WeightedSampled_%s_%s_vertically_stacked.png",
                        combo$Settlement_Dist,
                        combo$BP_Dist)
    
    ggsave(filename, p, width = 20, height = 20, dpi = 300)
  }
}

create_original_vs_sampled_plots(agg_data_1)
create_original_vs_resampled_plots(agg_data_1)   
create_original_vs_weighted_plots(agg_data_1)    

create_original_vs_resampled_plots_wide(agg_data_1)    
create_original_vs_weighted_plots_wide(agg_data_1)      

### Merge all into one
create_combined_comparison_plots_wide <- function(agg_data_1) {
  # 50 distinct, vibrant colors
  distinct_colors <- c(
    "#FF4136", "#FF851B", "#FFDC00", "#2ECC40", "#0074D9", "#B10DC9", "#85144b", 
    "#3D9970", "#39CCCC", "#01FF70", "#F012BE", "#7FDBFF", "#870C25", "#FF6B6B", 
    "#4ECDC4", "#45B7D1", "#FDBB5A", "#6A5ACD", "#20B2AA", "#FF69B4", "#8A2BE2", 
    "#00CED1", "#FF4500", "#32CD32", "#DC143C", "#FF7F50", "#00FF7F", "#1E90FF", 
    "#FF00FF", "#00FF00", "#FF1493", "#00BFFF", "#1E1E1E", "#8B4513", "#4682B4", 
    "#D2691E", "#5F9EA0", "#7B68EE", "#00FA9A", "#FF4500", "#2F4F4F", "#DAA520", 
    "#008B8B", "#9400D3", "#FF6347", "#40E0D0", "#EE82EE", "#F0E68C", "#D8BFD8", 
    "#FF9999", "#66CDAA"
  )
  
  # Get unique combinations of Settlement_Dist and BP_Dist
  combinations <- agg_data_1 %>%
    select(Settlement_Dist, BP_Dist) %>%
    distinct()
  
  # For each combination of Settlement_Dist and BP_Dist
  for(i in 1:nrow(combinations)) {
    combo <- combinations[i,]
    
    # Get all data for this combination
    plot_data <- agg_data_1 %>%
      filter(Settlement_Dist == combo$Settlement_Dist,
             BP_Dist == combo$BP_Dist,
             Dataset %in% c("Original", "Sampled", "Weighted Sampled", "Resampled30"),
             Sampling_Pct %in% c(20, 30, 40)) %>%
      group_by(calBP, Dataset, Sample_Seed, Sampling_Pct) %>%
      summarise(plot_density = sum(Adjusted_PrDens), .groups = "drop")
    
    # Create three separate plots
    # 1. Original vs Sampled
    p1 <- ggplot() +
      geom_line(data = plot_data %>% filter(Dataset == "Sampled"),
                aes(x = calBP, y = plot_density, 
                    group = interaction(Sample_Seed, Sampling_Pct),
                    color = factor(Sample_Seed)),
                linewidth = 1.25,
                alpha = 0.7,
                show.legend = FALSE) +
      geom_line(data = plot_data %>% filter(Dataset == "Original"),
                aes(x = calBP, y = plot_density),
                color = "black",
                linewidth = 3) +
      facet_grid(rows = vars("Original-Sampled"),
                 cols = vars(Sampling_Pct),
                 labeller = labeller(
                   Sampling_Pct = function(x) paste0(x, "%"),
                   "Original-Sampled" = as_labeller(function(x) "Original-Sampled"))) +
      scale_color_manual(values = distinct_colors) +
      scale_x_reverse(limits = c(5500, 3000),
                      breaks = seq(5500, 3000, -200))
    
    # 2. Original vs Resampled30
    p2 <- ggplot() +
      geom_line(data = plot_data %>% filter(Dataset == "Resampled30"),
                aes(x = calBP, y = plot_density, 
                    group = interaction(Sample_Seed, Sampling_Pct),
                    color = factor(Sample_Seed)),
                linewidth = 1.25,
                alpha = 0.7,
                show.legend = FALSE) +
      geom_line(data = plot_data %>% filter(Dataset == "Original"),
                aes(x = calBP, y = plot_density),
                color = "black",
                linewidth = 3) +
      facet_grid(rows = vars("Original-Resampled30"),
                 cols = vars(Sampling_Pct),
                 labeller = labeller(
                   Sampling_Pct = function(x) paste0(x, "%"),
                   "Original-Resampled30" = as_labeller(function(x) "Original-Resampled30"))) +
      scale_color_manual(values = distinct_colors) +
      scale_x_reverse(limits = c(5500, 3000),
                      breaks = seq(5500, 3000, -200))
    
    # 3. Original vs Weighted Sampled
    p3 <- ggplot() +
      geom_line(data = plot_data %>% filter(Dataset == "Weighted Sampled"),
                aes(x = calBP, y = plot_density, 
                    group = interaction(Sample_Seed, Sampling_Pct),
                    color = factor(Sample_Seed)),
                linewidth = 1.25,
                alpha = 0.7,
                show.legend = FALSE) +
      geom_line(data = plot_data %>% filter(Dataset == "Original"),
                aes(x = calBP, y = plot_density),
                color = "black",
                linewidth = 3) +
      facet_grid(rows = vars("Original-Weighted Sampled"),
                 cols = vars(Sampling_Pct),
                 labeller = labeller(
                   Sampling_Pct = function(x) paste0(x, "%"),
                   "Original-Weighted Sampled" = as_labeller(function(x) "Original-Weighted Sampled"))) +
      scale_color_manual(values = distinct_colors) +
      scale_x_reverse(limits = c(5500, 3000),
                      breaks = seq(5500, 3000, -200))
    
    plots <- list(p1, p2, p3) %>%
      map(~ . +
            labs(x = "Years BP",
                 y = "Probability Density") +
            theme_bw() +
            theme(
              # Axis text
              axis.text.x = element_text(size = 35, face = "bold", angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_text(size = 35, face = "bold"),
              
              # Axis titles
              axis.title.x = element_text(size = 35,
                                          face = "bold",
                                          margin = margin(t = 20, r = 0, b = 0, l = 0)),
              axis.title.y = element_text(size = 35,
                                          face = "bold",
                                          margin = margin(t = 0, r = 20, b = 0, l = 0)),
              
              strip.text.x = element_text(size = 40, face = "bold"), 
              strip.text.y = element_text(size = 35, face = "bold"), 
              strip.background.x = element_rect(fill = "white"),
              strip.background.y = element_rect(fill = "white"),
              
              legend.position = "none",
              panel.grid.minor = element_blank(),
              plot.margin = margin(t = 30, r = 30, b = 30, l = 30)
            )
      )
    
    combined_plot <- wrap_plots(plots, ncol = 1) +
      plot_annotation(
        title = sprintf("SPD Comparisons\n%s - %s",
                        combo$Settlement_Dist,
                        combo$BP_Dist),
        theme = theme(
          plot.title = element_text(size = 40, face = "bold", hjust = 0.5)
        )
      )
      
    
    # Save plot
    filename <- sprintf("SPD_plots/Combined_Comparisons_%s_%s_wide.png",
                        combo$Settlement_Dist,
                        combo$BP_Dist)
    
    ggsave(filename, combined_plot, width = 40, height = 33, dpi = 300)
  }
}

create_combined_comparison_plots_wide(agg_data_1)




############ STEP 5 (Figure 6a&b): rolling mean SPDs by settlement###############
guides(
  fill = guide_legend(nrow = 1),
  color = guide_legend(nrow = 1)
)


create_settlement_plot_from_file_smooth_RM50 <- function(seed, settlement_dist, bp_dist, file_path = "", smooth_k = 50) {
  # Format file name with path
  filename <- file.path(file_path, sprintf("normalized_S%d_RS1_%s_%s.csv", 
                                           seed, 
                                           settlement_dist,
                                           bp_dist))
  
  cat(sprintf("\nReading file: %s\n", filename))
  
  # Read data
  all_plot_df <- read.csv(filename)
  
  # Create plot
  n_settlements <- length(unique(all_plot_df$Settlement))
  distinct_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
                       "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5")
  
  if(n_settlements > length(distinct_colors)) {
    distinct_colors <- colorRampPalette(distinct_colors)(n_settlements)
  }
  
  # --- Apply smoothing using rolling mean ---
  all_plot_df <- all_plot_df %>%
    filter(Dataset %in% c("Original", "Sampled", "Weighted Sampled", "Resampled30")) %>%
    mutate(Dataset = factor(Dataset, 
                            levels = c("Original", "Sampled", 
                                       "Weighted Sampled", "Resampled30"))) %>%
    group_by(Settlement, Dataset) %>%
    arrange(calBP, .by_group = TRUE) %>%
    mutate(
      Adjusted_PrDens_smooth = zoo::rollmean(Adjusted_PrDens, k = smooth_k, fill = 0, align = "center"),
      Adjusted_PrDens_smooth = Adjusted_PrDens_smooth / sum(Adjusted_PrDens_smooth) * sum(Adjusted_PrDens)
    ) %>%
    ungroup()
  
  # Create plot
  p <- all_plot_df %>%
    ggplot(aes(x = calBP, y = Adjusted_PrDens_smooth, fill = factor(Settlement), 
               color = factor(Settlement))) +
    geom_area(alpha = 0.4, position = "identity") +
    geom_line(linewidth = 0.5) +
    facet_wrap(~Dataset, nrow = 2,
               labeller = as_labeller(c(
                 "Original" = "Hypothetical Population",
                 "Sampled" = "Random Sample",
                 "Weighted Sampled" = "Weighted",
                 "Resampled30" = "Bootstrapped"
               )))+
    scale_fill_manual(name = "Settlement", values = distinct_colors) +
    scale_color_manual(name = "Settlement", values = distinct_colors) +
    guides(
      fill = guide_legend(nrow = 1),
      color = guide_legend(nrow = 1)
    ) +
    labs(title = "Settlement SPD Comparison",
         subtitle = sprintf("%s - %s (Seed %d)", 
                            settlement_dist, bp_dist, seed),
         x = "Cal BP",
         y = "Probability Density") +
    scale_x_reverse() +
    theme_bw() +
    theme(
      plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 40, face = "bold", hjust = 0.5, 
                                   margin = margin(t = 10, b = 20)),
      axis.title.x = element_text(size = 35, face = "bold", margin = margin(t = 15)),
      axis.title.y = element_text(size = 35, face = "bold", margin = margin(r = 15)),
      axis.text.x = element_text(size = 35, face = "bold"),
      axis.text.y = element_text(size = 35, face = "bold"),
      strip.text = element_text(size = 40, face = "bold"),
      legend.title = element_text(size = 40, face = "bold"),
      legend.text = element_text(size = 40, face="bold"),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = margin(t = 5, b = 0, l = 0, r = 0),
      legend.spacing.x = unit(1.5, 'cm'),
      legend.key.width = unit(1.5, "cm"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      panel.spacing = unit(1, "cm"),
      plot.margin = margin(t = 20, r = 40, b = 20, l = 20)
    )
  
  # Define output filename with path
  output_filename <- sprintf("Settlement_SPD_S%d_%s_%s.png",
                                                  seed,
                                                  settlement_dist,
                                                  bp_dist)
  
  # Save plot
  ggsave(output_filename, p, width = 30, height = 16)
  
  return(p)
}

  



# Create plot for specific combination(Figure 6a)
plot6a <- create_settlement_plot_from_file_smooth_RM50(
  seed = 45,
  settlement_dist = "power_law",
  bp_dist = "uniform",
  file_path= ""
)

# Create plot for specific combination(Figure 6a)
plot6b <- create_settlement_plot_from_file_smooth_RM50(
  seed = 50,
  settlement_dist = "power_law",
  bp_dist = "normal",
  file_path= ""
)



# Process multiple seeds
seeds <- c(1:50)
for(seed in seeds) {
  plot <- create_settlement_plot_from_file_smooth_RM50(
    seed = seed,
    settlement_dist = "power_law",
    bp_dist = "skewed"
  )
}


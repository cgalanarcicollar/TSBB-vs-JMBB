library(dplyr)
library(ggplot2)
library(purrr)

# --- Load files  ---
phi_labels <- c("005", "05", "1")
file_paths <- sprintf("Results/S2_p%sa_0.rds", phi_labels)

plots <- list()
summary_df_list <- list()

# --- Summary of longitudinal frequencies and survival events ---

for (i in seq_along(file_paths)) {
  path <- file_paths[i]
  phi <- phi_labels[i]
  
  result_list <- readRDS(path)

  summary_df <- expand.grid(i = 1:50, j = 1:10) %>%
    pmap_dfr(function(i, j) {
      data_long <- result_list$Data[[i]][[j]]$longitudinal
      data_surv <- result_list$Data[[i]][[j]]$survival
      freqs <- tabulate(data_long$id)
      
      tibble(
        freq_1 = sum(freqs == 1),
        freq_2 = sum(freqs == 2),
        freq_3 = sum(freqs == 3),
        freq_4 = sum(freqs == 4),
        total_meas = nrow(data_long),
        events = sum(data_surv$event)
      )
    })
  
  summary_df_list[[phi]] <- summary_df
  
  mean_freqs <- colMeans(summary_df[, c("freq_1", "freq_2", "freq_3", "freq_4")])
  record_meas_mean <- mean(summary_df$total_meas)
  events_summary <- summary_df$events
  
  meas_df <- tibble(
    Measurements = 1:4,
    Frequency = mean_freqs
  )
  
  q <- ggplot(meas_df, aes(x = Measurements, y = Frequency)) +
    geom_bar(stat = "identity", fill = "#3399FF", color = "#e9ecef", alpha = 0.9) +
    scale_x_continuous(breaks = 1:4) +
    ylim(0, 420) +
    theme_bw() +
    xlab("measurement") +
    ylab("frequency") +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))

  plots[[phi]] <- q
}

print(plots[["005"]])  # Plot for phi = 005 (It is the same for all phi)


# --- Summary of data for each phi ---
combined_summary <- bind_rows(
  lapply(names(summary_df_list), function(phi) {
    summary_df <- summary_df_list[[phi]]
    tibble(
      phi = phi,
      mean_events = mean(summary_df$events),
      mean_recorded_measurements = mean(summary_df$total_meas),
      mean_freq_1 = mean(summary_df$freq_1)
    )
  })
)


print(combined_summary)

col_means <- colMeans(combined_summary[, sapply(combined_summary, is.numeric)])

print(col_means)

library(dplyr)
library(ggplot2)
library(purrr)

# --- Load S1 files ---
phi_labels <- c("005", "05", "1")
file_paths <- sprintf("Results/S1_p%sa0.rds", phi_labels)

# --- Collect the first simulated data and plot baseline measurement ---
plots <- map2(file_paths, phi_labels, function(path, phi_label) {
  result_list <- readRDS(path)
  
  data_long <- result_list$Data[[1]][[1]]$longitudinal
  baseline <- data_long[!duplicated(data_long$id), ]
  y_vals <- as.numeric(unlist(baseline$y))
  
  ggplot(data.frame(y = y_vals), aes(x = y)) +
    geom_histogram(binwidth = 1, fill = "#69b3a2", color = "#e9ecef", alpha = 0.9) +
    scale_x_continuous(breaks = seq(0, 24, by = 1)) +
    theme_bw() +
    labs(
      x = "Measurement",
      y = "Frequency"
    ) +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    )
})

print(plots[[1]])
print(plots[[2]])
print(plots[[3]])


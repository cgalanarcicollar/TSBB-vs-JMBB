library(dplyr)
library(ggplot2)
library(purrr)
library(patchwork)

# For loading and processing
load_sim_data <- function(prefix, folder, alpha_values, phi_values, alpha_map = identity, filename_alpha_format = "a%s") {
  param_grid <- expand.grid(alpha = alpha_values, phi = phi_values, stringsAsFactors = FALSE)
  
  pmap_dfr(param_grid, function(alpha, phi) {
    phi_label <- ifelse(phi == 0.05, "005", ifelse(phi == 0.5, "05", "1"))
    alpha_file_val <- alpha_map(alpha)  # e.g., 15 for -1.5
    file_path <- sprintf("%s/%s_p%s%s.rds", folder, prefix, phi_label, sprintf(filename_alpha_format, alpha_file_val))
    
    if (!file.exists(file_path)) {
      stop(paste("File does not exist:", file_path))
    }
    
    data <- readRDS(file_path)
    
    bind_rows(
      tibble(est = data$result_JM$surv$Alpha[,1], alpha = alpha, phi = phi, model = "JMBB"),
      tibble(est = data$result_TS$surv$Alpha[,1], alpha = alpha, phi = phi, model = "TSBB")
    )
  })
}


# For S1
S1 <- load_sim_data(
  prefix = "S1",
  folder = "Results",
  alpha_values = c(0, 2, 4),
  phi_values = c(0.05, 0.5, 1),
  filename_alpha_format = "a%s"
)

# For S2
alpha_map_s2 <- function(a) {
  if (a == -1.5) return("15")
  if (a == -3) return("3")
  return("0")
}

S2 <- load_sim_data(
  prefix = "S2",
  folder = "Results",
  alpha_values = c(0, -1.5, -3),
  phi_values = c(0.05, 0.5, 1),
  alpha_map = alpha_map_s2,
  filename_alpha_format = "a_%s"
)

# ---- %Bias calculation   ----
process_data <- function(df, alpha_levels, alpha_labels) {
  df %>%
    mutate(
      phi = as.factor(phi),
      est = as.numeric(est),
      bias = if_else(alpha == 0, est - alpha, (est - alpha) / alpha),
      Alpha = factor(alpha, levels = alpha_levels),
      model = factor(model, levels = c("JMBB", "TSBB"))
    )
}


S1_proc <- process_data(S1, c(0, 2, 4), c("0", "2", "4"))
S2_proc <- process_data(S2, c(0, -1.5, -3), c("0", "-1.5", "-3"))

# ---- Plot S1 ----
q1 <- S1_proc %>%
  mutate(Alpha = factor(Alpha,
                        levels = c(0, 2, 4),
                        labels = c(expression(paste(alpha, " = 0")), 
                                   expression(paste(alpha, " = 2")),
                                   expression(paste(alpha, " = 4"))))) %>%
  ggplot(aes(x = phi, y = bias, fill = model)) +
  geom_boxplot() +
  facet_grid(cols = vars(Alpha), labeller = label_parsed) +
  geom_hline(yintercept = 0, color = "red") +
  scale_fill_manual(values = c(JMBB = 'grey74', TSBB = 'white')) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  xlab(expression(Phi)) 

# ---- Plot S2 ----
q2 <- S2_proc %>%
  mutate(Alpha = factor(Alpha,
                        levels = c(0, -1.5, -3),
                        labels = c(expression(paste(alpha, " = 0")), 
                                   expression(paste(alpha, " = -1.5")),
                                   expression(paste(alpha, " = -3"))))) %>%
  ggplot(aes(x = phi, y = bias, fill = model)) +
  geom_boxplot() +
  facet_grid(cols = vars(Alpha), labeller = label_parsed) +
  geom_hline(yintercept = 0, color = "red") +
  scale_fill_manual(values = c(JMBB = 'grey74', TSBB = 'white')) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  xlab(expression(Phi)) 

# ---- Combine plots ----
q1 + q2 + plot_layout(guides = 'collect')

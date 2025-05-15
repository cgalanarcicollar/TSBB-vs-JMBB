library(xtable)
library(dplyr)


proportion <- function(b, se, true, level = .95, df = Inf) {
  qtile <- level + (1 - level)/2
  lower.bound <- b - qt(qtile, df = df) * se
  upper.bound <- b + qt(qtile, df = df) * se
  in.ci <- ifelse(true >= lower.bound & true <= upper.bound, 1, 0)
  cp <- mean(in.ci)
  return(cp * 100)
}

# --- Load S2 files to get results ---
rds_files <- list.files("Results", pattern = "^S2_.*\\.rds$", full.names = TRUE)
summary_tables <- list()

for (file in rds_files) {
  out <- readRDS(file)
  dataset_name <- gsub("\\.rds$", "", basename(file))
  alpha <- out$parameters[[1]]
  
  # --- JM  ---
  Alpha_JM <- apply(out$result_JM$surv$Alpha, 2, mean)
  asd_alphJM <- Alpha_JM[3]
  esd_alphJM <- sd(out$result_JM$surv$Alpha[, 1])
  bias_alphJM <- if (alpha == 0) Alpha_JM[1] - alpha else (Alpha_JM[1] - alpha) / alpha
  CP_JM <- mean(out$result_JM$surv$Alpha[, 4] < alpha & alpha < out$result_JM$surv$Alpha[, 8]) * 100
  
  # --- TS  ---
  Alpha_TS <- apply(out$result_TS$survival$Alpha, 2, mean, na.rm = TRUE)
  asd_alphTS <- Alpha_TS[2]
  esd_alphTS <- sd(out$result_TS$survival$Alpha[, 1], na.rm = TRUE)
  bias_alphTS <- if (alpha == 0) Alpha_TS[1] - alpha else (Alpha_TS[1] - alpha) / alpha
  CP_TS <- proportion(out$result_TS$survival$Alpha[, 1], out$result_TS$survival$Alpha[, 2], alpha)
  
  # --- Results ---
  result_table <- rbind(
    JM = c(bias_alphJM, esd_alphJM, asd_alphJM, CP_JM),
    TS = c(bias_alphTS, esd_alphTS, asd_alphTS, CP_TS)
  )
  colnames(result_table) <- c("Bias", "ESD", "ASD", "CP")
  summary_tables[[dataset_name]] <- result_table
}

# --- Combine all results ---
final_table <- do.call(rbind, lapply(names(summary_tables), function(name) {
  cbind(Dataset = name, Method = rownames(summary_tables[[name]]), summary_tables[[name]])
}))

final_table_df <- as.data.frame(final_table, stringsAsFactors = FALSE)
colnames(final_table_df) <- c("Dataset", "Method", "Bias", "ESD", "ASD", "CP")

final_table_clean <- final_table_df %>%
  mutate(
    Bias = round(as.numeric(Bias), 2),
    ESD = round(as.numeric(ESD), 3),
    ASD = round(as.numeric(ASD), 3),
    CP = round(as.numeric(CP), 1)
  )

# --- Table ---
xtable(final_table_clean, include.rownames = FALSE)

                     
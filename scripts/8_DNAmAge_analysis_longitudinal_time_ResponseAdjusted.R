

################################################################################
# SETUP
################################################################################

work_dir <- "S://Project/WP-epigenetics/06_DNAmAge/"


################################################################################
# IMPORT / LOAD
################################################################################

# SET WORKING DIRECTORY
setwd(work_dir)

# LIBRARIES
library(pkgload, lib.loc = "C:/Users/kihof7027/Downloads/")
library(DBI, lib.loc = "Z:/Bioconductorpackages/319")
library(tidyselect, lib.loc = "Z:/Bioconductorpackages/319")
library(shiny, lib.loc = "Z:/Bioconductorpackages/319")
library(readxl)
library(ggplot2)
library(lme4)
library(lmerTest)
library(writexl)
library(tidyverse)
library(robustlmm)

# LOAD SAMPLESHEET
samplesheet <- as.data.frame(read_xlsx("output/Tables/Samplesheet_OCD.xlsx"))

# LOAD THEME
load("Resources/theme.RData")
load("Resources/theme_transparent.Rdata")

# LOAD DNAmAge
DNAmAge <- as.data.frame(read_xlsx("output/Tables/DNAmAge_All.xlsx"))

# nr of independent tests
load("output/RData/Meff_OCD.RData")


################################################################################
# PREPARE
################################################################################

clocks <- c("PC_Horvath", "PC_GrimAge", "PC_PhenoAge", "DunedinPACE")
DNAmAge <- DNAmAge[which(DNAmAge$Clock %in% clocks), ]


# SAMPLES: At least 2 time points
time_per_indiv <- data.frame(indiv = unique(samplesheet$Resp_nr), D1 = NA, D4 = NA, M3 = NA)

for (i in 1:nrow(time_per_indiv)){
  RespNr <- time_per_indiv$indiv[i]
  timepoints <- samplesheet$Time_point[which(samplesheet$Resp_nr == RespNr)]
  
  if("D1" %in% timepoints){
    time_per_indiv$D1[i] <- "Y"
  }
  if("D4" %in% timepoints){
    time_per_indiv$D4[i] <- "Y"
  }
  if("M3" %in% timepoints){
    time_per_indiv$M3[i] <- "Y"
  }
}


time_per_indiv$Status <- apply(time_per_indiv, 1, function(row){
  if(sum(row == "Y", na.rm = TRUE) >= 2){
    return("include")
  } else {
    return(NA)
  }
})

included_samples <- time_per_indiv$indiv[which(time_per_indiv$Status == "include")]

samples_long <- samplesheet[which(samplesheet$Resp_nr %in% included_samples), ]
table(samples_long$Time_point)
# D1  D4  M3 
# 413 405 401 
length(unique(samples_long$Resp_nr)) #419


# VARIABELS
Basename <- samples_long$Basename
Sex <- as.factor(samples_long$Sex)
Age <- scale(samples_long$Age)
Smoking <- scale(samples_long$Smokingresid)
BMI <- scale(samples_long$BMIresid)
Time_point <- as.factor(samples_long$Time_point)
Resp_nr <- as.factor(samples_long$Resp_nr)
Response <- scale(samples_long$Percent_improvement)


# make comorbiditiy df
comorbidity_df <- data.frame(
  depression = samples_long$Depression_combined_CURRENT,
  bipolar = samples_long$Bipolar_combined_CURRENT,
  panic = samples_long$PanicDisorder_CURRENTMONTH,
  agora = samples_long$Agoraphobia_CURRENT,
  socialphobia = samples_long$SocialPhobia_combined_CURRENTMONTH,
  PTSD = samples_long$PTSD_CURRENTMONTH,
  anorexia = samples_long$Anorexia_CURRENT3MONTHS,
  bulimia = samples_long$Bulimia_CURRENT3MONTHS,
  GAD = samples_long$GAD_CURRENT6MONTHS,
  burden = NA,
  Basename = samples_long$Basename)

# burden of comorbidities
for (i in c(1:nrow(comorbidity_df))){
  comorbidities <- unlist(as.vector(comorbidity_df[i, ]))
  length_Y <- sum(comorbidities == "Y", na.rm = TRUE)
  comorbidity_df$burden[i] <- length_Y
}


medication_df <- data.frame(
  Psychoactive_medicine = samples_long$Psychoactive_medicine,
  Antidepressants = samples_long$Antidepressants,
  Stimulants = samples_long$Stimulants,
  Antipsychotics = samples_long$Antipsychotics,
  Mood_stabilizers = samples_long$Mood_stabilizers,
  Sedatives = samples_long$Sedatives,
  Anxiolytics = samples_long$Anxiolytics,
  Basename = samples_long$Basename)


################################################################################
# Time only
################################################################################

results_TP <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 6))
colnames(results_TP) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  mAge <- DNAmAge[which(DNAmAge$Sample %in% samples_long$Basename & DNAmAge$Clock == clock), "AgeAccelResid"]

  model <- rlmer(formula = mAge ~ Time_point + Age + Sex + Response + (1|Resp_nr))
  sum <- summary(model)
  
  #compute p values (robust Wald test)
  coefs <- summary(model)$coefficients
  test_stat <- coefs[, "Estimate"] / coefs[, "Std. Error"]
  p_values <- 2 * pnorm(abs(test_stat), lower.tail = FALSE)
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(results) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- "Post"
  results[1, 3] <- sum$coefficients["Time_pointD4", "Estimate"] 
  results[1, 4] <- sum$coefficients["Time_pointD4", "Std. Error"]
  results[1, 5] <- sum$coefficients["Time_pointD4", "t value"]
  results[1, 6] <- p_values[names(p_values) == "Time_pointD4"]
  
  results_TP <- rbind(results_TP, results)
  
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(results) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- "Follow-up"
  results[1, 3] <- sum$coefficients["Time_pointM3", "Estimate"]
  results[1, 4] <- sum$coefficients["Time_pointM3", "Std. Error"]
  results[1, 5] <- sum$coefficients["Time_pointM3", "t value"]
  results[1, 6] <- p_values[names(p_values) == "Time_pointM3"]
  
  results_TP <- rbind(results_TP, results)
}


# adjust p values
results_TP$adj_p <- results_TP$p_value * Meff
results_TP$adj_p <- ifelse(results_TP$adj_p > 1, 1, results_TP$adj_p)
write_xlsx(results_TP, path = "output/Tables/Results_OCD_TimeOnly_ResponseAdjusted.xlsx")


results_TP$Time_point <- factor(results_TP$Time_point, levels = c("Post", "Follow-up"))
results_TP$fill_color <- ifelse(results_TP$adj_p <= 0.1, "#d4f3fc", "white")
results_TP$fill_color <- ifelse(results_TP$adj_p <= 0.05, "#18a6b9", results_TP$fill_color)
results_TP$rounded_p <- format(results_TP$adj_p, scientific = TRUE, digits = 3)
results_TP$clock <- factor(results_TP$clock, levels = sort(unique(results_TP$clock), decreasing = TRUE))

TP_plot <- ggplot(results_TP, aes(x=Time_point, y=clock, fill=fill_color))+
  geom_tile(color = "black") +
  geom_text(aes(label = rounded_p), color = "black", size = 3) +
  scale_fill_identity() +
  th + th_transparent + 
  theme(
    axis.text.x =  element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "none",
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("AgeAccelResid ~ Time\n(adj. p values)")
TP_plot

ggsave("output/Figures/Results/Results_OCD_TPonly_ResponseAdjusted.pdf", TP_plot, units = "cm", width = 8, height = 8, dpi = 400)
ggsave("output/Figures/Results/Results_OCD_TPonly_ResponseAdjusted.svg", TP_plot, device = "svg", units = "cm", width = 8, height = 8, dpi = 400)
save(TP_plot, file = "output/Figures/Results/Results_OCD_TPonly_ResponseAdjusted.RData")


################################################################################
# Time only (INCL SMOKING/BMI)
################################################################################

results_TP_BMIsmoking <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 6))
colnames(results_TP_BMIsmoking) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  mAge <- DNAmAge[which(DNAmAge$Sample %in% samples_long$Basename & DNAmAge$Clock == clock), "AgeAccelResid"]
  
  model <- rlmer(formula = mAge ~ Time_point + Age + Sex + BMI + Smoking + Response + (1|Resp_nr))
  sum <- summary(model)
  
  #compute p values (robust Wald test)
  coefs <- summary(model)$coefficients
  test_stat <- coefs[, "Estimate"] / coefs[, "Std. Error"]
  p_values <- 2 * pnorm(abs(test_stat), lower.tail = FALSE)
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(results) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- "Post"
  results[1, 3] <- sum$coefficients["Time_pointD4", "Estimate"] 
  results[1, 4] <- sum$coefficients["Time_pointD4", "Std. Error"]
  results[1, 5] <- sum$coefficients["Time_pointD4", "t value"]
  results[1, 6] <- p_values[names(p_values) == "Time_pointD4"]
  
  results_TP_BMIsmoking <- rbind(results_TP_BMIsmoking, results)
  
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(results) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- "Follow-up"
  results[1, 3] <- sum$coefficients["Time_pointM3", "Estimate"]
  results[1, 4] <- sum$coefficients["Time_pointM3", "Std. Error"]
  results[1, 5] <- sum$coefficients["Time_pointM3", "t value"]
  results[1, 6] <- p_values[names(p_values) == "Time_pointM3"]
  
  results_TP_BMIsmoking <- rbind(results_TP_BMIsmoking, results)
}


# adjust p values
results_TP_BMIsmoking$adj_p <- results_TP_BMIsmoking$p_value * Meff
results_TP_BMIsmoking$adj_p <- ifelse(results_TP_BMIsmoking$adj_p > 1, 1, results_TP_BMIsmoking$adj_p)
write_xlsx(results_TP_BMIsmoking, path = "output/Tables/Results_OCD_TimeOnly_BMIsmoking_ResponseAdjusted.xlsx")


results_TP_BMIsmoking$Time_point <- factor(results_TP_BMIsmoking$Time_point, levels = c("Post", "Follow-up"))
results_TP_BMIsmoking$fill_color <- ifelse(results_TP_BMIsmoking$adj_p <= 0.1, "#d4f3fc", "white")
results_TP_BMIsmoking$fill_color <- ifelse(results_TP_BMIsmoking$adj_p <= 0.05, "#18a6b9", results_TP_BMIsmoking$fill_color)
results_TP_BMIsmoking$rounded_p <- round(results_TP_BMIsmoking$adj_p, 3)
results_TP_BMIsmoking$clock <- factor(results_TP_BMIsmoking$clock, levels = sort(unique(results_TP_BMIsmoking$clock), decreasing = TRUE))

TP_plot_Smoking_BMI <- ggplot(results_TP_BMIsmoking, aes(x=Time_point, y=clock, fill=fill_color))+
  geom_tile(color = "black") +
  geom_text(aes(label = rounded_p), color = "black", size = 3) +
  scale_fill_identity() +
  th + th_transparent + 
  theme(
    axis.text.x =  element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "none",
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("AgeAccelResid ~ Time\n(adj. p values)")
TP_plot_Smoking_BMI

ggsave("output/Figures/Results/Results_OCD_TP_BMIsmoking_ResponseAdjusted.pdf", TP_plot_Smoking_BMI, units = "cm", width = 8, height = 8, dpi = 400)
ggsave("output/Figures/Results/Results_OCD_TP_BMIsmoking_ResponseAdjusted.svg", TP_plot_Smoking_BMI, device = "svg", units = "cm", width = 8, height = 8, dpi = 400)
save(TP_plot_Smoking_BMI, file = "output/Figures/Results/Results_OCD_TP_BMIsmoking_ResponseAdjusted.RData")



################################################################################
# Time only (INCL SMOKING ONLY)
################################################################################

results_TP_SmokingOnly <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 6))
colnames(results_TP_SmokingOnly) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  mAge <- DNAmAge[which(DNAmAge$Sample %in% samples_long$Basename & DNAmAge$Clock == clock), "AgeAccelResid"]
  
  model <- rlmer(formula = mAge ~ Time_point + Age + Sex + Smoking + Response + (1|Resp_nr))
  sum <- summary(model)
  
  #compute p values (robust Wald test)
  coefs <- summary(model)$coefficients
  test_stat <- coefs[, "Estimate"] / coefs[, "Std. Error"]
  p_values <- 2 * pnorm(abs(test_stat), lower.tail = FALSE)
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(results) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- "Post"
  results[1, 3] <- sum$coefficients["Time_pointD4", "Estimate"] 
  results[1, 4] <- sum$coefficients["Time_pointD4", "Std. Error"]
  results[1, 5] <- sum$coefficients["Time_pointD4", "t value"]
  results[1, 6] <- p_values[names(p_values) == "Time_pointD4"]
  
  results_TP_SmokingOnly <- rbind(results_TP_SmokingOnly, results)
  
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(results) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- "Follow-up"
  results[1, 3] <- sum$coefficients["Time_pointM3", "Estimate"]
  results[1, 4] <- sum$coefficients["Time_pointM3", "Std. Error"]
  results[1, 5] <- sum$coefficients["Time_pointM3", "t value"]
  results[1, 6] <- p_values[names(p_values) == "Time_pointM3"]
  
  results_TP_SmokingOnly <- rbind(results_TP_SmokingOnly, results)
}


# adjust p values
results_TP_SmokingOnly$adj_p <- results_TP_SmokingOnly$p_value * Meff
results_TP_SmokingOnly$adj_p <- ifelse(results_TP_SmokingOnly$adj_p > 1, 1, results_TP_SmokingOnly$adj_p)
write_xlsx(results_TP_SmokingOnly, path = "output/Tables/Results_OCD_TimeOnly_SmokingOnly_ResponseAdjusted.xlsx")


results_TP_SmokingOnly$Time_point <- factor(results_TP_SmokingOnly$Time_point, levels = c("Post", "Follow-up"))
results_TP_SmokingOnly$fill_color <- ifelse(results_TP_SmokingOnly$adj_p <= 0.1, "#d4f3fc", "white")
results_TP_SmokingOnly$fill_color <- ifelse(results_TP_SmokingOnly$adj_p <= 0.05, "#18a6b9", results_TP_SmokingOnly$fill_color)
results_TP_SmokingOnly$rounded_p <- round(results_TP_SmokingOnly$adj_p, 3)
results_TP_SmokingOnly$clock <- factor(results_TP_SmokingOnly$clock, levels = sort(unique(results_TP_SmokingOnly$clock), decreasing = TRUE))

TP_plot_SmokingOnly <- ggplot(results_TP_SmokingOnly, aes(x=Time_point, y=clock, fill=fill_color))+
  geom_tile(color = "black") +
  geom_text(aes(label = rounded_p), color = "black", size = 3) +
  scale_fill_identity() +
  th + th_transparent + 
  theme(
    axis.text.x =  element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "none",
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("AgeAccelResid ~ Time\n(adj. p values)")
TP_plot_SmokingOnly

ggsave("output/Figures/Results/Results_OCD_TP_SmokingOnly_ResponseAdjusted.pdf", TP_plot_SmokingOnly, units = "cm", width = 8, height = 8, dpi = 400)
ggsave("output/Figures/Results/Results_OCD_TP_SmokingOnly_ResponseAdjusted.svg", TP_plot_SmokingOnly, device = "svg", units = "cm", width = 8, height = 8, dpi = 400)
save(TP_plot_SmokingOnly, file = "output/Figures/Results/Results_OCD_TP_SmokingOnly_ResponseAdjusted.RData")




################################################################################
# Time only (INCL ONLY BMI)
################################################################################

results_TP_BMIOnly <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 6))
colnames(results_TP_BMIOnly) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  mAge <- DNAmAge[which(DNAmAge$Sample %in% samples_long$Basename & DNAmAge$Clock == clock), "AgeAccelResid"]
  
  model <- rlmer(formula = mAge ~ Time_point + Age + Sex + BMI + Response + (1|Resp_nr))
  sum <- summary(model)
  
  #compute p values (robust Wald test)
  coefs <- summary(model)$coefficients
  test_stat <- coefs[, "Estimate"] / coefs[, "Std. Error"]
  p_values <- 2 * pnorm(abs(test_stat), lower.tail = FALSE)
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(results) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- "Post"
  results[1, 3] <- sum$coefficients["Time_pointD4", "Estimate"] 
  results[1, 4] <- sum$coefficients["Time_pointD4", "Std. Error"]
  results[1, 5] <- sum$coefficients["Time_pointD4", "t value"]
  results[1, 6] <- p_values[names(p_values) == "Time_pointD4"]
  
  results_TP_BMIOnly <- rbind(results_TP_BMIOnly, results)
  
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(results) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- "Follow-up"
  results[1, 3] <- sum$coefficients["Time_pointM3", "Estimate"]
  results[1, 4] <- sum$coefficients["Time_pointM3", "Std. Error"]
  results[1, 5] <- sum$coefficients["Time_pointM3", "t value"]
  results[1, 6] <- p_values[names(p_values) == "Time_pointM3"]
  
  results_TP_BMIOnly <- rbind(results_TP_BMIOnly, results)
}


# adjust p values
results_TP_BMIOnly$adj_p <- results_TP_BMIOnly$p_value * Meff
results_TP_BMIOnly$adj_p <- ifelse(results_TP_BMIOnly$adj_p > 1, 1, results_TP_BMIOnly$adj_p)
write_xlsx(results_TP_BMIOnly, path = "output/Tables/Results_OCD_TimeOnly_BMIOnly_ResponseAdjusted.xlsx")


results_TP_BMIOnly$Time_point <- factor(results_TP_BMIOnly$Time_point, levels = c("Post", "Follow-up"))
results_TP_BMIOnly$fill_color <- ifelse(results_TP_BMIOnly$adj_p <= 0.1, "#d4f3fc", "white")
results_TP_BMIOnly$fill_color <- ifelse(results_TP_BMIOnly$adj_p <= 0.05, "#18a6b9", results_TP_BMIOnly$fill_color)
results_TP_BMIOnly$rounded_p <- round(results_TP_BMIOnly$adj_p, 3)
results_TP_BMIOnly$clock <- factor(results_TP_BMIOnly$clock, levels = sort(unique(results_TP_BMIOnly$clock), decreasing = TRUE))

TP_plot_BMIOnly <- ggplot(results_TP_BMIOnly, aes(x=Time_point, y=clock, fill=fill_color))+
  geom_tile(color = "black") +
  geom_text(aes(label = rounded_p), color = "black", size = 3) +
  scale_fill_identity() +
  th + th_transparent + 
  theme(
    axis.text.x =  element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "none",
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("AgeAccelResid ~ Time\n(adj. p values)")
TP_plot_BMIOnly

ggsave("output/Figures/Results/Results_OCD_TP_BMIOnly_ResponseAdjusted.pdf", TP_plot_BMIOnly, units = "cm", width = 8, height = 8, dpi = 400)
ggsave("output/Figures/Results/Results_OCD_TP_BMIOnly_ResponseAdjusted.svg", TP_plot_BMIOnly, device = "svg", units = "cm", width = 8, height = 8, dpi = 400)
save(TP_plot_BMIOnly, file = "output/Figures/Results/Results_OCD_TP_BMIOnly_ResponseAdjusted.RData")





################################################################################
# Interaction with comorbidities
################################################################################

results_TP_comorb <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 7))
colnames(results_TP_comorb) <- c("clock", "Comorbidity", "Time_point", "Estimate", "SE", "t_value", "p_value")

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  mAge <- DNAmAge[which(DNAmAge$Sample %in% samples_long$Basename & DNAmAge$Clock == clock), "AgeAccelResid"]
  
  # BURDEN
  comorb = "burden"
  comorb_data <- comorbidity_df$burden
  
  model <- rlmer(mAge ~ Time_point * comorb_data + Age + Sex + Response + (1|Resp_nr))
  sum <- summary(model)
  
  #compute p values (robust Wald test)
  coefs <- summary(model)$coefficients
  test_stat <- coefs[, "Estimate"] / coefs[, "Std. Error"]
  p_values <- 2 * pnorm(abs(test_stat), lower.tail = FALSE)
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 7))
  colnames(results) <- c("clock", "Comorbidity", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- comorb
  results[1, 3] <- "Post"
  results[1, 4] <- sum$coefficients["Time_pointD4:comorb_data", "Estimate"] 
  results[1, 5] <- sum$coefficients["Time_pointD4:comorb_data", "Std. Error"]
  results[1, 6] <- sum$coefficients["Time_pointD4:comorb_data", "t value"]
  results[1, 7] <- p_values[names(p_values) == "Time_pointD4:comorb_data"]
  
  results_TP_comorb <- rbind(results_TP_comorb, results)
  
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 7))
  colnames(results) <- c("clock", "Comorbidity", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- comorb
  results[1, 3] <- "Follow-up"
  results[1, 4] <- sum$coefficients["Time_pointM3:comorb_data", "Estimate"]
  results[1, 5] <- sum$coefficients["Time_pointM3:comorb_data", "Std. Error"]
  results[1, 6] <- sum$coefficients["Time_pointM3:comorb_data", "t value"]
  results[1, 7] <- p_values[names(p_values) == "Time_pointM3:comorb_data"]
  
  results_TP_comorb <- rbind(results_TP_comorb, results)
  
  
  # SINGLE COMORBIDITIES
  for (k in 1:9){
    comorb <- colnames(comorbidity_df)[k]
    comorb_data <- comorbidity_df[, k]
    comorb_data_table <- as.data.frame(table(comorb_data))
    
    if(nrow(comorb_data_table) > 1){
      
      model <- rlmer(mAge ~ Time_point * comorb_data + Age + Sex + Response + (1|Resp_nr))
      sum <- summary(model)
      
      #compute p values (robust Wald test)
      coefs <- summary(model)$coefficients
      test_stat <- coefs[, "Estimate"] / coefs[, "Std. Error"]
      p_values <- 2 * pnorm(abs(test_stat), lower.tail = FALSE)
      
      results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 7))
      colnames(results) <- c("clock", "Comorbidity", "Time_point", "Estimate", "SE", "t_value", "p_value")
      
      results[1, 1] <- clock
      results[1, 2] <- comorb
      results[1, 3] <- "Post"
      results[1, 4] <- sum$coefficients["Time_pointD4:comorb_dataY", "Estimate"] 
      results[1, 5] <- sum$coefficients["Time_pointD4:comorb_dataY", "Std. Error"]
      results[1, 6] <- sum$coefficients["Time_pointD4:comorb_dataY", "t value"]
      results[1, 7] <- p_values[names(p_values) == "Time_pointD4:comorb_dataY"]
      
      results_TP_comorb <- rbind(results_TP_comorb, results)
      
      
      results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 7))
      colnames(results) <- c("clock", "Comorbidity", "Time_point", "Estimate", "SE", "t_value", "p_value")
      
      results[1, 1] <- clock
      results[1, 2] <- comorb
      results[1, 3] <- "Follow-up"
      results[1, 4] <- sum$coefficients["Time_pointM3:comorb_dataY", "Estimate"]
      results[1, 5] <- sum$coefficients["Time_pointM3:comorb_dataY", "Std. Error"]
      results[1, 6] <- sum$coefficients["Time_pointM3:comorb_dataY", "t value"]
      results[1, 7] <- p_values[names(p_values) == "Time_pointM3:comorb_dataY"]
      
      results_TP_comorb <- rbind(results_TP_comorb, results)
    }
  }
}


results_TP_sign <- results_TP[which(results_TP$adj_p <= 0.05), ]
results_TP_sign$TPClock <- paste0(results_TP_sign$Time_point, "_", results_TP_sign$clock)

results_TP_comorb$TPClock <- paste0(results_TP_comorb$Time_point, "_", results_TP_comorb$clock)
results_TP_comorb_selected <- results_TP_comorb[which(results_TP_comorb$TPClock %in% results_TP_sign$TPClock), ]

length(unique(results_TP_comorb_selected$Comorbidity)) #9
results_TP_comorb_selected$adj_p <- results_TP_comorb_selected$p_value * 9
results_TP_comorb_selected$adj_p <- ifelse(results_TP_comorb_selected$adj_p > 1, 1, results_TP_comorb_selected$adj_p)

write_xlsx(results_TP_comorb_selected, path = "output/Tables/Results_OCD_TimeOnly_comorb_ResponseAdjusted.xlsx")


################################################################################
# Interaction with medication
################################################################################


results_TP_medic <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 7))
colnames(results_TP_medic) <- c("clock", "medication", "Time_point", "Estimate", "SE", "t_value", "p_value")

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  mAge <- DNAmAge[which(DNAmAge$Sample %in% samples_long$Basename & DNAmAge$Clock == clock), "AgeAccelResid"]

  # SINGLE medication
  for (k in 1:7){
    medic <- colnames(medication_df)[k]
    medic_data <- medication_df[, k]
    medic_data_table <- as.data.frame(table(medic_data))
    
    if(nrow(medic_data_table) > 1){
      
      model <- rlmer(mAge ~ Time_point * medic_data + Age + Sex + Response + (1|Resp_nr))
      sum <- summary(model)
      
      #compute p values (robust Wald test)
      coefs <- summary(model)$coefficients
      test_stat <- coefs[, "Estimate"] / coefs[, "Std. Error"]
      p_values <- 2 * pnorm(abs(test_stat), lower.tail = FALSE)
      
      results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 7))
      colnames(results) <- c("clock", "medication", "Time_point", "Estimate", "SE", "t_value", "p_value")
      
      results[1, 1] <- clock
      results[1, 2] <- medic
      results[1, 3] <- "Post"
      results[1, 4] <- sum$coefficients["Time_pointD4:medic_dataY", "Estimate"] 
      results[1, 5] <- sum$coefficients["Time_pointD4:medic_dataY", "Std. Error"]
      results[1, 6] <- sum$coefficients["Time_pointD4:medic_dataY", "t value"]
      results[1, 7] <- p_values[names(p_values) == "Time_pointD4:medic_dataY"]
      
      results_TP_medic <- rbind(results_TP_medic, results)
      
      
      results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 7))
      colnames(results) <- c("clock", "medication", "Time_point", "Estimate", "SE", "t_value", "p_value")
      
      results[1, 1] <- clock
      results[1, 2] <- medic
      results[1, 3] <- "Follow-up"
      results[1, 4] <- sum$coefficients["Time_pointM3:medic_dataY", "Estimate"]
      results[1, 5] <- sum$coefficients["Time_pointM3:medic_dataY", "Std. Error"]
      results[1, 6] <- sum$coefficients["Time_pointM3:medic_dataY", "t value"]
      results[1, 7] <- p_values[names(p_values) == "Time_pointM3:medic_dataY"]
      
      results_TP_medic <- rbind(results_TP_medic, results)
    }
  }
}


results_TP_medic$TPClock <- paste0(results_TP_medic$Time_point, "_", results_TP_medic$clock)
results_TP_medic_selected <- results_TP_medic[which(results_TP_medic$TPClock %in% results_TP_sign$TPClock), ]

length(unique(results_TP_medic_selected$medication)) #7
results_TP_medic_selected$adj_p <- results_TP_medic_selected$p_value * 7
results_TP_medic_selected$adj_p <- ifelse(results_TP_medic_selected$adj_p > 1, 1, results_TP_medic_selected$adj_p)
write_xlsx(results_TP_medic_selected, path = "output/Tables/Results_OCD_TimeOnly_medic_ResponseAdjusted.xlsx")



################################################################################
# SENSITIVITY: NON-PC CLOCKS
################################################################################

# LOAD DNAmAge
DNAmAge <- as.data.frame(read_xlsx("output/Tables/DNAmAge_All.xlsx"))
clocks <- c("GrimAge2", "PhenoAge")
DNAmAge <- DNAmAge[which(DNAmAge$Clock %in% clocks), ]

# BASIC MODEL
results_TP_Sensitivity <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 6))
colnames(results_TP_Sensitivity) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  mAge <- DNAmAge[which(DNAmAge$Sample %in% samples_long$Basename & DNAmAge$Clock == clock), "AgeAccelResid"]
  
  model <- rlmer(formula = mAge ~ Time_point + Age + Sex + Response + (1|Resp_nr))
  sum <- summary(model)
  
  #compute p values (robust Wald test)
  coefs <- summary(model)$coefficients
  test_stat <- coefs[, "Estimate"] / coefs[, "Std. Error"]
  p_values <- 2 * pnorm(abs(test_stat), lower.tail = FALSE)
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(results) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- "Post"
  results[1, 3] <- sum$coefficients["Time_pointD4", "Estimate"] 
  results[1, 4] <- sum$coefficients["Time_pointD4", "Std. Error"]
  results[1, 5] <- sum$coefficients["Time_pointD4", "t value"]
  results[1, 6] <- p_values[names(p_values) == "Time_pointD4"]
  
  results_TP_Sensitivity <- rbind(results_TP_Sensitivity, results)
  
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(results) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- "Follow-up"
  results[1, 3] <- sum$coefficients["Time_pointM3", "Estimate"]
  results[1, 4] <- sum$coefficients["Time_pointM3", "Std. Error"]
  results[1, 5] <- sum$coefficients["Time_pointM3", "t value"]
  results[1, 6] <- p_values[names(p_values) == "Time_pointM3"]
  
  results_TP_Sensitivity <- rbind(results_TP_Sensitivity, results)
}

write_xlsx(results_TP_Sensitivity, path = "output/Tables/Sensitivity_OCD_TimeOnly_ResponseAdjusted.xlsx")

# BMI/Smoking Adjusted

results_TP_Sensitivity_BMIsmoking <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 6))
colnames(results_TP_Sensitivity_BMIsmoking) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  mAge <- DNAmAge[which(DNAmAge$Sample %in% samples_long$Basename & DNAmAge$Clock == clock), "AgeAccelResid"]
  
  model <- rlmer(formula = mAge ~ Time_point + Age + Sex + BMI + Smoking + Response + (1|Resp_nr))
  sum <- summary(model)
  
  #compute p values (robust Wald test)
  coefs <- summary(model)$coefficients
  test_stat <- coefs[, "Estimate"] / coefs[, "Std. Error"]
  p_values <- 2 * pnorm(abs(test_stat), lower.tail = FALSE)
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(results) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- "Post"
  results[1, 3] <- sum$coefficients["Time_pointD4", "Estimate"] 
  results[1, 4] <- sum$coefficients["Time_pointD4", "Std. Error"]
  results[1, 5] <- sum$coefficients["Time_pointD4", "t value"]
  results[1, 6] <- p_values[names(p_values) == "Time_pointD4"]
  
  results_TP_Sensitivity_BMIsmoking <- rbind(results_TP_Sensitivity_BMIsmoking, results)
  
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(results) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- "Follow-up"
  results[1, 3] <- sum$coefficients["Time_pointM3", "Estimate"]
  results[1, 4] <- sum$coefficients["Time_pointM3", "Std. Error"]
  results[1, 5] <- sum$coefficients["Time_pointM3", "t value"]
  results[1, 6] <- p_values[names(p_values) == "Time_pointM3"]
  
  results_TP_Sensitivity_BMIsmoking <- rbind(results_TP_Sensitivity_BMIsmoking, results)
}


# adjust p values
write_xlsx(results_TP_Sensitivity_BMIsmoking, path = "output/Tables/Sensitivity_OCD_TimeOnly_ResponseAdjusted_BMIsmoking.xlsx")


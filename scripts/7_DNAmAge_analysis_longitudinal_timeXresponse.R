

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

samples_long <- samplesheet[which(samplesheet$Resp_nr %in% included_samples & !is.na(samplesheet$Percent_improvement)), ]
table(samples_long$Time_point)
#D1  D4  M3 
#395 388 383 
length(unique(samples_long$Resp_nr)) #401


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
# Time x Response
################################################################################

results_TPXResp <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 6))
colnames(results_TPXResp) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  mAge <- DNAmAge[which(DNAmAge$Sample %in% samples_long$Basename & DNAmAge$Clock == clock), "AgeAccelResid"]

  model <- rlmer(formula = mAge ~ Response * Time_point + Age + Sex + (1|Resp_nr))
  sum <- summary(model)
  
  #compute p values (robust Wald test)
  coefs <- summary(model)$coefficients
  test_stat <- coefs[, "Estimate"] / coefs[, "Std. Error"]
  p_values <- 2 * pnorm(abs(test_stat), lower.tail = FALSE)
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(results) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- "Post"
  results[1, 3] <- sum$coefficients["Response:Time_pointD4", "Estimate"] 
  results[1, 4] <- sum$coefficients["Response:Time_pointD4", "Std. Error"]
  results[1, 5] <- sum$coefficients["Response:Time_pointD4", "t value"]
  results[1, 6] <- p_values[names(p_values) == "Response:Time_pointD4"]
  
  results_TPXResp <- rbind(results_TPXResp, results)
  
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(results) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- "Follow-up"
  results[1, 3] <- sum$coefficients["Response:Time_pointM3", "Estimate"]
  results[1, 4] <- sum$coefficients["Response:Time_pointM3", "Std. Error"]
  results[1, 5] <- sum$coefficients["Response:Time_pointM3", "t value"]
  results[1, 6] <- p_values[names(p_values) == "Response:Time_pointM3"]
  
  results_TPXResp <- rbind(results_TPXResp, results)
}


# adjust p values
results_TPXResp$adj_p <- results_TPXResp$p_value * Meff
results_TPXResp$adj_p <- ifelse(results_TPXResp$adj_p > 1, 1, results_TPXResp$adj_p)
write_xlsx(results_TPXResp, path = "output/Tables/Results_OCD_TimeXResponse.xlsx")


results_TPXResp$Time_point <- factor(results_TPXResp$Time_point, levels = c("Post", "Follow-up"))
results_TPXResp$fill_color <- ifelse(results_TPXResp$adj_p <= 0.1, "#d4f3fc", "white")
results_TPXResp$fill_color <- ifelse(results_TPXResp$adj_p <= 0.05, "#18a6b9", results_TPXResp$fill_color)
results_TPXResp$rounded_p <- round(results_TPXResp$adj_p, 3)
results_TPXResp$clock <- factor(results_TPXResp$clock, levels = sort(unique(results_TPXResp$clock), decreasing = TRUE))

TPXResp_plot <- ggplot(results_TPXResp, aes(x=Time_point, y=clock, fill=fill_color))+
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
  ggtitle("AgeAccelResid ~ Time x Resp\n(adj. p values)")
TPXResp_plot

ggsave("output/Figures/Results/Results_OCD_TPXResp.pdf", TPXResp_plot, units = "cm", width = 8, height = 8, dpi = 400)
ggsave("output/Figures/Results/Results_OCD_TPXResp.svg", TPXResp_plot, device = "svg", units = "cm", width = 8, height = 8, dpi = 400)
save(TPXResp_plot, file = "output/Figures/Results/Results_OCD_TPXResp.RData")


################################################################################
# Time x Resp (INCL SMOKING/BMI)
################################################################################

results_TPXResp_BMIsmoking <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 6))
colnames(results_TPXResp_BMIsmoking) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  mAge <- DNAmAge[which(DNAmAge$Sample %in% samples_long$Basename & DNAmAge$Clock == clock), "AgeAccelResid"]
  
  model <- rlmer(formula = mAge ~ Response*Time_point + Age + Sex + BMI + Smoking + (1|Resp_nr))
  sum <- summary(model)
  
  #compute p values (robust Wald test)
  coefs <- summary(model)$coefficients
  test_stat <- coefs[, "Estimate"] / coefs[, "Std. Error"]
  p_values <- 2 * pnorm(abs(test_stat), lower.tail = FALSE)
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(results) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- "Post"
  results[1, 3] <- sum$coefficients["Response:Time_pointD4", "Estimate"] 
  results[1, 4] <- sum$coefficients["Response:Time_pointD4", "Std. Error"]
  results[1, 5] <- sum$coefficients["Response:Time_pointD4", "t value"]
  results[1, 6] <- p_values[names(p_values) == "Response:Time_pointD4"]
  
  results_TPXResp_BMIsmoking <- rbind(results_TPXResp_BMIsmoking, results)
  
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(results) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- "Follow-up"
  results[1, 3] <- sum$coefficients["Response:Time_pointM3", "Estimate"]
  results[1, 4] <- sum$coefficients["Response:Time_pointM3", "Std. Error"]
  results[1, 5] <- sum$coefficients["Response:Time_pointM3", "t value"]
  results[1, 6] <- p_values[names(p_values) == "Response:Time_pointM3"]
  
  results_TPXResp_BMIsmoking <- rbind(results_TPXResp_BMIsmoking, results)
}


# adjust p values
results_TPXResp_BMIsmoking$adj_p <- results_TPXResp_BMIsmoking$p_value * Meff
results_TPXResp_BMIsmoking$adj_p <- ifelse(results_TPXResp_BMIsmoking$adj_p > 1, 1, results_TPXResp_BMIsmoking$adj_p)
write_xlsx(results_TPXResp_BMIsmoking, path = "output/Tables/Results_OCD_ResponsexTime_BMIsmoking.xlsx")


results_TPXResp_BMIsmoking$Time_point <- factor(results_TPXResp_BMIsmoking$Time_point, levels = c("Post", "Follow-up"))
results_TPXResp_BMIsmoking$fill_color <- ifelse(results_TPXResp_BMIsmoking$adj_p <= 0.1, "#d4f3fc", "white")
results_TPXResp_BMIsmoking$fill_color <- ifelse(results_TPXResp_BMIsmoking$adj_p <= 0.05, "#18a6b9", results_TPXResp_BMIsmoking$fill_color)
results_TPXResp_BMIsmoking$rounded_p <- round(results_TPXResp_BMIsmoking$adj_p, 3)
results_TPXResp_BMIsmoking$clock <- factor(results_TPXResp_BMIsmoking$clock, levels = sort(unique(results_TPXResp_BMIsmoking$clock), decreasing = TRUE))

TPXResp_plot_Smoking_BMI <- ggplot(results_TPXResp_BMIsmoking, aes(x=Time_point, y=clock, fill=fill_color))+
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
TPXResp_plot_Smoking_BMI

ggsave("output/Figures/Results/Results_OCD_TPXResp_BMIsmoking.pdf", TPXResp_plot_Smoking_BMI, units = "cm", width = 8, height = 8, dpi = 400)
ggsave("output/Figures/Results/Results_OCD_TPXResp_BMIsmoking.svg", TPXResp_plot_Smoking_BMI, device = "svg", units = "cm", width = 8, height = 8, dpi = 400)
save(TPXResp_plot_Smoking_BMI, file = "output/Figures/Results/Results_OCD_TPXResp_BMIsmoking.RData")




################################################################################
# Time x Resp (INCL SMOKING ONLY)
################################################################################

results_TPXResp_SmokingOnly <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 6))
colnames(results_TPXResp_SmokingOnly) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  mAge <- DNAmAge[which(DNAmAge$Sample %in% samples_long$Basename & DNAmAge$Clock == clock), "AgeAccelResid"]
  
  model <- rlmer(formula = mAge ~ Response*Time_point + Age + Sex + Smoking + (1|Resp_nr))
  sum <- summary(model)
  
  #compute p values (robust Wald test)
  coefs <- summary(model)$coefficients
  test_stat <- coefs[, "Estimate"] / coefs[, "Std. Error"]
  p_values <- 2 * pnorm(abs(test_stat), lower.tail = FALSE)
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(results) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- "Post"
  results[1, 3] <- sum$coefficients["Response:Time_pointD4", "Estimate"] 
  results[1, 4] <- sum$coefficients["Response:Time_pointD4", "Std. Error"]
  results[1, 5] <- sum$coefficients["Response:Time_pointD4", "t value"]
  results[1, 6] <- p_values[names(p_values) == "Response:Time_pointD4"]
  
  results_TPXResp_SmokingOnly <- rbind(results_TPXResp_SmokingOnly, results)
  
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(results) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- "Follow-up"
  results[1, 3] <- sum$coefficients["Response:Time_pointM3", "Estimate"]
  results[1, 4] <- sum$coefficients["Response:Time_pointM3", "Std. Error"]
  results[1, 5] <- sum$coefficients["Response:Time_pointM3", "t value"]
  results[1, 6] <- p_values[names(p_values) == "Response:Time_pointM3"]
  
  results_TPXResp_SmokingOnly <- rbind(results_TPXResp_SmokingOnly, results)
}


# adjust p values
results_TPXResp_SmokingOnly$adj_p <- results_TPXResp_SmokingOnly$p_value * Meff
results_TPXResp_SmokingOnly$adj_p <- ifelse(results_TPXResp_SmokingOnly$adj_p > 1, 1, results_TPXResp_SmokingOnly$adj_p)
write_xlsx(results_TPXResp_SmokingOnly, path = "output/Tables/Results_OCD_ResponsexTime_SmokingOnly.xlsx")


results_TPXResp_SmokingOnly$Time_point <- factor(results_TPXResp_SmokingOnly$Time_point, levels = c("Post", "Follow-up"))
results_TPXResp_SmokingOnly$fill_color <- ifelse(results_TPXResp_SmokingOnly$adj_p <= 0.1, "#d4f3fc", "white")
results_TPXResp_SmokingOnly$fill_color <- ifelse(results_TPXResp_SmokingOnly$adj_p <= 0.05, "#18a6b9", results_TPXResp_SmokingOnly$fill_color)
results_TPXResp_SmokingOnly$rounded_p <- round(results_TPXResp_SmokingOnly$adj_p, 3)
results_TPXResp_SmokingOnly$clock <- factor(results_TPXResp_SmokingOnly$clock, levels = sort(unique(results_TPXResp_SmokingOnly$clock), decreasing = TRUE))

TPXResp_plot_SmokingOnly <- ggplot(results_TPXResp_SmokingOnly, aes(x=Time_point, y=clock, fill=fill_color))+
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
TPXResp_plot_SmokingOnly

ggsave("output/Figures/Results/Results_OCD_TPXResp_SmokingOnly.pdf", TPXResp_plot_SmokingOnly, units = "cm", width = 8, height = 8, dpi = 400)
ggsave("output/Figures/Results/Results_OCD_TPXResp_SmokingOnly.svg", TPXResp_plot_SmokingOnly, device = "svg", units = "cm", width = 8, height = 8, dpi = 400)
save(TPXResp_plot_SmokingOnly, file = "output/Figures/Results/Results_OCD_TPXResp_SmokingOnly.RData")





################################################################################
# Time x Resp (INCL BMI ONLY)
################################################################################

results_TPXResp_BMIOnly <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 6))
colnames(results_TPXResp_BMIOnly) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  mAge <- DNAmAge[which(DNAmAge$Sample %in% samples_long$Basename & DNAmAge$Clock == clock), "AgeAccelResid"]
  
  model <- rlmer(formula = mAge ~ Response*Time_point + Age + Sex + BMI + (1|Resp_nr))
  sum <- summary(model)
  
  #compute p values (robust Wald test)
  coefs <- summary(model)$coefficients
  test_stat <- coefs[, "Estimate"] / coefs[, "Std. Error"]
  p_values <- 2 * pnorm(abs(test_stat), lower.tail = FALSE)
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(results) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- "Post"
  results[1, 3] <- sum$coefficients["Response:Time_pointD4", "Estimate"] 
  results[1, 4] <- sum$coefficients["Response:Time_pointD4", "Std. Error"]
  results[1, 5] <- sum$coefficients["Response:Time_pointD4", "t value"]
  results[1, 6] <- p_values[names(p_values) == "Response:Time_pointD4"]
  
  results_TPXResp_BMIOnly <- rbind(results_TPXResp_BMIOnly, results)
  
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(results) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- "Follow-up"
  results[1, 3] <- sum$coefficients["Response:Time_pointM3", "Estimate"]
  results[1, 4] <- sum$coefficients["Response:Time_pointM3", "Std. Error"]
  results[1, 5] <- sum$coefficients["Response:Time_pointM3", "t value"]
  results[1, 6] <- p_values[names(p_values) == "Response:Time_pointM3"]
  
  results_TPXResp_BMIOnly <- rbind(results_TPXResp_BMIOnly, results)
}


# adjust p values
results_TPXResp_BMIOnly$adj_p <- results_TPXResp_BMIOnly$p_value * Meff
results_TPXResp_BMIOnly$adj_p <- ifelse(results_TPXResp_BMIOnly$adj_p > 1, 1, results_TPXResp_BMIOnly$adj_p)
write_xlsx(results_TPXResp_BMIOnly, path = "output/Tables/Results_OCD_ResponsexTime_BMIOnly.xlsx")


results_TPXResp_BMIOnly$Time_point <- factor(results_TPXResp_BMIOnly$Time_point, levels = c("Post", "Follow-up"))
results_TPXResp_BMIOnly$fill_color <- ifelse(results_TPXResp_BMIOnly$adj_p <= 0.1, "#d4f3fc", "white")
results_TPXResp_BMIOnly$fill_color <- ifelse(results_TPXResp_BMIOnly$adj_p <= 0.05, "#18a6b9", results_TPXResp_BMIOnly$fill_color)
results_TPXResp_BMIOnly$rounded_p <- round(results_TPXResp_BMIOnly$adj_p, 3)
results_TPXResp_BMIOnly$clock <- factor(results_TPXResp_BMIOnly$clock, levels = sort(unique(results_TPXResp_BMIOnly$clock), decreasing = TRUE))

TPXResp_plot_BMI <- ggplot(results_TPXResp_BMIOnly, aes(x=Time_point, y=clock, fill=fill_color))+
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
TPXResp_plot_BMI

ggsave("output/Figures/Results/Results_OCD_TPXResp_BMIOnly.pdf", TPXResp_plot_BMI, units = "cm", width = 8, height = 8, dpi = 400)
ggsave("output/Figures/Results/Results_OCD_TPXResp_BMIOnly.svg", TPXResp_plot_BMI, device = "svg", units = "cm", width = 8, height = 8, dpi = 400)
save(TPXResp_plot_BMI, file = "output/Figures/Results/Results_OCD_TPXResp_BMIOnly.RData")




################################################################################
# Interaction with comorbidities
################################################################################

results_TPXResp_comorb <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 7))
colnames(results_TPXResp_comorb) <- c("clock", "Comorbidity", "Time_point", "Estimate", "SE", "t_value", "p_value")

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  mAge <- DNAmAge[which(DNAmAge$Sample %in% samples_long$Basename & DNAmAge$Clock == clock), "AgeAccelResid"]
  
  # BURDEN
  comorb = "burden"
  comorb_data <- comorbidity_df$burden
  
  model <- rlmer(mAge ~ Response * Time_point * comorb_data + Age + Sex + (1|Resp_nr))
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
  results[1, 4] <- sum$coefficients["Response:Time_pointD4:comorb_data", "Estimate"] 
  results[1, 5] <- sum$coefficients["Response:Time_pointD4:comorb_data", "Std. Error"]
  results[1, 6] <- sum$coefficients["Response:Time_pointD4:comorb_data", "t value"]
  results[1, 7] <- p_values[names(p_values) == "Response:Time_pointD4:comorb_data"]
  
  results_TPXResp_comorb <- rbind(results_TPXResp_comorb, results)
  
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 7))
  colnames(results) <- c("clock", "Comorbidity", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- comorb
  results[1, 3] <- "Follow-up"
  results[1, 4] <- sum$coefficients["Response:Time_pointM3:comorb_data", "Estimate"]
  results[1, 5] <- sum$coefficients["Response:Time_pointM3:comorb_data", "Std. Error"]
  results[1, 6] <- sum$coefficients["Response:Time_pointM3:comorb_data", "t value"]
  results[1, 7] <- p_values[names(p_values) == "Response:Time_pointM3:comorb_data"]
  
  results_TPXResp_comorb <- rbind(results_TPXResp_comorb, results)
  
  
  # SINGLE COMORBIDITIES
  for (k in 1:9){
    comorb <- colnames(comorbidity_df)[k]
    comorb_data <- comorbidity_df[, k]
    comorb_data_table <- as.data.frame(table(comorb_data))
    
    if(nrow(comorb_data_table) > 1){
      
      model <- rlmer(mAge ~ Response * Time_point * comorb_data + Age + Sex + (1|Resp_nr))
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
      results[1, 4] <- sum$coefficients["Response:Time_pointD4:comorb_dataY", "Estimate"] 
      results[1, 5] <- sum$coefficients["Response:Time_pointD4:comorb_dataY", "Std. Error"]
      results[1, 6] <- sum$coefficients["Response:Time_pointD4:comorb_dataY", "t value"]
      results[1, 7] <- p_values[names(p_values) == "Response:Time_pointD4:comorb_dataY"]
      
      results_TPXResp_comorb <- rbind(results_TPXResp_comorb, results)
      
      
      results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 7))
      colnames(results) <- c("clock", "Comorbidity", "Time_point", "Estimate", "SE", "t_value", "p_value")
      
      results[1, 1] <- clock
      results[1, 2] <- comorb
      results[1, 3] <- "Follow-up"
      results[1, 4] <- sum$coefficients["Response:Time_pointM3:comorb_dataY", "Estimate"]
      results[1, 5] <- sum$coefficients["Response:Time_pointM3:comorb_dataY", "Std. Error"]
      results[1, 6] <- sum$coefficients["Response:Time_pointM3:comorb_dataY", "t value"]
      results[1, 7] <- p_values[names(p_values) == "Response:Time_pointM3:comorb_dataY"]
      
      results_TPXResp_comorb <- rbind(results_TPXResp_comorb, results)
    }
  }
}


results_TPXResp_sign <- results_TPXResp[which(results_TPXResp$adj_p <= 0.05), ]
results_TPXResp_sign$TPClock <- paste0(results_TPXResp_sign$Time_point, "_", results_TPXResp_sign$clock)

results_TPXResp_comorb$TPClock <- paste0(results_TPXResp_comorb$Time_point, "_", results_TPXResp_comorb$clock)
results_TPXResp_comorb_selected <- results_TPXResp_comorb[which(results_TPXResp_comorb$TPClock %in% results_TPXResp_sign$TPClock), ]

length(unique(results_TPXResp_comorb_selected$Comorbidity)) #9
results_TPXResp_comorb_selected$adj_p <- results_TPXResp_comorb_selected$p_value * 9
results_TPXResp_comorb_selected$adj_p <- ifelse(results_TPXResp_comorb_selected$adj_p > 1, 1, results_TPXResp_comorb_selected$adj_p)

write_xlsx(results_TPXResp_comorb_selected, path = "output/Tables/Results_OCD_ResponsexTime_comorb.xlsx")


################################################################################
# Interaction with medication
################################################################################


results_TPXResp_medic <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 7))
colnames(results_TPXResp_medic) <- c("clock", "medication", "Time_point", "Estimate", "SE", "t_value", "p_value")

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  mAge <- DNAmAge[which(DNAmAge$Sample %in% samples_long$Basename & DNAmAge$Clock == clock), "AgeAccelResid"]
  
  # SINGLE medication
  for (k in 1:7){
    medic <- colnames(medication_df)[k]
    medic_data <- medication_df[, k]
    medic_data_table <- as.data.frame(table(medic_data))
    
    if(nrow(medic_data_table) > 1){
      
      model <- rlmer(mAge ~ Response * Time_point * medic_data + Age + Sex + (1|Resp_nr))
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
      results[1, 4] <- sum$coefficients["Response:Time_pointD4:medic_dataY", "Estimate"] 
      results[1, 5] <- sum$coefficients["Response:Time_pointD4:medic_dataY", "Std. Error"]
      results[1, 6] <- sum$coefficients["Response:Time_pointD4:medic_dataY", "t value"]
      results[1, 7] <- p_values[names(p_values) == "Response:Time_pointD4:medic_dataY"]
      
      results_TPXResp_medic <- rbind(results_TPXResp_medic, results)
      
      
      results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 7))
      colnames(results) <- c("clock", "medication", "Time_point", "Estimate", "SE", "t_value", "p_value")
      
      results[1, 1] <- clock
      results[1, 2] <- medic
      results[1, 3] <- "Follow-up"
      results[1, 4] <- sum$coefficients["Response:Time_pointM3:medic_dataY", "Estimate"]
      results[1, 5] <- sum$coefficients["Response:Time_pointM3:medic_dataY", "Std. Error"]
      results[1, 6] <- sum$coefficients["Response:Time_pointM3:medic_dataY", "t value"]
      results[1, 7] <- p_values[names(p_values) == "Response:Time_pointM3:medic_dataY"]
      
      results_TPXResp_medic <- rbind(results_TPXResp_medic, results)
    }
  }
}


results_TPXResp_medic$TPClock <- paste0(results_TPXResp_medic$Time_point, "_", results_TPXResp_medic$clock)
results_TPXResp_medic_selected <- results_TPXResp_medic[which(results_TPXResp_medic$TPClock %in% results_TPXResp_sign$TPClock), ]

length(unique(results_TPXResp_medic_selected$medication)) #7
results_TPXResp_medic_selected$adj_p <- results_TPXResp_medic_selected$p_value * 7
results_TPXResp_medic_selected$adj_p <- ifelse(results_TPXResp_medic_selected$adj_p > 1, 1, results_TPXResp_medic_selected$adj_p)
write_xlsx(results_TPXResp_medic_selected, path = "output/Tables/Results_OCD_ResponsexTime_medic.xlsx")



################################################################################
# SENSITIVITY: NON-PC CLOCKS
################################################################################

# LOAD DNAmAge
DNAmAge <- as.data.frame(read_xlsx("output/Tables/DNAmAge_All.xlsx"))
clocks <- c("GrimAge2")
DNAmAge <- DNAmAge[which(DNAmAge$Clock %in% clocks), ]

# BASIC MODEL
results_TPXResp_Sensitivity <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 6))
colnames(results_TPXResp_Sensitivity) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  mAge <- DNAmAge[which(DNAmAge$Sample %in% samples_long$Basename & DNAmAge$Clock == clock), "AgeAccelResid"]
  
  model <- rlmer(formula = mAge ~ Response * Time_point + Age + Sex + (1|Resp_nr))
  sum <- summary(model)
  
  #compute p values (robust Wald test)
  coefs <- summary(model)$coefficients
  test_stat <- coefs[, "Estimate"] / coefs[, "Std. Error"]
  p_values <- 2 * pnorm(abs(test_stat), lower.tail = FALSE)
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(results) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- "Post"
  results[1, 3] <- sum$coefficients["Response:Time_pointD4", "Estimate"] 
  results[1, 4] <- sum$coefficients["Response:Time_pointD4", "Std. Error"]
  results[1, 5] <- sum$coefficients["Response:Time_pointD4", "t value"]
  results[1, 6] <- p_values[names(p_values) == "Response:Time_pointD4"]
  
  results_TPXResp_Sensitivity <- rbind(results_TPXResp_Sensitivity, results)
  
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(results) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- "Follow-up"
  results[1, 3] <- sum$coefficients["Response:Time_pointM3", "Estimate"]
  results[1, 4] <- sum$coefficients["Response:Time_pointM3", "Std. Error"]
  results[1, 5] <- sum$coefficients["Response:Time_pointM3", "t value"]
  results[1, 6] <- p_values[names(p_values) == "Response:Time_pointM3"]
  
  results_TPXResp_Sensitivity <- rbind(results_TPXResp_Sensitivity, results)
}

write_xlsx(results_TPXResp_Sensitivity, path = "output/Tables/Sensitivity_OCD_ResponsexTime.xlsx")

#

results_TPXResp_Sensitivity_BMIsmoking <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 6))
colnames(results_TPXResp_Sensitivity_BMIsmoking) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  mAge <- DNAmAge[which(DNAmAge$Sample %in% samples_long$Basename & DNAmAge$Clock == clock), "AgeAccelResid"]
  
  model <- rlmer(formula = mAge ~ Response * Time_point + Age + Sex + BMI + Smoking + (1|Resp_nr))
  sum <- summary(model)
  
  #compute p values (robust Wald test)
  coefs <- summary(model)$coefficients
  test_stat <- coefs[, "Estimate"] / coefs[, "Std. Error"]
  p_values <- 2 * pnorm(abs(test_stat), lower.tail = FALSE)
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(results) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- "Post"
  results[1, 3] <- sum$coefficients["Response:Time_pointD4", "Estimate"] 
  results[1, 4] <- sum$coefficients["Response:Time_pointD4", "Std. Error"]
  results[1, 5] <- sum$coefficients["Response:Time_pointD4", "t value"]
  results[1, 6] <- p_values[names(p_values) == "Response:Time_pointD4"]
  
  results_TPXResp_Sensitivity_BMIsmoking <- rbind(results_TPXResp_Sensitivity_BMIsmoking, results)
  
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(results) <- c("clock", "Time_point", "Estimate", "SE", "t_value", "p_value")
  
  results[1, 1] <- clock
  results[1, 2] <- "Follow-up"
  results[1, 3] <- sum$coefficients["Response:Time_pointM3", "Estimate"]
  results[1, 4] <- sum$coefficients["Response:Time_pointM3", "Std. Error"]
  results[1, 5] <- sum$coefficients["Response:Time_pointM3", "t value"]
  results[1, 6] <- p_values[names(p_values) == "Response:Time_pointM3"]
  
  results_TPXResp_Sensitivity_BMIsmoking <- rbind(results_TPXResp_Sensitivity_BMIsmoking, results)
}


# adjust p values
write_xlsx(results_TPXResp_Sensitivity_BMIsmoking, path = "output/Tables/Sensitivity_OCD_ResponsexTime_BMIsmoking.xlsx")






################################################################################
# Visualisation
################################################################################

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  samples_clin <- samplesheet[which(samplesheet$Basename %in% samples_long$Basename), c("Basename", "Resp_nr", "Response_status")]
  
  DNAmAge_red <- DNAmAge[which(DNAmAge$Clock == clock), ]
  
  merged_DNAmAge <- merge(DNAmAge_red, samples_clin, by.x = "Sample", by.y = "Basename", all.y = TRUE)
  merged_DNAmAge <- merged_DNAmAge[which(!is.na(merged_DNAmAge$Response_status)), ]
  
  
  df_avg <- as.data.frame(merged_DNAmAge %>%
                            group_by(Time_point, Response_status) %>%
                            summarise(
                              avg = mean(AgeAccelResid, na.rm = TRUE),
                              se = sd(AgeAccelResid, na.rm = TRUE) / sqrt(n())))
  
  df_avg$Time_point[which(df_avg$Time_point == "D1")] <- "pre"
  df_avg$Time_point[which(df_avg$Time_point == "D4")] <- "post"
  df_avg$Time_point[which(df_avg$Time_point == "M3")] <- "follow-up"
  
  df_avg$Time_point <- factor(df_avg$Time_point, levels = c("pre", "post", "follow-up"))
  
  Age_plot <- ggplot(df_avg, aes(x = Time_point, y = avg, color = Response_status, group = Response_status)) +
    geom_point(size = 3) +
    geom_line() +
    geom_errorbar(aes(ymin=avg-se, ymax=avg+se), width = 0.2) +
    th + th_transparent +
    labs(
      x = "Time point",
      y = "Age acceleration (residuals)"
    ) + 
    ggtitle(clock) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    theme(legend.position = "top", legend.title = element_blank(), legend.box.background = element_blank(), legend.key = element_blank()) +
    scale_color_manual(values=c("responder"= "#084081", "non-responder" = "#891a1a")) +
    theme(plot.title = element_text(hjust = 0.5))
  Age_plot
  
  ggsave(paste0("output/Figures/ResponseTime/LongitudinalCateg_AgeAccel_ResponseStatus_", clock, "_categ.pdf"), Age_plot, units = "cm", width = 9, height = 10, dpi = 400)
  save(Age_plot, file = paste0("output/Figures/ResponseTime/LongitudinalCateg_AgeAccel_ResponseStatus_", clock, "_categ.RData"))
  
}





for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  samples_clin <- samplesheet[which(samplesheet$Basename %in% samples_long$Basename), c("Basename", "Resp_nr", "Percent_improvement")]
  
  DNAmAge_red <- DNAmAge[which(DNAmAge$Clock == clock), ]
  
  merged_DNAmAge <- merge(DNAmAge_red, samples_clin, by.x = "Sample", by.y = "Basename", all.y = TRUE)
  merged_DNAmAge <- merged_DNAmAge[which(!is.na(merged_DNAmAge$Percent_improvement)), ]
  
  unique_Resp_nr <- unique(merged_DNAmAge$Resp_nr)
  
  plot_df <- data.frame(Resp_nr = unique_Resp_nr, Percent_improvement = NA, FU_Pre = NA)
  
  for (i in 1:nrow(plot_df)){
    Resp_nr <- plot_df$Resp_nr[i]
    plot_df$Percent_improvement[i] <- unique(samplesheet$Percent_improvement[which(samplesheet$Resp_nr == Resp_nr)])
    FU <- merged_DNAmAge$AgeAccelResid[which(merged_DNAmAge$Resp_nr == Resp_nr & merged_DNAmAge$Time_point == "M3")]
    Pre <- merged_DNAmAge$AgeAccelResid[which(merged_DNAmAge$Resp_nr == Resp_nr & merged_DNAmAge$Time_point == "D1")]
    if(length(Pre) > 0 & length(FU) > 0){
      plot_df$FU_Pre[i] <- FU - Pre
    }
  }
  
  plot_df <- na.omit(plot_df)
  
  corr <- round(cor(plot_df$Percent_improvement, plot_df$FU_Pre), 3)
  
  Corr_plot <- ggplot(plot_df, aes(x = Percent_improvement, y = FU_Pre)) +
    geom_point(color = "#084081", alpha = 0.15) + 
    geom_smooth(color = "#084081", method = "lm", linewidth=1.2, alpha = 0.35) +
    scale_color_identity() +
    ggtitle(clock) +
    th + th_transparent +
    xlab("% clinical improvement") + ylab("Change Age Acceleration (FU-Pre)") +
    theme(plot.title = element_text(hjust = 0.5)) + 
    annotate("text", x=15, y = min(plot_df$FU_Pre), label = paste0("r = ", corr), color = "black")
  
  ggsave(paste0("output/Figures/ResponseTime//FollowUpChange_AgeAccel_Response_", clock, ".pdf"), Corr_plot, units = "cm", width = 7, height = 10, dpi = 400)
  save(Corr_plot, file = paste0("output/Figures/ResponseTime/FollowUpChange_AgeAccel_Response_", clock, ".RData"))
}




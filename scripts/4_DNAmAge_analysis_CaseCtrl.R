

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
library(writexl)
library(robustbase)

# LOAD SAMPLESHEET
samplesheet <- as.data.frame(read_xlsx("output/Tables/Samplesheet_CaseCtrl.xlsx"))

# LOAD THEME
load("Resources/theme.RData")
load("Resources/theme_transparent.Rdata")
source("Resources/Function_star_pvalue.R")

# LOAD DNAmAge
DNAmAge <- as.data.frame(read_xlsx("output/Tables/DNAmAge_All.xlsx"))

# nr of independent tests
load("output/RData/Meff_CaseCtrl.RData")



################################################################################
# PREPARE
################################################################################

# SELECT CLOCKS TO TEST
clocks <- c("PC_Horvath", "PC_GrimAge", "PC_PhenoAge", "DunedinPACE")
DNAmAge <- DNAmAge[which(DNAmAge$Clock %in% clocks & DNAmAge$Sample %in% samplesheet$Basename), ]

Diagnosis <- as.factor(samplesheet$Diagnosis)
Sex <- as.factor(samplesheet$Sex)
Age <- scale(samplesheet$Age)
BMI <- scale(samplesheet$BMIresid)
Smoking <- scale(samplesheet$Smokingresid)

################################################################################
# BASIC MODEL
################################################################################

results_CaseCtrl <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 5))
colnames(results_CaseCtrl) <- c("Estimate", "SE", "t_value", "p_value", "clock")

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  AgeAccelResid <- DNAmAge[which(DNAmAge$Clock == clock), "AgeAccelResid"]
  
  model <- lmrob(AgeAccelResid ~ Diagnosis + Age + Sex)
  sum <- summary(model)
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 5))
  colnames(results) <- c("Estimate", "SE", "t_value", "p_value", "clock")
  
  results[1, 1] <- sum$coefficients["DiagnosisOCD", "Estimate"] #Response_statusresponder
  results[1, 2] <- sum$coefficients["DiagnosisOCD", "Std. Error"]
  results[1, 3] <- sum$coefficients["DiagnosisOCD", "t value"]
  results[1, 4] <- sum$coefficients["DiagnosisOCD", "Pr(>|t|)"]
  results[1, 5] <- clock
  
  
  results_CaseCtrl <- rbind(results_CaseCtrl, results)
}

# adjust p values
results_CaseCtrl$adj_p <- results_CaseCtrl$p_value * Meff
results_CaseCtrl$adj_p <- ifelse(results_CaseCtrl$adj_p > 1, 1, results_CaseCtrl$adj_p)
write_xlsx(results_CaseCtrl, path = "output/Tables/Results_CaseCtrl.xlsx")




################################################################################
# INCL SMOKING/BMI TO THE MODEL
################################################################################

Results_CaseCtrl_BMIsmoking <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 5))
colnames(Results_CaseCtrl_BMIsmoking) <- c("Estimate", "SE", "t_value", "p_value", "clock")

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  AgeAccelResid <- DNAmAge[which(DNAmAge$Clock == clock), "AgeAccelResid"]
  
  
  model <- lmrob(AgeAccelResid ~ Diagnosis + Age + Sex + BMI + Smoking)
  sum <- summary(model)
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 5))
  colnames(results) <- c("Estimate", "SE", "t_value", "p_value", "clock")
  
  results[1, 1] <- sum$coefficients["DiagnosisOCD", "Estimate"] #Response_statusresponder
  results[1, 2] <- sum$coefficients["DiagnosisOCD", "Std. Error"]
  results[1, 3] <- sum$coefficients["DiagnosisOCD", "t value"]
  results[1, 4] <- sum$coefficients["DiagnosisOCD", "Pr(>|t|)"]
  results[1, 5] <- clock
  
  
  Results_CaseCtrl_BMIsmoking <- rbind(Results_CaseCtrl_BMIsmoking, results)
}

# adjust p values
Results_CaseCtrl_BMIsmoking$adj_p <- Results_CaseCtrl_BMIsmoking$p_value * Meff
Results_CaseCtrl_BMIsmoking$adj_p <- ifelse(Results_CaseCtrl_BMIsmoking$adj_p > 1, 1, Results_CaseCtrl_BMIsmoking$adj_p)
write_xlsx(Results_CaseCtrl_BMIsmoking, path = "output/Tables/Results_CaseCtrl_BMIsmoking.xlsx")



################################################################################
# INCL SMOKING/BMI TO THE MODEL
################################################################################

Results_CaseCtrl_SmokingOnly <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 5))
colnames(Results_CaseCtrl_SmokingOnly) <- c("Estimate", "SE", "t_value", "p_value", "clock")

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  AgeAccelResid <- DNAmAge[which(DNAmAge$Clock == clock), "AgeAccelResid"]
  
  
  model <- lmrob(AgeAccelResid ~ Diagnosis + Age + Sex + Smoking)
  sum <- summary(model)
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 5))
  colnames(results) <- c("Estimate", "SE", "t_value", "p_value", "clock")
  
  results[1, 1] <- sum$coefficients["DiagnosisOCD", "Estimate"] #Response_statusresponder
  results[1, 2] <- sum$coefficients["DiagnosisOCD", "Std. Error"]
  results[1, 3] <- sum$coefficients["DiagnosisOCD", "t value"]
  results[1, 4] <- sum$coefficients["DiagnosisOCD", "Pr(>|t|)"]
  results[1, 5] <- clock
  
  
  Results_CaseCtrl_SmokingOnly <- rbind(Results_CaseCtrl_SmokingOnly, results)
}

# adjust p values
Results_CaseCtrl_SmokingOnly$adj_p <- Results_CaseCtrl_SmokingOnly$p_value * Meff
Results_CaseCtrl_SmokingOnly$adj_p <- ifelse(Results_CaseCtrl_SmokingOnly$adj_p > 1, 1, Results_CaseCtrl_SmokingOnly$adj_p)
write_xlsx(Results_CaseCtrl_SmokingOnly, path = "output/Tables/Results_CaseCtrl_SmokingOnly.xlsx")



################################################################################
# INCL BMI TO THE MODEL
################################################################################

Results_CaseCtrl_BMIOnly <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 5))
colnames(Results_CaseCtrl_BMIOnly) <- c("Estimate", "SE", "t_value", "p_value", "clock")

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  AgeAccelResid <- DNAmAge[which(DNAmAge$Clock == clock), "AgeAccelResid"]
  
  
  model <- lmrob(AgeAccelResid ~ Diagnosis + Age + Sex + BMI)
  sum <- summary(model)
  
  results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 5))
  colnames(results) <- c("Estimate", "SE", "t_value", "p_value", "clock")
  
  results[1, 1] <- sum$coefficients["DiagnosisOCD", "Estimate"] #Response_statusresponder
  results[1, 2] <- sum$coefficients["DiagnosisOCD", "Std. Error"]
  results[1, 3] <- sum$coefficients["DiagnosisOCD", "t value"]
  results[1, 4] <- sum$coefficients["DiagnosisOCD", "Pr(>|t|)"]
  results[1, 5] <- clock
  
  
  Results_CaseCtrl_BMIOnly <- rbind(Results_CaseCtrl_BMIOnly, results)
}

# adjust p values
Results_CaseCtrl_BMIOnly$adj_p <- Results_CaseCtrl_BMIOnly$p_value * Meff
Results_CaseCtrl_BMIOnly$adj_p <- ifelse(Results_CaseCtrl_BMIOnly$adj_p > 1, 1, Results_CaseCtrl_BMIOnly$adj_p)
write_xlsx(Results_CaseCtrl_BMIOnly, path = "output/Tables/Results_CaseCtrl_BMIOnly.xlsx")




################################################################################
# Visualisation
################################################################################

# Adjusted

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  DNAmAge_red <- DNAmAge[which(DNAmAge$Clock == clock), ]
  
  merged_DNAmAge <- merge(DNAmAge_red, samplesheet, by.x = "Sample", by.y = "Basename", all.x = TRUE)
  
  clockp <- star_p_value(results_CaseCtrl[which(results_CaseCtrl$clock == clock), "adj_p"])
  diff <- abs(max(merged_DNAmAge$AgeAccelResid) - min(merged_DNAmAge$AgeAccelResid))
  
  Age_plot <- ggplot(merged_DNAmAge, aes(x = Diagnosis, y = AgeAccelResid, fill = Diagnosis)) +
    geom_boxplot(alpha = 0.7) +
    th + th_transparent +
    labs(
      x = NULL,
      y = "AgeAccelResid"
    ) + 
    ggtitle(clock) + 
    theme(legend.position = "top", legend.title = element_blank(), legend.box.background = element_blank(), legend.key = element_blank()) +
    scale_fill_manual(values=c("CTRL"= "#084081", "OCD" = "#891a1a")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    annotate("segment", x = 1, xend = 2, y = max(merged_DNAmAge$AgeAccelResid) + 0.05*diff, 
             yend = max(merged_DNAmAge$AgeAccelResid)  + 0.05*diff, color = "black") +
    annotate("text", x = 1.5, y = max(merged_DNAmAge$AgeAccelResid)  + 0.095*diff, label = clockp, size = 5)
  Age_plot
  
  ggsave(paste0("output/Figures/CaseCtrl/CaseCtrl_AgeAccel_", clock, ".pdf"), Age_plot, units = "cm", width = 9, height = 10, dpi = 400)
  save(Age_plot, file = paste0("output/Figures/CaseCtrl//CaseCtrl_AgeAccel_", clock, ".RData"))
}

#raw
for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  DNAmAge_red <- DNAmAge[which(DNAmAge$Clock == clock), ]
  
  merged_DNAmAge <- merge(DNAmAge_red, samplesheet, by.x = "Sample", by.y = "Basename", all.x = TRUE)
  
  clockp <- star_p_value(results_CaseCtrl[which(results_CaseCtrl$clock == clock), "adj_p"])
  diff <- abs(max(merged_DNAmAge$AgeAccelResid) - min(merged_DNAmAge$Age_Acceleration))
  
  Age_plot <- ggplot(merged_DNAmAge, aes(x = Diagnosis, y = Age_Acceleration, fill = Diagnosis)) +
    geom_boxplot(alpha = 0.7) +
    th + th_transparent +
    labs(
      x = NULL,
      y = "AgeAccelResid"
    ) + 
    ggtitle(clock) + 
    theme(legend.position = "top", legend.title = element_blank(), legend.box.background = element_blank(), legend.key = element_blank()) +
    scale_fill_manual(values=c("CTRL"= "#084081", "OCD" = "#891a1a")) +
    theme(plot.title = element_text(hjust = 0.5))  +
    annotate("segment", x = 1, xend = 2, y = max(merged_DNAmAge$Age_Acceleration) + 0.05*diff, 
             yend = max(merged_DNAmAge$Age_Acceleration)  + 0.05*diff, color = "black") +
    annotate("text", x = 1.5, y = max(merged_DNAmAge$Age_Acceleration)  + 0.095*diff, label = clockp, size = 5)
  Age_plot
  
  ggsave(paste0("output/Figures/CaseCtrl/CaseCtrl_AgeAccel_", clock, "raw.pdf"), Age_plot, units = "cm", width = 9, height = 10, dpi = 400)
  save(Age_plot, file = paste0("output/Figures/CaseCtrl//CaseCtrl_AgeAccel_", clock, "raw.RData"))
}


################################################################################
# Sensitivity - GrimAge2
################################################################################

# LOAD DNAmAge
DNAmAge <- as.data.frame(read_xlsx("output/Tables/DNAmAge_All.xlsx"))
clocks <- c("GrimAge2")
DNAmAge <- DNAmAge[which(DNAmAge$Clock %in% clocks & DNAmAge$Sample %in% samplesheet$Basename), ]

AgeAccelResid <- DNAmAge[which(DNAmAge$Clock == "GrimAge2"), "AgeAccelResid"]
model <- lmrob(AgeAccelResid ~ Diagnosis + Age + Sex)
sum <- as.data.frame(summary(model)$coefficients)
sum$Variables <- rownames(sum)
write_xlsx(sum, path = "output/Tables/Sensitivity_CaseCtrl_GrimAge2_basic.xlsx")
sum

# Estimate Std. Error    t value     Pr(>|t|)    Variables
# (Intercept)   0.7879827  0.1949535  4.0418996 5.836796e-05  (Intercept)
# DiagnosisOCD -1.5545754  0.2223221 -6.9924468 5.882469e-12 DiagnosisOCD



model <- lmrob(AgeAccelResid ~ Diagnosis + Age + Sex + BMI + Smoking)
sum <- as.data.frame(summary(model)$coefficients)
sum$Variables <- rownames(sum)
write_xlsx(sum, path = "output/Tables/Sensitivity_CaseCtrl_GrimAge2_BMIsmoking.xlsx")
sum
# Estimate Std. Error   t value      Pr(>|t|)    Variables
# (Intercept)   0.73477296 0.14273045  5.147976  3.351147e-07  (Intercept)
# DiagnosisOCD -0.76585537 0.18412880 -4.159346  3.553178e-05 DiagnosisOCD
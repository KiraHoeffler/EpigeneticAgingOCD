

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
samplesheet <- as.data.frame(read_xlsx("output/Tables/Samplesheet_CaseFollowUpCtrl.xlsx"))

# LOAD THEME
load("Resources/theme.RData")
load("Resources/theme_transparent.Rdata")
source("Resources/Function_star_pvalue.R")

# LOAD DNAmAge
DNAmAge <- as.data.frame(read_xlsx("output/Tables/DNAmAge_All.xlsx"))

# nr of independent tests
load("output/RData/Meff_CaseFollowUpCtrl.RData")



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

results_CaseFollowUpCtrl <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 5))
colnames(results_CaseFollowUpCtrl) <- c("Estimate", "SE", "t_value", "p_value", "clock")

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
  
  
  results_CaseFollowUpCtrl <- rbind(results_CaseFollowUpCtrl, results)
}

# adjust p values
results_CaseFollowUpCtrl$adj_p <- results_CaseFollowUpCtrl$p_value * Meff
results_CaseFollowUpCtrl$adj_p <- ifelse(results_CaseFollowUpCtrl$adj_p > 1, 1, results_CaseFollowUpCtrl$adj_p)
write_xlsx(results_CaseFollowUpCtrl, path = "output/Tables/Results_CaseFollowUpCtrl.xlsx")




################################################################################
# INCL SMOKING/BMI TO THE MODEL
################################################################################

Results_CaseFollowUpCtrl_BMIsmoking <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 5))
colnames(Results_CaseFollowUpCtrl_BMIsmoking) <- c("Estimate", "SE", "t_value", "p_value", "clock")

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
  
  
  Results_CaseFollowUpCtrl_BMIsmoking <- rbind(Results_CaseFollowUpCtrl_BMIsmoking, results)
}

# adjust p values
Results_CaseFollowUpCtrl_BMIsmoking$adj_p <- Results_CaseFollowUpCtrl_BMIsmoking$p_value * Meff
Results_CaseFollowUpCtrl_BMIsmoking$adj_p <- ifelse(Results_CaseFollowUpCtrl_BMIsmoking$adj_p > 1, 1, Results_CaseFollowUpCtrl_BMIsmoking$adj_p)
write_xlsx(Results_CaseFollowUpCtrl_BMIsmoking, path = "output/Tables/Results_CaseFollowUpCtrl_BMIsmoking.xlsx")



################################################################################
# INCL SMOKING TO THE MODEL
################################################################################

Results_CaseFollowUpCtrl_SmokingOnly <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 5))
colnames(Results_CaseFollowUpCtrl_SmokingOnly) <- c("Estimate", "SE", "t_value", "p_value", "clock")

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
  
  
  Results_CaseFollowUpCtrl_SmokingOnly <- rbind(Results_CaseFollowUpCtrl_SmokingOnly, results)
}

# adjust p values
Results_CaseFollowUpCtrl_SmokingOnly$adj_p <- Results_CaseFollowUpCtrl_SmokingOnly$p_value * Meff
Results_CaseFollowUpCtrl_SmokingOnly$adj_p <- ifelse(Results_CaseFollowUpCtrl_SmokingOnly$adj_p > 1, 1, Results_CaseFollowUpCtrl_SmokingOnly$adj_p)
write_xlsx(Results_CaseFollowUpCtrl_SmokingOnly, path = "output/Tables/Results_CaseFollowUpCtrl_SmokingOnly.xlsx")



################################################################################
# INCL BMI TO THE MODEL
################################################################################

Results_CaseFollowUpCtrl_BMIOnly <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 5))
colnames(Results_CaseFollowUpCtrl_BMIOnly) <- c("Estimate", "SE", "t_value", "p_value", "clock")

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
  
  
  Results_CaseFollowUpCtrl_BMIOnly <- rbind(Results_CaseFollowUpCtrl_BMIOnly, results)
}

# adjust p values
Results_CaseFollowUpCtrl_BMIOnly$adj_p <- Results_CaseFollowUpCtrl_BMIOnly$p_value * Meff
Results_CaseFollowUpCtrl_BMIOnly$adj_p <- ifelse(Results_CaseFollowUpCtrl_BMIOnly$adj_p > 1, 1, Results_CaseFollowUpCtrl_BMIOnly$adj_p)
write_xlsx(Results_CaseFollowUpCtrl_BMIOnly, path = "output/Tables/Results_CaseFollowUpCtrl_BMIOnly.xlsx")



################################################################################
# Visualisation
################################################################################

# Adjusted

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  DNAmAge_red <- DNAmAge[which(DNAmAge$Clock == clock), ]
  
  merged_DNAmAge <- merge(DNAmAge_red, samplesheet, by.x = "Sample", by.y = "Basename", all.x = TRUE)
  
  clockp <- star_p_value(results_CaseFollowUpCtrl[which(results_CaseFollowUpCtrl$clock == clock), "adj_p"])
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
    theme(plot.title = element_text(hjust = 0.5))+
    annotate("segment", x = 1, xend = 2, y = max(merged_DNAmAge$AgeAccelResid) + 0.05*diff, 
             yend = max(merged_DNAmAge$AgeAccelResid)  + 0.05*diff, color = "black") +
    annotate("text", x = 1.5, y = max(merged_DNAmAge$AgeAccelResid)  + 0.095*diff, label = clockp, size = 5)
  Age_plot
  
  ggsave(paste0("output/Figures/CaseFollowUpCtrl/CaseFollowUpCtrl_AgeAccel_", clock, ".pdf"), Age_plot, units = "cm", width = 9, height = 10, dpi = 400)
  save(Age_plot, file = paste0("output/Figures/CaseFollowUpCtrl//CaseFollowUpCtrl_AgeAccel_", clock, ".RData"))
}

#raw
for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  DNAmAge_red <- DNAmAge[which(DNAmAge$Clock == clock), ]
  
  merged_DNAmAge <- merge(DNAmAge_red, samplesheet, by.x = "Sample", by.y = "Basename", all.x = TRUE)
  
  clockp <- star_p_value(results_CaseFollowUpCtrl[which(results_CaseFollowUpCtrl$clock == clock), "adj_p"])
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
    theme(plot.title = element_text(hjust = 0.5)) +
    annotate("segment", x = 1, xend = 2, y = max(merged_DNAmAge$Age_Acceleration) + 0.05*diff, 
             yend = max(merged_DNAmAge$Age_Acceleration)  + 0.05*diff, color = "black") +
    annotate("text", x = 1.5, y = max(merged_DNAmAge$Age_Acceleration)  + 0.095*diff, label = clockp, size = 5)
  Age_plot
  
  ggsave(paste0("output/Figures/CaseFollowUpCtrl/CaseFollowUpCtrl_AgeAccel_", clock, "raw.pdf"), Age_plot, units = "cm", width = 9, height = 10, dpi = 400)
  save(Age_plot, file = paste0("output/Figures/CaseFollowUpCtrl//CaseFollowUpCtrl_AgeAccel_", clock, "raw.RData"))
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
write_xlsx(sum, path = "output/Tables/Sensitivity_CaseFollowUpCtrl_GrimAge2_basic.xlsx")
sum

# Estimate Std. Error    t value     Pr(>|t|)    Variables
# (Intercept)   0.89499453  0.1984700  4.5094708 7.525207e-06  (Intercept)
# DiagnosisOCD -1.30852677  0.2259433 -5.7913951 1.021591e-08 DiagnosisOCD



model <- lmrob(AgeAccelResid ~ Diagnosis + Age + Sex + BMI + Smoking)
sum <- as.data.frame(summary(model)$coefficients)
sum$Variables <- rownames(sum)
write_xlsx(sum, path = "output/Tables/Sensitivity_CaseFollowUpCtrl_GrimAge2_BMIsmoking.xlsx")
sum
# Estimate Std. Error    t value      Pr(>|t|)    Variables
# (Intercept)   0.88082233 0.14311944  6.1544561  1.220035e-09  (Intercept)
# DiagnosisOCD -0.60375614 0.18233985 -3.3111585  9.730973e-04 DiagnosisOCD

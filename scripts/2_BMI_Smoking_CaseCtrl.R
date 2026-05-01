
################################################################################
# SETUP
################################################################################

work_dir <- "S://Project/WP-epigenetics/06_DNAmAge/"
samples_path <- "S:/Project/WP-epigenetics/04_Pipeline_combinedPhases/samplesheet_BMI.xlsx"
ctrl_path <- "S://Project/WP-epigenetics/04_Pipeline_combinedPhases/output/RData/Positive_ctrlprobe_intensities.RData"
SNP_path <- "S://Project/WP-epigenetics/04_Pipeline_combinedPhases/output/RData/comb_SNPs.RData"

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
library(tidyverse)
library(dplyr)
library(ggplot2)
library(writexl)
library(tictoc)
library(robustbase)


# LOAD SAMPLESHEET
samplesheet <- as.data.frame(read_xlsx(samples_path))

# LOAD
load("Resources/theme.RData") #ggplot theme
load("Resources/theme_transparent.Rdata") #ggplot theme
source("Resources/Function_star_pvalue.R")

# LOAD RESIDUALISED SNP BETAS
load(SNP_path)

# LOAD CONTROL PROBES
load(ctrl_path)


################################################################################
# PREPARE
################################################################################

# FILTER
samplesheet <- samplesheet[which((samplesheet$Diagnosis == "CTRL" | samplesheet$Time_point == "D1")), ]
Slide_Diagnosis <- as.data.frame(table(samplesheet$Sentrix_ID, samplesheet$Diagnosis))
Slide_Diagnosis_wide <- as.data.frame(Slide_Diagnosis %>%
  pivot_wider(names_from = "Var2", values_from = "Freq"))
Slide_Diagnosis_wide <- Slide_Diagnosis[which(Slide_Diagnosis_wide$CTRL > 0 & Slide_Diagnosis_wide$OCD > 0), ]
samplesheet <- samplesheet[which(samplesheet$Sentrix_ID %in% Slide_Diagnosis_wide$Var1), ]
table(samplesheet$Diagnosis)
# CTRL  OCD 
# 384  387


################################################################################
# CALCULATE RESIDUALIZED AGE ACCELERATION
# (correcting for covariates)
################################################################################

# COVARIATES TO ADJUST FOR
Age <- scale(samplesheet$Age)
Sex <- as.factor(samplesheet$Sex)
AMP_Plate <- as.factor(samplesheet$AMP_Plate)

# ancestry
red_anc <- comb_SNPs[, samplesheet$Basename]
anc_pca <- prcomp(t(red_anc))
anc_pcs <- as.data.frame(anc_pca$x[,1:5])
rm(red_anc, anc_pca)

anc_pc1 <- scale(anc_pcs$PC1)

# control probes (not PC4 as multi-colinearity with sex)
ctrl_red <- ctrl[, samplesheet$Basename]
ctrl_pca <- prcomp(t(ctrl_red))
ctrl_pcs <- as.data.frame(ctrl_pca$x[,1:15])
rm(ctrl_red, ctrl_pca)

ctrl_pc1 <- scale(ctrl_pcs$PC1)
ctrl_pc2 <- scale(ctrl_pcs$PC2)
ctrl_pc3 <- scale(ctrl_pcs$PC3)
ctrl_pc5 <- scale(ctrl_pcs$PC5)
ctrl_pc6 <- scale(ctrl_pcs$PC6)
ctrl_pc7 <- scale(ctrl_pcs$PC7)
ctrl_pc8 <- scale(ctrl_pcs$PC8)
ctrl_pc9 <- scale(ctrl_pcs$PC9)
ctrl_pc10 <- scale(ctrl_pcs$PC10)

# cell type proportions
samples_red <- samplesheet[, c("Epi", "Fib", "comb_ICs")]
ct_pca <- prcomp(samples_red)
ct_pcs <- as.data.frame(ct_pca$x[,1:2])
rm(samples_red, ct_pca)
cell_type1 <- scale(ct_pcs$PC1)
cell_type2 <- scale(ct_pcs$PC2)


################################################################################
# RESIDUALIZE BMI + SMOKING SCORES
################################################################################

#BMI
BMI <- samplesheet$epigenetic_bmi_score_do_2023
model <- lm(BMI ~ Age + Sex + cell_type1 + cell_type2 + AMP_Plate + 
              ctrl_pc1 + ctrl_pc2 + ctrl_pc3 + ctrl_pc5 + ctrl_pc6 + ctrl_pc7 + ctrl_pc8 + ctrl_pc9 + ctrl_pc10 + anc_pc1)
samplesheet$BMIresid <- residuals(model)

#Smoking
smoking <- samplesheet$smoking
model <- lm(smoking ~ Age + Sex + cell_type1 + cell_type2 + AMP_Plate + 
              ctrl_pc1 + ctrl_pc2 + ctrl_pc3 + ctrl_pc5 + ctrl_pc6 + ctrl_pc7 + ctrl_pc8 + ctrl_pc9 + ctrl_pc10 + anc_pc1)
samplesheet$Smokingresid <- residuals(model)

#EXPORT
write_xlsx(samplesheet, path = "output/Tables/Samplesheet_CaseCtrl.xlsx")


################################################################################
# COMPARE BMI + SMOKING
################################################################################

Diagnosis <- factor(samplesheet$Diagnosis)
BMIresid <- samplesheet$BMIresid
Smokingresid <- samplesheet$Smokingresid

#BMI
model <- lmrob(BMIresid ~ Diagnosis + Age + Sex)
BMI_sum <- as.data.frame(summary(model)$coefficients)
BMI_sum$Variable <- rownames(BMI_sum)
write_xlsx(BMI_sum, path = "output/Tables/BMI_CaseCtrl.xlsx")
# Estimate   Std. Error    t value  Pr(>|t|)     Variable
# (Intercept)   4.409188e-05 1.212119e-04  0.3637588 0.7161384  (Intercept)
# DiagnosisOCD -2.381285e-04 1.585981e-04 -1.5014585 0.1336486 DiagnosisOCD
BMI_sum


#Smoking
model <- lmrob(smoking ~ Diagnosis + Age + Sex)
Smoking_sum <- as.data.frame(summary(model)$coefficients)
Smoking_sum$Variable <- rownames(Smoking_sum)
write_xlsx(Smoking_sum, path = "output/Tables/Smoking_CaseCtrl.xlsx")
Smoking_sum

# Estimate Std. Error    t value      Pr(>|t|)     Variable
# (Intercept)  -1.70107459 0.03662188 -46.449681 4.067456e-225  (Intercept)
# DiagnosisOCD -0.24865644 0.04161569  -5.975064  3.519708e-09 DiagnosisOCD

# Visualize

BMIp <- star_p_value(BMI_sum["DiagnosisOCD", "Pr(>|t|)"])

BMI_plot <- ggplot(samplesheet, aes(x = Diagnosis, y = BMIresid, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.7) +
  th + th_transparent +
  labs(
    x = NULL,
    y = "epigenetic BMI (residuals)"
  ) + 
  theme(legend.position = "top", legend.title = element_blank(), legend.box.background = element_blank(), legend.key = element_blank()) +
  scale_fill_manual(values=c("CTRL"= "#084081", "OCD" = "#891a1a")) +
  annotate("segment", x = 1, xend = 2, y = max(samplesheet$BMIresid) + 0.002, 
           yend = max(samplesheet$BMIresid) + 0.002, color = "black") +
  annotate("text", x = 1.5, y = max(samplesheet$BMIresid) + 0.003, label = BMIp, size = 5)
BMI_plot

ggsave("output/Figures//BMIsmoking/BMI_CaseCtrl.pdf", BMI_plot, units = "cm", width = 9, height = 10, dpi = 400)
save(BMI_plot, file = "output/Figures/BMIsmoking/BMI_CaseCtrl.RData")







Smokingp <- star_p_value(Smoking_sum["DiagnosisOCD", "Pr(>|t|)"])

Smoking_plot <- ggplot(samplesheet, aes(x = Diagnosis, y = Smokingresid, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.7) +
  th + th_transparent +
  labs(
    x = NULL,
    y = "epigenetic Smoking (residuals)"
  ) + 
  theme(legend.position = "top", legend.title = element_blank(), legend.box.background = element_blank(), legend.key = element_blank()) +
  scale_fill_manual(values=c("CTRL"= "#084081", "OCD" = "#891a1a")) +
  annotate("segment", x = 1, xend = 2, y = max(samplesheet$Smokingresid) + 0.4, 
           yend = max(samplesheet$Smokingresid) + 0.4, color = "black") +
  annotate("text", x = 1.5, y = max(samplesheet$Smokingresid) + 0.7, label = Smokingp, size = 5)
Smoking_plot

ggsave("output/Figures/BMIsmoking//Smoking_CaseCtrl.pdf", Smoking_plot, units = "cm", width = 9, height = 10, dpi = 400)
save(Smoking_plot, file = "output/Figures/BMIsmoking//Smoking_CaseCtrl.RData")

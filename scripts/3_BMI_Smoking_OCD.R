
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
library(lme4)
library(lmerTest)
library(robustlmm)


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
samplesheet <- samplesheet[which(samplesheet$Diagnosis == "OCD"), ]
table(samplesheet$Time_point)
# D1  D4  M3 
# 899 405 404
length(unique(samplesheet$Resp_nr)) #908


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
write_xlsx(samplesheet, path = "output/Tables/Samplesheet_OCD.xlsx")



################################################################################
# LONGITUDINAL SAMPLES
################################################################################


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



#### COVARIATES

# COVARIATES TO ADJUST FOR
Age <- scale(samples_long$Age)
Sex <- as.factor(samples_long$Sex)
Smoking <- samples_long$Smokingresid
BMI <- samples_long$BMIresid
Time_point <- as.factor(samples_long$Time_point)
Resp_nr <- as.factor(samples_long$Resp_nr)



################################################################################
# BMI ANALYSIS
################################################################################

model <- rlmer(formula = BMI ~ Time_point + Age + Sex + (1|Resp_nr))
BMI_sum <- as.data.frame(summary(model)$coefficients)
BMI_sum$Variable <- rownames(BMI_sum)

#compute p values (robust Wald test)
coefs <- summary(model)$coefficients
test_stat <- coefs[, "Estimate"] / coefs[, "Std. Error"]
p_values <- 2 * pnorm(abs(test_stat), lower.tail = FALSE)

BMI_sum$p_value <- p_values

write_xlsx(BMI_sum, path = "output/Tables/BMI_OCD.xlsx")
BMI_sum

# Estimate   Std. Error    t value     Variable   p_values
# (Intercept)  -1.319719e-04 1.260499e-04 -1.0469820  (Intercept) 0.29510787
# Time_pointD4  1.130030e-04 7.348568e-05  1.5377554 Time_pointD4 0.12410843
# Time_pointM3  1.641924e-04 7.379110e-05  2.2250978 Time_pointM3 0.02607467



################################################################################
# Smoking ANALYSIS
################################################################################

model <- rlmer(formula = Smoking ~ Time_point + Age + Sex + (1|Resp_nr))
Smoking_sum <- as.data.frame(summary(model)$coefficients)
Smoking_sum$Variable <- rownames(Smoking_sum)


#compute p values (robust Wald test)
coefs <- summary(model)$coefficients
test_stat <- coefs[, "Estimate"] / coefs[, "Std. Error"]
p_values <- 2 * pnorm(abs(test_stat), lower.tail = FALSE)

Smoking_sum$p_value <- p_values

write_xlsx(Smoking_sum, path = "output/Tables/Smoking_OCD.xlsx")
Smoking_sum

# Estimate Std. Error     t value     Variable     p_values
# (Intercept)  -0.0580069708 0.01361802 -4.25957508  (Intercept) 2.048159e-05
# Time_pointD4  0.0036001436 0.01100673  0.32708564 Time_pointD4 7.436031e-01
# Time_pointM3  0.0005414354 0.01104693  0.04901229 Time_pointM3 9.609095e-01



################################################################################
# Visualisation
################################################################################

# BMI

BMIp_D4 <- star_p_value(BMI_sum["Time_pointD4", "p_value"])
BMIp_M3 <- star_p_value(BMI_sum["Time_pointM3", "p_value"])

df_avg <- as.data.frame(samples_long %>%
                          group_by(Time_point) %>%
                          summarise(
                            avg = mean(BMIresid, na.rm = TRUE),
                            se = sd(BMIresid, na.rm = TRUE) / sqrt(n())))

df_avg$Time_point[which(df_avg$Time_point == "D1")] <- "pre"
df_avg$Time_point[which(df_avg$Time_point == "D4")] <- "post"
df_avg$Time_point[which(df_avg$Time_point == "M3")] <- "follow-up"

df_avg$Time_point <- factor(df_avg$Time_point, levels = c("pre", "post", "follow-up"))

BMI_plot <- ggplot(df_avg, aes(x = Time_point, y = avg, group = 1)) +
  geom_point(size = 3, color = "black") +
  geom_line(color = "black") +
  geom_errorbar(aes(ymin=avg-se, ymax=avg+se), width = 0.2) +
  th + th_transparent +
  labs(
    x = "Time point",
    y = "epigenetic BMI (residuals)"
  ) + 
  ggtitle("Epigenetic BMI") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme(legend.position = "top", legend.title = element_blank(), legend.box.background = element_blank(), legend.key = element_blank())  + 
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("segment", x = 1, xend = 2, y = max(df_avg$avg) + 0.00013, 
           yend = max(df_avg$avg) + 0.00013, color = "black") +
  annotate("text", x = 1.5, y = max(df_avg$avg) + 0.00016, label = BMIp_D4, size = 5) +
  annotate("segment", x = 1, xend = 3, y = max(df_avg$avg) + 0.00018, 
           yend = max(df_avg$avg) + 0.00018, color = "black") +
  annotate("text", x = 2, y = max(df_avg$avg) + 0.0002, label = BMIp_M3, size = 5)
BMI_plot

ggsave("output/Figures/BMIsmoking/BMI_plot_OCD.pdf", BMI_plot, units = "cm", width = 10, height = 10, dpi = 400)


# Smoking

Smokingp_D4 <- star_p_value(Smoking_sum["Time_pointD4", "p_value"])
Smokingp_M3 <- star_p_value(Smoking_sum["Time_pointM3", "p_value"])

df_avg <- as.data.frame(samples_long %>%
                          group_by(Time_point) %>%
                          summarise(
                            avg = mean(Smokingresid, na.rm = TRUE),
                            se = sd(Smokingresid, na.rm = TRUE) / sqrt(n())))

df_avg$Time_point[which(df_avg$Time_point == "D1")] <- "pre"
df_avg$Time_point[which(df_avg$Time_point == "D4")] <- "post"
df_avg$Time_point[which(df_avg$Time_point == "M3")] <- "follow-up"

df_avg$Time_point <- factor(df_avg$Time_point, levels = c("pre", "post", "follow-up"))

smoking_plot <- ggplot(df_avg, aes(x = Time_point, y = avg, group = 1)) +
  geom_point(size = 3, color = "black") +
  geom_line(color = "black") +
  geom_errorbar(aes(ymin=avg-se, ymax=avg+se), width = 0.2) +
  th + th_transparent +
  labs(
    x = "Time point",
    y = "epigenetic smoking (residuals)"
  ) + 
  ggtitle("Epigenetic smoking") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme(legend.position = "top", legend.title = element_blank(), legend.box.background = element_blank(), legend.key = element_blank())  + 
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("segment", x = 1, xend = 2, y = max(df_avg$avg) + 0.023, 
           yend = max(df_avg$avg) + 0.023, color = "black") +
  annotate("text", x = 1.5, y = max(df_avg$avg) + 0.025, label = Smokingp_D4, size = 5) +
  annotate("segment", x = 1, xend = 3, y = max(df_avg$avg) + 0.029, 
           yend = max(df_avg$avg) + 0.029, color = "black") +
  annotate("text", x = 2, y = max(df_avg$avg) + 0.031, label = Smokingp_M3, size = 5)
smoking_plot

ggsave("output/Figures/BMIsmoking/Smoking_OCD.pdf", smoking_plot, units = "cm", width = 10, height = 10, dpi = 400)







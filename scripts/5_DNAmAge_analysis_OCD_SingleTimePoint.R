

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
samplesheet <- as.data.frame(read_xlsx("output/Tables/Samplesheet_OCD.xlsx"))
samplesheet <- samplesheet[which(!is.na(samplesheet$Percent_improvement)), ]
table(samplesheet$Time_point)
#D1  D4  M3 
#881 388 385 
length(unique(samplesheet$Resp_nr)) #889

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

# SELECT CLOCKS TO TEST
clocks <- c("PC_Horvath", "PC_GrimAge", "PC_PhenoAge", "DunedinPACE")
DNAmAge <- DNAmAge[which(DNAmAge$Clock %in% clocks), ]


# SAMPLES FOR DIFFERENT TIME POINTS
samples_D1 <- samplesheet$Basename[which(samplesheet$Time_point == "D1")]
samples_D4 <- samplesheet$Basename[which(samplesheet$Time_point == "D4")]
samples_M3 <- samplesheet$Basename[which(samplesheet$Time_point == "M3")]

samples_list <- list(samples_D1, samples_D4, samples_M3)
sample_names <- c("pre", "post", "follow-up")


################################################################################
# Response / Single Time point (BASIC MODEL)
################################################################################

results_singleTP <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 6))
colnames(results_singleTP) <- c("sample", "Estimate", "SE", "t_value", "p_value", "clock")

for (i in c(1:3)){
  sample_name <- sample_names[i]
  
  for (j in c(1:length(clocks))){
    clock <- clocks[j]
    
    samples_red <- samplesheet[which(samplesheet$Basename %in% samples_list[[i]]),]
    
    AgeAccelResid <- DNAmAge[which(DNAmAge$Sample %in% samples_list[[i]] & DNAmAge$Clock == clock), "AgeAccelResid"]
    Response <- scale(samples_red$Percent_improvement)
    Sex <- as.factor(samples_red$Sex)
    Age <- scale(samples_red$Age)
    
    model <- lmrob(AgeAccelResid ~ Response + Age + Sex)
    sum <- summary(model)
    
    results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
    colnames(results) <- c("sample", "Estimate", "SE", "t_value", "p_value", "clock")
    
    results[1, 1] <- sample_name
    results[1, 2] <- sum$coefficients["Response", "Estimate"] #Response_statusresponder
    results[1, 3] <- sum$coefficients["Response", "Std. Error"]
    results[1, 4] <- sum$coefficients["Response", "t value"]
    results[1, 5] <- sum$coefficients["Response", "Pr(>|t|)"]
    results[1, 6] <- clock
    
    
    results_singleTP <- rbind(results_singleTP, results)
  }
}

# adjust p values
results_singleTP$adj_p <- results_singleTP$p_value * Meff
results_singleTP$adj_p <- ifelse(results_singleTP$adj_p > 1, 1, results_singleTP$adj_p)
write_xlsx(results_singleTP, path = "output/Tables/Results_OCD_SingleTimePoints_Response.xlsx")



results_singleTP$sample <- factor(results_singleTP$sample, levels = c("pre", "post", "follow-up"))
results_singleTP$fill_color <- ifelse(results_singleTP$adj_p <= 0.1, "#d4f3fc", "white")
results_singleTP$fill_color <- ifelse(results_singleTP$adj_p <= 0.05, "#18a6b9", results_singleTP$fill_color)
results_singleTP$rounded_p <- round(results_singleTP$adj_p, 3)
results_singleTP$clock <- factor(results_singleTP$clock, levels = sort(unique(results_singleTP$clock), decreasing = TRUE))

singleTP_plot <- ggplot(results_singleTP, aes(x=sample, y=clock, fill=fill_color))+
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
  ggtitle("AgeAccelResid ~ Resp\nat single time point\n(adj. p values)")
singleTP_plot

ggsave("output/Figures/Results/Results_OCD_SingleTP_Response.pdf", singleTP_plot, units = "cm", width = 8, height = 8, dpi = 400)
ggsave("output/Figures/Results/Results_OCD_SingleTP_Response.svg", singleTP_plot, device = "svg", units = "cm", width = 8, height = 8, dpi = 400)
save(singleTP_plot, file = "output/Figures/Results/Results_OCD_SingleTP_Response.RData")



################################################################################
# Response / Single Time point -> INCL SMOKING/BMI TO THE MODEL
################################################################################


results_singleTP_BMIsmoking <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 6))
colnames(results_singleTP_BMIsmoking) <- c("sample", "Estimate", "SE", "t_value", "p_value", "clock")

for (i in c(1:3)){
  sample_name <- sample_names[i]
  
  for (j in c(1:length(clocks))){
    clock <- clocks[j]
    
    samples_red <- samplesheet[which(samplesheet$Basename %in% samples_list[[i]]),]
    
    AgeAccelResid <- DNAmAge[which(DNAmAge$Sample %in% samples_list[[i]] & DNAmAge$Clock == clock), "AgeAccelResid"]
    Response <- scale(samples_red$Percent_improvement)
    Sex <- as.factor(samples_red$Sex)
    Age <- scale(samples_red$Age)
    BMI <- scale(samples_red$BMIresid)
    Smoking <- scale(samples_red$Smokingresid)
    
    model <- lmrob(AgeAccelResid ~ Response + Age + Sex + BMI + Smoking)
    sum <- summary(model)
    
    results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
    colnames(results) <- c("sample", "Estimate", "SE", "t_value", "p_value", "clock")
    
    results[1, 1] <- sample_name
    results[1, 2] <- sum$coefficients["Response", "Estimate"] #Response_statusresponder
    results[1, 3] <- sum$coefficients["Response", "Std. Error"]
    results[1, 4] <- sum$coefficients["Response", "t value"]
    results[1, 5] <- sum$coefficients["Response", "Pr(>|t|)"]
    results[1, 6] <- clock
    
    
    results_singleTP_BMIsmoking <- rbind(results_singleTP_BMIsmoking, results)
  }
}

# adjust p values
results_singleTP_BMIsmoking$adj_p <- results_singleTP_BMIsmoking$p_value * Meff
results_singleTP_BMIsmoking$adj_p <- ifelse(results_singleTP_BMIsmoking$adj_p > 1, 1, results_singleTP_BMIsmoking$adj_p)
write_xlsx(results_singleTP_BMIsmoking, path = "output/Tables/Results_OCD_SingleTimePoints_Response_BMI_Smoking.xlsx")




results_singleTP_BMIsmoking$sample <- factor(results_singleTP_BMIsmoking$sample, levels = c("pre", "post", "follow-up"))
results_singleTP_BMIsmoking$fill_color <- ifelse(results_singleTP_BMIsmoking$adj_p <= 0.1, "#d4f3fc", "white")
results_singleTP_BMIsmoking$fill_color <- ifelse(results_singleTP_BMIsmoking$adj_p <= 0.05, "#18a6b9", results_singleTP_BMIsmoking$fill_color)
results_singleTP_BMIsmoking$rounded_p <- round(results_singleTP_BMIsmoking$adj_p, 3)
results_singleTP_BMIsmoking$clock <- factor(results_singleTP_BMIsmoking$clock, levels = sort(unique(results_singleTP_BMIsmoking$clock), decreasing = TRUE))

singleTP_plot_Smoking_BMI <- ggplot(results_singleTP_BMIsmoking, aes(x=sample, y=clock, fill=fill_color))+
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
  ggtitle("AgeAccelResid ~ Resp\nat single time point\n(adj. p values)")
singleTP_plot_Smoking_BMI

ggsave("output/Figures/Results/Results_OCD_singleTP_Response_BMIsmoking.pdf", singleTP_plot_Smoking_BMI, units = "cm", width = 8, height = 8, dpi = 400)
ggsave("output/Figures/Results/Results_OCD_singleTP_Response_BMIsmoking.svg", singleTP_plot_Smoking_BMI, device = "svg", units = "cm", width = 8, height = 8, dpi = 400)
save(singleTP_plot_Smoking_BMI, file = "output/Figures/Results/Results_OCD_singleTP_Response_BMIsmoking.RData")


################################################################################
# Response / Single Time point -> INCL ONLY SMOKING
################################################################################


results_singleTP_Smoking <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 6))
colnames(results_singleTP_Smoking) <- c("sample", "Estimate", "SE", "t_value", "p_value", "clock")

for (i in c(1:3)){
  sample_name <- sample_names[i]
  
  for (j in c(1:length(clocks))){
    clock <- clocks[j]
    
    samples_red <- samplesheet[which(samplesheet$Basename %in% samples_list[[i]]),]
    
    AgeAccelResid <- DNAmAge[which(DNAmAge$Sample %in% samples_list[[i]] & DNAmAge$Clock == clock), "AgeAccelResid"]
    Response <- scale(samples_red$Percent_improvement)
    Sex <- as.factor(samples_red$Sex)
    Age <- scale(samples_red$Age)
    Smoking <- scale(samples_red$Smokingresid)
    
    model <- lmrob(AgeAccelResid ~ Response + Age + Sex + Smoking)
    sum <- summary(model)
    
    results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
    colnames(results) <- c("sample", "Estimate", "SE", "t_value", "p_value", "clock")
    
    results[1, 1] <- sample_name
    results[1, 2] <- sum$coefficients["Response", "Estimate"] #Response_statusresponder
    results[1, 3] <- sum$coefficients["Response", "Std. Error"]
    results[1, 4] <- sum$coefficients["Response", "t value"]
    results[1, 5] <- sum$coefficients["Response", "Pr(>|t|)"]
    results[1, 6] <- clock
    
    
    results_singleTP_Smoking <- rbind(results_singleTP_Smoking, results)
  }
}

# adjust p values
results_singleTP_Smoking$adj_p <- results_singleTP_Smoking$p_value * Meff
results_singleTP_Smoking$adj_p <- ifelse(results_singleTP_Smoking$adj_p > 1, 1, results_singleTP_Smoking$adj_p)
write_xlsx(results_singleTP_Smoking, path = "output/Tables/Results_OCD_SingleTimePoints_Response_SmokingOnly.xlsx")




results_singleTP_Smoking$sample <- factor(results_singleTP_Smoking$sample, levels = c("pre", "post", "follow-up"))
results_singleTP_Smoking$fill_color <- ifelse(results_singleTP_Smoking$adj_p <= 0.1, "#d4f3fc", "white")
results_singleTP_Smoking$fill_color <- ifelse(results_singleTP_Smoking$adj_p <= 0.05, "#18a6b9", results_singleTP_Smoking$fill_color)
results_singleTP_Smoking$rounded_p <- round(results_singleTP_Smoking$adj_p, 3)
results_singleTP_Smoking$clock <- factor(results_singleTP_Smoking$clock, levels = sort(unique(results_singleTP_Smoking$clock), decreasing = TRUE))

singleTP_plot_Smoking <- ggplot(results_singleTP_Smoking, aes(x=sample, y=clock, fill=fill_color))+
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
  ggtitle("AgeAccelResid ~ Resp\nat single time point\n(adj. p values)")
singleTP_plot_Smoking

ggsave("output/Figures/Results/Results_OCD_singleTP_Response_SmokingOnly.pdf", singleTP_plot_Smoking, units = "cm", width = 8, height = 8, dpi = 400)
ggsave("output/Figures/Results/Results_OCD_singleTP_Response_SmokingOnly.svg", singleTP_plot_Smoking, device = "svg", units = "cm", width = 8, height = 8, dpi = 400)
save(singleTP_plot_Smoking, file = "output/Figures/Results/Results_OCD_singleTP_Response_SmokingOnly.RData")


################################################################################
# Response / Single Time point -> INCL ONLY BMI
################################################################################


results_singleTP_BMI <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 6))
colnames(results_singleTP_BMI) <- c("sample", "Estimate", "SE", "t_value", "p_value", "clock")

for (i in c(1:3)){
  sample_name <- sample_names[i]
  
  for (j in c(1:length(clocks))){
    clock <- clocks[j]
    
    samples_red <- samplesheet[which(samplesheet$Basename %in% samples_list[[i]]),]
    
    AgeAccelResid <- DNAmAge[which(DNAmAge$Sample %in% samples_list[[i]] & DNAmAge$Clock == clock), "AgeAccelResid"]
    Response <- scale(samples_red$Percent_improvement)
    Sex <- as.factor(samples_red$Sex)
    Age <- scale(samples_red$Age)
    BMI <- scale(samples_red$BMIresid)
    
    model <- lmrob(AgeAccelResid ~ Response + Age + Sex + BMI)
    sum <- summary(model)
    
    results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
    colnames(results) <- c("sample", "Estimate", "SE", "t_value", "p_value", "clock")
    
    results[1, 1] <- sample_name
    results[1, 2] <- sum$coefficients["Response", "Estimate"] #Response_statusresponder
    results[1, 3] <- sum$coefficients["Response", "Std. Error"]
    results[1, 4] <- sum$coefficients["Response", "t value"]
    results[1, 5] <- sum$coefficients["Response", "Pr(>|t|)"]
    results[1, 6] <- clock
    
    
    results_singleTP_BMI <- rbind(results_singleTP_BMI, results)
  }
}

# adjust p values
results_singleTP_BMI$adj_p <- results_singleTP_BMI$p_value * Meff
results_singleTP_BMI$adj_p <- ifelse(results_singleTP_BMI$adj_p > 1, 1, results_singleTP_BMI$adj_p)
write_xlsx(results_singleTP_BMI, path = "output/Tables/Results_OCD_SingleTimePoints_Response_BMIOnly.xlsx")




results_singleTP_BMI$sample <- factor(results_singleTP_BMI$sample, levels = c("pre", "post", "follow-up"))
results_singleTP_BMI$fill_color <- ifelse(results_singleTP_BMI$adj_p <= 0.1, "#d4f3fc", "white")
results_singleTP_BMI$fill_color <- ifelse(results_singleTP_BMI$adj_p <= 0.05, "#18a6b9", results_singleTP_BMI$fill_color)
results_singleTP_BMI$rounded_p <- round(results_singleTP_BMI$adj_p, 3)
results_singleTP_BMI$clock <- factor(results_singleTP_BMI$clock, levels = sort(unique(results_singleTP_BMI$clock), decreasing = TRUE))

singleTP_plot_BMI <- ggplot(results_singleTP_BMI, aes(x=sample, y=clock, fill=fill_color))+
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
  ggtitle("AgeAccelResid ~ Resp\nat single time point\n(adj. p values)")
singleTP_plot_BMI

ggsave("output/Figures/Results/Results_OCD_singleTP_Response_BMIOnly.pdf", singleTP_plot_BMI, units = "cm", width = 8, height = 8, dpi = 400)
ggsave("output/Figures/Results/Results_OCD_singleTP_Response_BMIOnly.svg", singleTP_plot_BMI, device = "svg", units = "cm", width = 8, height = 8, dpi = 400)
save(singleTP_plot_BMI, file = "output/Figures/Results/Results_OCD_singleTP_Response_BMIOnly.RData")



################################################################################
# Visualisation
################################################################################

### RESPONSE STATUS AT SINGLE TIME POINT (BOXPLOT)###

# PRE

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  samples_clin <- samplesheet[which(samplesheet$Case_Control == "Case"), c("Basename", "Resp_nr", "Response_status")]
  
  DNAmAge_red <- DNAmAge[which(DNAmAge$Clock == clock), ]
  DNAmAge_red <- DNAmAge_red[which(DNAmAge_red$Time_point == "D1"), ]
  
  
  merged_DNAmAge <- merge(DNAmAge_red, samples_clin, by.x = "Sample", by.y = "Basename", all.y = TRUE)
  merged_DNAmAge <- merged_DNAmAge[which(!is.na(merged_DNAmAge$Response_status)), ]
  merged_DNAmAge <- merged_DNAmAge[which(!is.na(merged_DNAmAge$Time_point)), ]
  
  Age_plot <- ggplot(merged_DNAmAge, aes(x = Response_status, y = AgeAccelResid, fill = Response_status)) +
    geom_boxplot(alpha = 0.7) +
    th + th_transparent +
    labs(
      x = NULL,
      y = "Age acceleration (residuals)"
    ) + 
    ggtitle(clock) + 
    theme(legend.position = "top", legend.title = element_blank(), legend.box.background = element_blank(), legend.key = element_blank()) +
    scale_fill_manual(values=c("responder"= "#084081", "non-responder" = "#891a1a")) +
    theme(plot.title = element_text(hjust = 0.5))
  Age_plot
  
  ggsave(paste0("output/Figures/SingleTimePoint/Baseline_AgeAccel_ResponseStatus_", clock, ".pdf"), Age_plot, units = "cm", width = 9, height = 10, dpi = 400)
  save(Age_plot, file = paste0("output/Figures/SingleTimePoint/Baseline_AgeAccel_ResponseStatus_", clock, ".RData"))
}

# POST

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  samples_clin <- samplesheet[which(samplesheet$Case_Control == "Case"), c("Basename", "Resp_nr", "Response_status")]
  
  DNAmAge_red <- DNAmAge[which(DNAmAge$Clock == clock), ]
  DNAmAge_red <- DNAmAge_red[which(DNAmAge_red$Time_point == "D4"), ]
  
  
  merged_DNAmAge <- merge(DNAmAge_red, samples_clin, by.x = "Sample", by.y = "Basename", all.y = TRUE)
  merged_DNAmAge <- merged_DNAmAge[which(!is.na(merged_DNAmAge$Response_status)), ]
  merged_DNAmAge <- merged_DNAmAge[which(!is.na(merged_DNAmAge$Time_point)), ]
  
  Age_plot <- ggplot(merged_DNAmAge, aes(x = Response_status, y = AgeAccelResid, fill = Response_status)) +
    geom_boxplot(alpha = 0.7) +
    th + th_transparent +
    labs(
      x = NULL,
      y = "Age acceleration (residuals)"
    ) + 
    ggtitle(clock) + 
    theme(legend.position = "top", legend.title = element_blank(), legend.box.background = element_blank(), legend.key = element_blank()) +
    scale_fill_manual(values=c("responder"= "#084081", "non-responder" = "#891a1a")) +
    theme(plot.title = element_text(hjust = 0.5))
  Age_plot
  
  ggsave(paste0("output/Figures/SingleTimePoint/Post_AgeAccel_ResponseStatus_", clock, ".pdf"), Age_plot, units = "cm", width = 9, height = 10, dpi = 400)
  save(Age_plot, file = paste0("output/Figures/SingleTimePoint/Post_AgeAccel_ResponseStatus_", clock, ".RData"))
}

# FOLLOW-UP

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  samples_clin <- samplesheet[which(samplesheet$Case_Control == "Case"), c("Basename", "Resp_nr", "Response_status")]
  
  DNAmAge_red <- DNAmAge[which(DNAmAge$Clock == clock), ]
  DNAmAge_red <- DNAmAge_red[which(DNAmAge_red$Time_point == "M3"), ]
  
  
  merged_DNAmAge <- merge(DNAmAge_red, samples_clin, by.x = "Sample", by.y = "Basename", all.y = TRUE)
  merged_DNAmAge <- merged_DNAmAge[which(!is.na(merged_DNAmAge$Response_status)), ]
  merged_DNAmAge <- merged_DNAmAge[which(!is.na(merged_DNAmAge$Time_point)), ]
  
  Age_plot <- ggplot(merged_DNAmAge, aes(x = Response_status, y = AgeAccelResid, fill = Response_status)) +
    geom_boxplot(alpha = 0.7) +
    th + th_transparent +
    labs(
      x = NULL,
      y = "Age acceleration (residuals)"
    ) + 
    ggtitle(clock) + 
    theme(legend.position = "top", legend.title = element_blank(), legend.box.background = element_blank(), legend.key = element_blank()) +
    scale_fill_manual(values=c("responder"= "#084081", "non-responder" = "#891a1a")) +
    theme(plot.title = element_text(hjust = 0.5))
  Age_plot
  
  ggsave(paste0("output/Figures/SingleTimePoint/FollowUp_AgeAccel_ResponseStatus_", clock, ".pdf"), Age_plot, units = "cm", width = 9, height = 10, dpi = 400)
  save(Age_plot, file = paste0("output/Figures/SingleTimePoint/FollowUp_AgeAccel_ResponseStatus_", clock, ".RData"))
}





### %IMPROVEMENT AT SINGLE TIME POINT (CORR) ###

# Pre

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  samples_clin <- samplesheet[which(samplesheet$Case_Control == "Case"), c("Basename", "Resp_nr", "Percent_improvement")]
  
  DNAmAge_red <- DNAmAge[which(DNAmAge$Clock == clock), ]
  DNAmAge_red <- DNAmAge_red[which(DNAmAge_red$Time_point == "D1"), ]
  
  
  merged_DNAmAge <- merge(DNAmAge_red, samples_clin, by.x = "Sample", by.y = "Basename", all.y = TRUE)
  merged_DNAmAge <- merged_DNAmAge[which(!is.na(merged_DNAmAge$Percent_improvement)), ]
  merged_DNAmAge <- merged_DNAmAge[which(!is.na(merged_DNAmAge$Time_point)), ]
  
  corr <- round(cor(merged_DNAmAge$Percent_improvement, merged_DNAmAge$AgeAccelResid), 3)
  
  Corr_plot <- ggplot(merged_DNAmAge, aes(x = Percent_improvement, y = AgeAccelResid)) +
    geom_point(color = "#084081", alpha = 0.15) + 
    geom_smooth(color = "#084081", method = "lm", linewidth=1.2, alpha = 0.35) +
    scale_color_identity() +
    ggtitle(clock) +
    th + th_transparent +
    xlab("% clinical improvement") + ylab("Age Acceleration (residuals)") +
    theme(plot.title = element_text(hjust = 0.5)) + 
    annotate("text", x=15, y = min(merged_DNAmAge$AgeAccelResid), label = paste0("r = ", corr), color = "black")
  
  ggsave(paste0("output/Figures/SingleTimePoint/Pre_AgeAccel_Response_", clock, ".pdf"), Corr_plot, units = "cm", width = 7, height = 10, dpi = 400)
  save(Corr_plot, file = paste0("output/Figures/SingleTimePoint/Pre_AgeAccel_Response_", clock, ".RData"))
}




# Post

for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  samples_clin <- samplesheet[which(samplesheet$Case_Control == "Case"), c("Basename", "Resp_nr", "Percent_improvement")]
  
  DNAmAge_red <- DNAmAge[which(DNAmAge$Clock == clock), ]
  DNAmAge_red <- DNAmAge_red[which(DNAmAge_red$Time_point == "D4"), ]
  
  
  merged_DNAmAge <- merge(DNAmAge_red, samples_clin, by.x = "Sample", by.y = "Basename", all.y = TRUE)
  merged_DNAmAge <- merged_DNAmAge[which(!is.na(merged_DNAmAge$Percent_improvement)), ]
  merged_DNAmAge <- merged_DNAmAge[which(!is.na(merged_DNAmAge$Time_point)), ]
  
  corr <- round(cor(merged_DNAmAge$Percent_improvement, merged_DNAmAge$AgeAccelResid), 3)
  
  Corr_plot <- ggplot(merged_DNAmAge, aes(x = Percent_improvement, y = AgeAccelResid)) +
    geom_point(color = "#084081", alpha = 0.15) + 
    geom_smooth(color = "#084081", method = "lm", linewidth=1.2, alpha = 0.35) +
    scale_color_identity() +
    ggtitle(clock) +
    th + th_transparent +
    xlab("% clinical improvement") + ylab("Age Acceleration (residuals)") +
    theme(plot.title = element_text(hjust = 0.5)) + 
    annotate("text", x=15, y = min(merged_DNAmAge$AgeAccelResid), label = paste0("r = ", corr), color = "black")
  
  ggsave(paste0("output/Figures/SingleTimePoint/Post_AgeAccel_Response_", clock, ".pdf"), Corr_plot, units = "cm", width = 7, height = 10, dpi = 400)
  save(Corr_plot, file = paste0("output/Figures/SingleTimePoint/Post_AgeAccel_Response_", clock, ".RData"))
}



# Follow-up
for (j in c(1:length(clocks))){
  clock <- clocks[j]
  
  samples_clin <- samplesheet[which(samplesheet$Case_Control == "Case"), c("Basename", "Resp_nr", "Percent_improvement")]
  
  DNAmAge_red <- DNAmAge[which(DNAmAge$Clock == clock), ]
  DNAmAge_red <- DNAmAge_red[which(DNAmAge_red$Time_point == "M3"), ]
  
  
  merged_DNAmAge <- merge(DNAmAge_red, samples_clin, by.x = "Sample", by.y = "Basename", all.y = TRUE)
  merged_DNAmAge <- merged_DNAmAge[which(!is.na(merged_DNAmAge$Percent_improvement)), ]
  merged_DNAmAge <- merged_DNAmAge[which(!is.na(merged_DNAmAge$Time_point)), ]
  
  corr <- round(cor(merged_DNAmAge$Percent_improvement, merged_DNAmAge$AgeAccelResid), 3)
  
  Corr_plot <- ggplot(merged_DNAmAge, aes(x = Percent_improvement, y = AgeAccelResid)) +
    geom_point(color = "#084081", alpha = 0.15) + 
    geom_smooth(color = "#084081", method = "lm", linewidth=1.2, alpha = 0.35) +
    scale_color_identity() +
    ggtitle(clock) +
    th + th_transparent +
    xlab("% clinical improvement") + ylab("Age Acceleration (residuals)") +
    theme(plot.title = element_text(hjust = 0.5)) + 
    annotate("text", x=15, y = min(merged_DNAmAge$AgeAccelResid), label = paste0("r = ", corr), color = "black")
  
  ggsave(paste0("output/Figures/SingleTimePoint/FollowUp_AgeAccel_Response_", clock, ".pdf"), Corr_plot, units = "cm", width = 7, height = 10, dpi = 400)
  save(Corr_plot, file = paste0("output/Figures/SingleTimePoint/FollowUp_AgeAccel_Response_", clock, ".RData"))
}


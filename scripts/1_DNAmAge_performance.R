
################################################################################
# SETUP
################################################################################

work_dir <- "S://Project/WP-epigenetics/06_DNAmAge/"
beta_path <- "S:/Project/WP-epigenetics/04_Pipeline_combinedPhases/output/RData/beta_final_autosomes.RData"
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
library(dnaMethyAge)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(writexl)
library(tictoc)
library(car)
library(ggcorrplot)

# LOAD BETA VALUES
load(beta_path)

# LOAD SAMPLESHEET
samplesheet <- as.data.frame(read_xlsx(samples_path))

# LOAD
load("Resources/theme.RData") #ggplot theme
load("Resources/theme_transparent.Rdata") #ggplot theme
load("Resources/methyage_function.RData") #adapted methyage function
load("Resources/preprocessDunedin.RData") # for adapted methyage function

# LOAD RESIDUALISED SNP BETAS
load(SNP_path)

# LOAD CONTROL PROBES
load(ctrl_path)

# ALL CLOCKS (GrimAge2 calculated separately)
clocks <- c("PCHorvathS2013", "PCGrimAge", "PCPhenoAge", "HorvathS2013", "LevineM2018", "DunedinPACE")
clock_names <- c("PC_Horvath", "PC_GrimAge", "PC_PhenoAge", "Horvath", "PhenoAge", "DunedinPACE")
clocks_test <- c("PC_Horvath", "PC_GrimAge", "PC_PhenoAge", "DunedinPACE")
clocks_performance <- c("PC_Horvath", "PC_GrimAge", "PC_PhenoAge")
all_clock_names <- c("PC_Horvath", "PC_GrimAge", "PC_PhenoAge", "Horvath", "PhenoAge", "DunedinPACE", "GrimAge2")

#GrimAGE
gold_path <- "S:/Project/WP-epigenetics/02_Import/DNAmGrimAgeGitHub/input/datMiniAnnotation3_Gold.csv" #by Ake Lu
GrimAge2_path <- "S:/Project/WP-epigenetics/02_Import/DNAmGrimAgeGitHub/input/DNAmGrimAge2_final.Rds"


################################################################################
# PREPARE
################################################################################

# FILTER
table(samplesheet$Time_point)
#D1  D4  M3 
#899 405 404 
length(unique(samplesheet$Resp_nr[which(samplesheet$Diagnosis == "OCD")])) #908
table(samplesheet$Diagnosis)
# CTRL  OCD 
# 384 1708 
beta <- beta[, samplesheet$Basename]

samples_red <- samplesheet[, c("Basename", "Age", "Sex")]
colnames(samples_red) <- c("Sample", "Age", "Sex")


#SET UP DATA FRAMES
clock_performance_all <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 5))
colnames(clock_performance_all) <- c("Clock_Name", "Pearson_corr", "RMSD", "MedianAbsoluteError", "Mean_diff")

DNAmAge_all <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 9))
colnames(DNAmAge_all) <- c("Age", "Sex", "mAge", "Age_Acceleration", "Diagnosis","Time_point", "abs_diff", "diff", "Clock")


# CALCULATE DNAmAGE + PERFORMANCE
for (j in c(1:length(clocks))){
  clock <- clocks[j]
  clock_name <- clock_names[j]
  DNAmAge <- methyAge(beta, clock = clock, inputation = FALSE, simple_mode = TRUE, age_info = samples_red, do_plot = FALSE)
  
  if (clock == "PCGrimAge"){
    DNAmAge <- DNAmAge[, c(1:4, 6)]
  }
  
  DNAmAge$Time_point <- NA
  DNAmAge$Diagnosis <- NA
  
  for (i in c(1:nrow(DNAmAge))){
    Basename <- DNAmAge$Sample[i]
    DNAmAge$Time_point[i] <- samplesheet$Time_point[which(samplesheet$Basename == Basename)]
    DNAmAge$Diagnosis[i] <- samplesheet$Diagnosis[which(samplesheet$Basename == Basename)]
  }
  
  clock_performance <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 5))
  colnames(clock_performance) <- c("Clock_Name", "Pearson_corr", "RMSD", "MeanAbsoluteError", "Mean_diff")
  clock_performance$Clock_Name <- clock_name
  

  DNAmAge_red <- DNAmAge
      
  cor_plot <- ggplot(DNAmAge_red, aes(x= Age, y = mAge)) +
    geom_point(color = "darkgrey", shape = 16) +
    #geom_abline(slope = 1, intercept = 0, color = "darkblue", linewidth = 1) +
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    xlab("Chronological Age") +
    ylab("DNAm Age") +
    ggtitle(clock_name) +
    th + th_transparent +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0("output/Figures/CorrPlots/CorrPlot_", clock_name, ".pdf"), cor_plot, units = "cm", width = 8, height = 10, dpi = 400)
  
  # Pearson correlation
  cor_value <- cor(DNAmAge_red$mAge, DNAmAge_red$Age)
  clock_performance$Pearson_corr[1] <- cor_value
    
  # Root mean square deviation
  RMSD <- sqrt(mean((DNAmAge_red$Age - DNAmAge_red$mAge)^2))
  clock_performance$RMSD[1] <- RMSD
    
  # Mean Absolute Error
  DNAmAge_red$abs_diff <- abs(DNAmAge_red$mAge - DNAmAge_red$Age)
  mae <- mean(DNAmAge_red$abs_diff)
  clock_performance$MeanAbsoluteError[1] <- mae
    
  # Mean Difference
  DNAmAge_red$diff <- DNAmAge_red$mAge - DNAmAge_red$Age
  mean_diff <- mean(DNAmAge_red$diff)
  clock_performance$Mean_diff[1] <- mean_diff
    
  DNAmAge_red$Clock <- clock_name
  DNAmAge_all <- rbind(DNAmAge_all, DNAmAge_red)
  
  clock_performance_all <- rbind(clock_performance_all, clock_performance)
}

# Age acceleration is calculated as the residual resulting from a linear regression model which DNAm age is regressed on chronological age






### GRIMAGE2 ###

#your input/outputOCD file names
grimage2=readRDS(GrimAge2_path)
gold_df <- read.csv(gold_path)

#some set up
Y.pred0.name=c('COX')
Y.pred.name=c('DNAmGrimAge2')
aa.name=c('AgeAccelGrim2')

cpgs=grimage2[[1]]
glmnet.final1=grimage2[[2]]
gold=grimage2[[3]]
cgs <- cpgs$var[grep("^cg", cpgs$var)]

F_scale<-function(INPUT0,Y.pred0.name,Y.pred.name,gold){
  out.para=subset(gold,var=='COX')
  out.para.age=subset(gold,var=='Age')
  m.age=out.para.age$mean
  sd.age=out.para.age$sd
  Y0=INPUT0[,Y.pred0.name]
  Y=(Y0-out.para$mean)/out.para$sd
  INPUT0[,Y.pred.name]=as.numeric((Y*sd.age)+m.age)
  return(INPUT0)
}


samples_red <- samplesheet[, c("Basename", "Age", "Sex")]
samples_red$Sex[which(samples_red$Sex == "F")] <- 1
samples_red$Sex[which(samples_red$Sex == "M")] <- 0
colnames(samples_red) <- c("SampleID", "Age", "Female")

beta_red <- beta[rownames(beta) %in% cpgs$var, samples_red$SampleID]


missing_cgs <- setdiff(cgs, rownames(beta_red)) 
missing_cgs_df <- as.data.frame(matrix(data = NA, nrow = length(missing_cgs), ncol = ncol(beta_red)))
colnames(missing_cgs_df) <- colnames(beta_red)
row.names(missing_cgs_df) <- missing_cgs

for (j in c(1:length(missing_cgs))){
  cg <- missing_cgs[j]
  value <- gold_df$gold[which(gold_df$CpG == cg)]
  missing_cgs_df[j, ] <- value
}

beta_comb <- rbind(beta_red, missing_cgs_df)
beta_add <- as.data.frame(t(beta_comb))
beta_add$Basename <- rownames(beta_add)
dat.meth <- merge(samples_red, beta_add, by.x = "SampleID", by.y = "Basename")
dat.meth$Intercept=1


#step1:generate DNAm Protien single variables

Ys=unique(cpgs$Y.pred)
Ys
for(k in 1:length(Ys)){
  cpgs1=subset(cpgs,Y.pred==Ys[k])
  Xs=subset(dat.meth,select=cpgs1$var)
  Y.pred=as.numeric(as.matrix(Xs)%*%cpgs1$beta)
  dat.meth[,Ys[k]]=Y.pred
  attr(dat.meth[,Ys[k]],"dimnames")<-NULL
}

#step2: generate a raw variable of DNAmGrimAge2(variable name=COX) and calibrate COX in units of years

vars=c('SampleID','Age','Female',Ys)
GrimAge2=subset(dat.meth,select=vars)
GrimAge2$Female <- as.numeric(GrimAge2$Female)
#

GrimAge2$COX=as.numeric(as.matrix(subset(GrimAge2,select=glmnet.final1$var))%*%glmnet.final1$beta)
GrimAge2=F_scale(GrimAge2,'COX',Y.pred.name,gold)
GrimAge2$DNAmtemp<-GrimAge2[,Y.pred.name]
GrimAge2[,aa.name]=residuals(lm(DNAmtemp~Age,data=GrimAge2,na.action = na.exclude))
GrimAge2$DNAmtemp<-GrimAge2$COX<-NULL
#
#rename

old.name=c('DNAmadm','DNAmCystatin_C','DNAmGDF_15','DNAmleptin','DNAmpai_1','DNAmTIMP_1','DNAmlog.CRP','DNAmlog.A1C')
new.name=c('DNAmADM','DNAmCystatinC','DNAmGDF15'  ,'DNAmLeptin','DNAmPAI1'  ,'DNAmTIMP1','DNAmlogCRP','DNAmlogA1C')           
for(k in 1:length(old.name)){
  id=which(names(GrimAge2)==old.name[k])
  names(GrimAge2)[id]=new.name[k]
}



GrimAge2$Time_point <- NA
GrimAge2$Diagnosis <- NA

for (i in c(1:nrow(GrimAge2))){
  Basename <- GrimAge2$SampleID[i]
  GrimAge2$Time_point[i] <- samplesheet$Time_point[which(samplesheet$Basename == Basename)]
  GrimAge2$Diagnosis[i] <- samplesheet$Diagnosis[which(samplesheet$Basename == Basename)]
}


clock_performance_GrimAge2 <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 5))
colnames(clock_performance_GrimAge2) <- c("Clock_Name", "Pearson_corr", "RMSD", "MeanAbsoluteError", "Mean_diff")
clock_performance_GrimAge2$Clock_Name <- "GrimAge2"


GrimAge2_red <- GrimAge2

cor_plot <- ggplot(GrimAge2_red, aes(x= Age, y = DNAmGrimAge2)) +
  geom_point(color = "darkgrey", shape = 16) +
  #geom_abline(slope = 1, intercept = 0, color = "darkblue", linewidth = 1) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  xlab("Chronological Age") +
  ylab("DNAm Age") +
  ggtitle(clock_name) +
  th + th_transparent +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("output/Figures/CorrPlots/CorrPlot_GrimAge2.pdf", cor_plot, units = "cm", width = 8, height = 10, dpi = 400)

# Pearson correlation
cor_value <- cor(GrimAge2_red$DNAmGrimAge2, GrimAge2_red$Age)
clock_performance_GrimAge2$Pearson_corr[1] <- cor_value

# Root mean square deviation
RMSD <- sqrt(mean((GrimAge2_red$Age - GrimAge2_red$DNAmGrimAge2)^2))
clock_performance_GrimAge2$RMSD[1] <- RMSD

# Mean Absolute Error
GrimAge2_red$abs_diff <- abs(GrimAge2_red$DNAmGrimAge2 - GrimAge2_red$Age)
mae <- mean(GrimAge2_red$abs_diff)
clock_performance_GrimAge2$MeanAbsoluteError[1] <- mae

# Mean Difference
GrimAge2_red$diff <- GrimAge2_red$DNAmGrimAge2 - GrimAge2_red$Age
mean_diff <- mean(GrimAge2_red$diff)
clock_performance_GrimAge2$Mean_diff[1] <- mean_diff

write_xlsx(GrimAge2_red, path = "output/Tables/GrimAge2.xlsx")

GrimAge2_red2 <- GrimAge2_red[, c(1:3, 14:19)]
GrimAge2_red2$Clock <- "GrimAge2"
colnames(GrimAge2_red2) <- colnames(DNAmAge_all)
GrimAge2_red2$Sex[which(GrimAge2_red2$Sex == 1)] <- "F"
GrimAge2_red2$Sex[which(GrimAge2_red2$Sex == 0)] <- "M"

DNAmAge_all <- rbind(DNAmAge_all, GrimAge2_red2)

clock_performance_all <- rbind(clock_performance_all, clock_performance_GrimAge2)



# EXPORT COMBINED OVERVIEWS
write_xlsx(clock_performance_all, path = "output/Tables/Clock_performances.xlsx")
write_xlsx(DNAmAge_all, path = "output/Tables/DNAmAge_All.xlsx")




################################################################################
# PERFORMANCE PLOTS
################################################################################

clock_performance_all <- clock_performance_all[which(clock_performance_all$Clock_Name %in% clocks_performance), ]

# PEARSON CORRELATION
clock_performance_all <- clock_performance_all[order(clock_performance_all$Pearson_corr, decreasing = TRUE), ]
clock_performance_all$Clock_Name <- factor(clock_performance_all$Clock_Name, levels = clock_performance_all$Clock_Name)

Pearson_plot_red <- ggplot(clock_performance_all, aes(x= Clock_Name, y = Pearson_corr)) +
  geom_bar(stat = "identity", position = "dodge", fill = "#00524b") +
  th + th_transparent +
  geom_hline(yintercept = 0, color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(legend.position = "top", legend.title = element_blank(), legend.box.background = element_blank(),
        panel.grid.major.y = element_line(color = "lightgrey")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1), breaks = seq(0,1,0.2)) +
  labs(x = "clock name", y = "Pearson r (predicted ~ chronol. age)")
Pearson_plot_red
ggsave("output/Figures/Performance/Pearson_all.pdf", Pearson_plot_red, units = "cm", width = 6, height = 12, dpi = 400)



# MeanAbsoluteError
clock_performance_all <- clock_performance_all[order(clock_performance_all$MeanAbsoluteError, decreasing = FALSE), ]
clock_performance_all$Clock_Name <- factor(clock_performance_all$Clock_Name, levels = clock_performance_all$Clock_Name)

MAE_plot_red <- ggplot(clock_performance_all, aes(x= Clock_Name, y = MeanAbsoluteError)) +
  geom_bar(stat = "identity", position = "dodge", fill = "#00524b") +
  th + th_transparent +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.major.y = element_line(color = "lightgrey")) +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "top", legend.title = element_blank(), legend.box.background = element_blank()) +
  labs(x = "clock name", y = "mean absolute error (in years)")
MAE_plot_red
ggsave("output/Figures/Performance/MeanAbsoluteError_all.pdf", MAE_plot_red, units = "cm", width = 6, height = 12, dpi = 400)


# Mean Error
clock_performance_all <- clock_performance_all[order(clock_performance_all$Mean_diff, decreasing = FALSE), ]
clock_performance_all$Clock_Name <- factor(clock_performance_all$Clock_Name, levels = clock_performance_all$Clock_Name)

ME_plot_red <- ggplot(clock_performance_all, aes(x= Clock_Name, y = Mean_diff)) +
  geom_bar(stat = "identity", position = "dodge", fill = "#00524b") +
  th + th_transparent +
  geom_hline(yintercept = 0, color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.major.y = element_line(color = "lightgrey")) +
  theme(legend.position = "top", legend.title = element_blank(), legend.box.background = element_blank()) +
  labs(x = "clock name", y = "mean error (in years)")
ME_plot_red

ggsave("output/Figures/Performance/Mean_error_all.pdf", ME_plot_red, units = "cm", width = 6, height = 12, dpi = 400)



# Root mean square deviation
clock_performance_all <- clock_performance_all[order(clock_performance_all$RMSD, decreasing = FALSE), ]
clock_performance_all$Clock_Name <- factor(clock_performance_all$Clock_Name, levels = clock_performance_all$Clock_Name)

RMSD_plot_red <- ggplot(clock_performance_all, aes(x= Clock_Name, y = RMSD, fill = Sample)) +
  geom_bar(stat = "identity", position = "dodge", fill = "#00524b") +
  th + th_transparent +
  geom_hline(yintercept = 0, color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "top", legend.title = element_blank(), legend.box.background = element_blank(),
        panel.grid.major.y = element_line(color = "lightgrey")) +
  labs(x = "clock name", y = "RMSD (in years)")
RMSD_plot_red

ggsave("output/Figures/Performance/RMSD_all.pdf", RMSD_plot_red, units = "cm", width = 6, height = 12, dpi = 400)


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

# new column
DNAmAge_all$AgeAccelResid <- NA


# residualize and scale DNAmAge for each clock -> AgeAccelResid
for (clock in all_clock_names){
  
  mAge <- DNAmAge_all$mAge[DNAmAge_all$Clock == clock]
  
  # residualize
  model <- lm(mAge ~ Age + Sex + cell_type1 + cell_type2 + AMP_Plate + 
                ctrl_pc1 + ctrl_pc2 + ctrl_pc3 + ctrl_pc5 + ctrl_pc6 + ctrl_pc7 + ctrl_pc8 + ctrl_pc9 + ctrl_pc10 + anc_pc1)
  AgeAccelResidvalues <- residuals(model)
  
  # add to existing df
  DNAmAge_all$AgeAccelResid[which(DNAmAge_all$Clock == clock)] <- AgeAccelResidvalues
}


write_xlsx(DNAmAge_all, path = "output/Tables/DNAmAge_All.xlsx")





################################################################################
# CORRELATION BETWEEN CLOCKS
################################################################################

DNAmAge_all_test <- DNAmAge_all[which(DNAmAge_all$Clock %in% clocks_test), ]

# ACCELERATION

# transform
DNAmAge_red <- DNAmAge_all_test[, c("Sample", "Age_Acceleration", "Clock")]
DNAmAge_red$Resp_nr <- as.factor(DNAmAge_red$Sample)

DNAmAge_wide <- as.data.frame(DNAmAge_red %>%
                                pivot_wider(
                                  id_cols = Sample, 
                                  names_from = Clock,
                                  values_from = Age_Acceleration 
                                )
)

# correlation matrix
cor_matrix <- cor(
  DNAmAge_wide[,-1],
  use = "pairwise.complete.obs",
  method = "pearson"
)

# heatmap
corr_heatmap <- ggcorrplot(cor_matrix, method = "square", type = "upper", lab = TRUE, lab_size = 3, tl.cex = 10, 
                           colors = c("darkblue", "white", "darkred"), 
                           title = "Age Acceleration",
                           ggtheme = th,
                           hc.order = TRUE
                           ) +
  theme(plot.title = element_text(hjust = 0.5))
corr_heatmap
ggsave(plot = corr_heatmap, filename = "output/Figures/Correlation_heatmap_Acceleration.pdf")


# DNAmAGE

# transform
DNAmAge_red <- DNAmAge_all_test[, c("Sample", "mAge", "Clock")]
DNAmAge_red$Resp_nr <- as.factor(DNAmAge_red$Sample)

DNAmAge_wide <- as.data.frame(DNAmAge_red %>%
                                pivot_wider(
                                  id_cols = Sample, 
                                  names_from = Clock,
                                  values_from = mAge 
                                )
)

# correlation matrix
cor_matrix <- cor(
  DNAmAge_wide[,-1],
  use = "pairwise.complete.obs",
  method = "pearson"
)

# heatmap
corr_heatmap <- ggcorrplot(cor_matrix, method = "square", type = "upper", lab = TRUE, lab_size = 3, tl.cex = 10, 
                           colors = c("darkblue", "white", "darkred"), 
                           title = "Epigenetic Age",
                           ggtheme = th,
                           hc.order = TRUE
) +
  theme(plot.title = element_text(hjust = 0.5))
corr_heatmap
ggsave(plot = corr_heatmap, filename = "output/Figures/Correlation_heatmap_DNAmAge.pdf")



# AgeAccelResid

# transform
DNAmAge_red <- DNAmAge_all_test[, c("Sample", "AgeAccelResid", "Clock")]
DNAmAge_red$Resp_nr <- as.factor(DNAmAge_red$Sample)

DNAmAge_wide <- as.data.frame(DNAmAge_red %>%
                                pivot_wider(
                                  id_cols = Sample, 
                                  names_from = Clock,
                                  values_from = AgeAccelResid 
                                )
)

# correlation matrix
cor_matrix <- cor(
  DNAmAge_wide[,-1],
  use = "pairwise.complete.obs",
  method = "pearson"
)

# heatmap
corr_heatmap <- ggcorrplot(cor_matrix, method = "square", type = "upper", lab = TRUE, lab_size = 3, tl.cex = 10, 
                           colors = c("darkblue", "white", "darkred"), 
                           title = "AgeAccelResid",
                           ggtheme = th,
                           hc.order = TRUE
) +
  theme(plot.title = element_text(hjust = 0.5))
corr_heatmap
ggsave(plot = corr_heatmap, filename = "output/Figures/Correlation_heatmap_AgeAccelResid.pdf")



################################################################################
# NR OF INDEPENDENT TESTS (Li/Ji method)
################################################################################

### CASE CONTROL ###

CaseCtrl_samples <- samplesheet$Basename[which(samplesheet$Diagnosis == "CTRL" | samplesheet$Time_point == "D1")]

DNAmAge_red <- DNAmAge_all[which(DNAmAge_all$Clock %in% clocks_test & DNAmAge_all$Sample %in% CaseCtrl_samples),
                           c("Sample", "AgeAccelResid", "Clock")]

DNAmAge_wide <- as.data.frame(DNAmAge_red %>%
                                pivot_wider(
                                  id_cols = Sample, 
                                  names_from = Clock,
                                  values_from = AgeAccelResid
                                )
)


# correlation matrix
cor_matrix <- cor(
  DNAmAge_wide[,-1],
  use = "pairwise.complete.obs",
  method = "pearson"
)



eigen_values <- eigen(cor_matrix, symmetric = TRUE)$values
Meff <- sum(eigen_values/ (eigen_values + 1))
print(Meff) #1.769173
save(Meff, file = "output/RData/Meff_CaseCtrl.RData")


### OCD ONLY ###

OCD_samples <- samplesheet$Basename[which(samplesheet$Diagnosis == "OCD")]

# transform
DNAmAge_red <- DNAmAge_all[which(DNAmAge_all$Clock %in% clocks_test & DNAmAge_all$Sample %in% OCD_samples), 
                           c("Sample", "AgeAccelResid", "Time_point", "Clock")]
DNAmAge_red$Resp_nr <- NA
for (i in 1:nrow(DNAmAge_red)){
  Sample <- DNAmAge_red$Sample[i]
  Resp_nr <- samplesheet$Resp_nr[which(samplesheet$Basename == Sample)]
  DNAmAge_red$Resp_nr[i] <- Resp_nr
}
DNAmAge_red <- DNAmAge_red[, c("Resp_nr", "AgeAccelResid", "Time_point", "Clock")]
DNAmAge_red$Resp_nr <- as.factor(DNAmAge_red$Resp_nr)

DNAmAge_wide <- as.data.frame(DNAmAge_red %>%
                                mutate(TimeClock = paste0(Time_point, "_", Clock)) %>%
                                pivot_wider(
                                  id_cols = Resp_nr, 
                                  names_from = TimeClock,
                                  values_from = AgeAccelResid
                                )
)

# correlation matrix
cor_matrix <- cor(
  DNAmAge_wide[,-1],
  use = "pairwise.complete.obs",
  method = "pearson"
)



eigen_values <- eigen(cor_matrix, symmetric = TRUE)$values
Meff <- sum(eigen_values/ (eigen_values + 1))
print(Meff) #3.434289
save(Meff, file = "output/RData/Meff_OCD.RData")




### CASE FOLLOW-UP CONTROL ###

CaseCtrl_samples <- samplesheet$Basename[which(samplesheet$Diagnosis == "CTRL" | samplesheet$Time_point == "M3")]

DNAmAge_red <- DNAmAge_all[which(DNAmAge_all$Clock %in% clocks_test & DNAmAge_all$Sample %in% CaseCtrl_samples),
                           c("Sample", "AgeAccelResid", "Clock")]

DNAmAge_wide <- as.data.frame(DNAmAge_red %>%
                                pivot_wider(
                                  id_cols = Sample, 
                                  names_from = Clock,
                                  values_from = AgeAccelResid
                                )
)


# correlation matrix
cor_matrix <- cor(
  DNAmAge_wide[,-1],
  use = "pairwise.complete.obs",
  method = "pearson"
)



eigen_values <- eigen(cor_matrix, symmetric = TRUE)$values
Meff <- sum(eigen_values/ (eigen_values + 1))
print(Meff) #1.755454
save(Meff, file = "output/RData/Meff_CaseFollowUpCtrl.RData")




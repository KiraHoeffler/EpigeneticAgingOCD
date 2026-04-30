
################################################################################
# SET-UP
################################################################################

# SET WORKING DIRECTORY
setwd("S:/Project/WP-epigenetics/06_DNAmAge/")

# LOAD PACKAGES
library(readxl)
library(writexl)
library(ggplot2)
library(tidyverse)

# LOAD THEME
load("Resources/theme.RData")

################################################################################
# SAMPLESHEET
################################################################################

# IMPORT SAMPLESHEET
samplesheet_D1 <- as.data.frame(read_xlsx("output/Tables/Samplesheet_CaseCtrl.xlsx")) #baseline
samplesheet_M3 <- as.data.frame(read_xlsx("output/Tables/Samplesheet_CaseFollowUpCtrl.xlsx")) #follow-up
samplesheet_OCD <- as.data.frame(read_xlsx("output/Tables/Samplesheet_OCD.xlsx")) #all OCD

################################################################################
# DEMOGRAPHICS
################################################################################

### SEX ###

samplesheet_resp <- unique(samplesheet_OCD[, c("Indiv_ID", "Response_status", "Percent_improvement", "Sex")])

table(samplesheet_resp$Sex)
# F   M 
# 630 278

table(samplesheet_D1$Diagnosis, samplesheet_D1$Sex)
# F   M
# CTRL 275 109
# OCD  281 106
table(samplesheet_M3$Diagnosis, samplesheet_M3$Sex)
# F   M
# CTRL 275 109
# OCD  274 106



### RESPONSE ###
mean(samplesheet_resp$Percent_improvement, na.rm = TRUE) #60.48
sd(samplesheet_resp$Percent_improvement, na.rm = TRUE) #27.25

table(is.na(samplesheet_resp$Response_status))
#FALSE  TRUE 
# 889    19 


### TIME POINT ###

table(samplesheet_OCD$Time_point)
# D1  D4  M3 
# 899 405 404 


### ANCESTRY ###

samplesheet_OCD_anc <- unique(samplesheet_OCD[, c("Resp_nr","Ancestry")])
table(samplesheet_OCD_anc$Ancestry)
table(is.na(samplesheet_OCD_anc$Ancestry))
# AFR     AMR     EAS     EUR     SAS unknown 
# 1       1       5     378       4      16 
# FALSE  TRUE 
# 405   503 

samplesheet_CTRL <- samplesheet_D1[which(samplesheet_D1$Diagnosis == "CTRL"), ]
table(samplesheet_CTRL$Ancestry)
#EUR     SAS unknown 
#132       1       3

table(is.na(samplesheet_CTRL$Ancestry))
#FALSE  TRUE 
#136   248 


### COMORBIDITIES AND MEDICATION ###

samplesheet_OCD_med_com <- unique(samplesheet_OCD[, c("Resp_nr","Sex", "Psychoactive_medicine", "Antidepressants", "Stimulants", "Antipsychotics", 
                                               "Mood_stabilizers",  "Sedatives", "Anxiolytics", "AnyComorbidity_CURRENT", "Depression_combined_CURRENT", 
                                               "ManicEpisode_CURRENT", "HypomanicEpisode_CURRENT", "Bipolar1_CURRENT", "Bipolar2_CURRENT",
                                               "BipolarMania_CURRENT", "Bipolar_combined_CURRENT", "PanicDisorder_CURRENTMONTH", "Agoraphobia_CURRENT", 
                                               "SocialPhobia_combined_CURRENTMONTH", "PTSD_CURRENTMONTH", "AlcoholDependence_12MONTH",
                                               "DrugDependence_12MONTH", "Psychosis_CURRENT", "PsychosisDuringMood_CURRENT", 
                                               "Anorexia_CURRENT3MONTHS", "Bulimia_CURRENT3MONTHS", "GAD_CURRENT6MONTHS", 
                                               "AntiSocialPersonalityDisorder_LIFE")])

length(unique(samplesheet_OCD_med_com$Resp_nr)) #908


### MEDICINE ###

table(samplesheet_OCD_med_com$Psychoactive_medicine)
# N   Y 
# 571 270 

table(is.na(samplesheet_OCD_med_com$Psychoactive_medicine))
# FALSE  TRUE 
# 841    67 

table(samplesheet_OCD_med_com$Antidepressants)
# N   Y 
# 618 225 

table(samplesheet_OCD_med_com$Stimulants)
# N   Y 
# 822  19

table(samplesheet_OCD_med_com$Antipsychotics)
# N   Y 
# 794  47 

table(samplesheet_OCD_med_com$Mood_stabilizers)
#N   Y 
#822  19

table(samplesheet_OCD_med_com$Sedatives)
#N   Y 
#831   10

table(samplesheet_OCD_med_com$Anxiolytics)
#N   Y 
#832   9




### COMORBIDITIES ###

table(samplesheet_OCD_med_com$AnyComorbidity_CURRENT)
#  N   Y 
# 521 336 

table(is.na(samplesheet_OCD_med_com$AnyComorbidity_CURRENT))
#FALSE  TRUE 
#857    51 

table(samplesheet_OCD_med_com$Depression_combined_CURRENT)
#N   Y 
#742 115

table(samplesheet_OCD_med_com$ManicEpisode_CURRENT)
#N 
#857

table(samplesheet_OCD_med_com$HypomanicEpisode_CURRENT)
#N 
#857

table(samplesheet_OCD_med_com$Bipolar_combined_CURRENT)
#N   Y 
#852   5

table(samplesheet_OCD_med_com$PanicDisorder_CURRENTMONTH)
#N   Y 
#788  69

table(samplesheet_OCD_med_com$Agoraphobia_CURRENT)
#N   Y 
#810  47

table(samplesheet_OCD_med_com$SocialPhobia_combined_CURRENTMONTH)
#  N   Y 
# 775  82

table(samplesheet_OCD_med_com$PTSD_CURRENTMONTH)
#N   Y 
#832  25 

table(samplesheet_OCD_med_com$AlcoholDependence_12MONTH)
#N 
#857

table(samplesheet_OCD_med_com$Psychosis_CURRENT)
#N 
#857

table(samplesheet_OCD_med_com$PsychosisDuringMood_CURRENT)
#N 
#857

table(samplesheet_OCD_med_com$Anorexia_CURRENT3MONTHS)
#N   Y 
#845  12

table(samplesheet_OCD_med_com$Bulimia_CURRENT3MONTHS)
#N   Y 
#852   5 

table(samplesheet_OCD_med_com$GAD_CURRENT6MONTHS)
#N   Y 
#719 138




### AGE (AND CASE CONTROL DIFFERENCES) ###

mean(samplesheet_OCD$Age)
# 30.83071
sd(samplesheet_OCD$Age)
#9.819882


mean(samplesheet_D1$Age[which(samplesheet_D1$Diagnosis == "CTRL")])
# 30.23997
sd(samplesheet_D1$Age[which(samplesheet_D1$Diagnosis == "CTRL")])
#6.812872

mean(samplesheet_D1$Age[which(samplesheet_D1$Diagnosis == "OCD")])
# 30.94264
sd(samplesheet_D1$Age[which(samplesheet_D1$Diagnosis == "OCD")])
#9.494401

mean(samplesheet_M3$Age[which(samplesheet_M3$Diagnosis == "OCD")])
#31.10395
sd(samplesheet_M3$Age[which(samplesheet_M3$Diagnosis == "OCD")])
#9.530886





wilcox.test(Age ~ Diagnosis, data = samplesheet_D1)
#data:  Age by Diagnosis
#W = 76306, p-value = 0.5174
#alternative hypothesis: true location shift is not equal to 0

chisq.test(table(samplesheet_D1$Diagnosis, samplesheet_D1$Sex))
# data:  table(samplesheet_D1$Diagnosis, samplesheet_D1$Sex)
# X-squared = 0.051896, df = 1, p-value = 0.8198


################################################################################
# TEST ASSOCIATION
################################################################################
# Kendells Tau because 1) Percent_improvement not normally distributed and 2) handles ties better than spearman

# Baseline severity
samplesheet_test <- unique(samplesheet_OCD[, c("Resp_nr", "Pre_YB_Sum", "Percent_improvement", "Age_Diagnosis")])
mean(samplesheet_test$Pre_YB_Sum, na.rm = TRUE) # 25.82359
sd(samplesheet_test$Pre_YB_Sum, na.rm = TRUE) #4.228278
cor.test(samplesheet_test$Pre_YB_Sum, samplesheet_test$Percent_improvement, method = "kendall")
# data:  samplesheet_test$Pre_YB_Sum and samplesheet_test$Percent_improvement
# z = 1.7749, p-value = 0.07592
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# 0.04119757 

# Age of onset
mean(samplesheet_OCD$Age_Diagnosis, na.rm = TRUE) # 21.74346
sd(samplesheet_OCD$Age_Diagnosis, na.rm = TRUE) #9.326973
table(is.na(samplesheet_test$Age_Diagnosis))
# FALSE  TRUE 
# 813    95
samplesheet_test <- unique(samplesheet_OCD[, c("Resp_nr", "Age_Diagnosis", "Percent_improvement")])
cor.test(samplesheet_test$Age_Diagnosis, samplesheet_test$Percent_improvement, method = "kendall")
# data:  samplesheet_test$Age_Diagnosis and samplesheet_test$Percent_improvement
# z = -0.75954, p-value = 0.4475
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.01834229 

# Age
samplesheet_test <- unique(samplesheet_OCD[, c("Resp_nr", "Age", "Percent_improvement")])
samplesheet_test2 <- na.omit(samplesheet_test[, c("Age", "Percent_improvement")])
cor.test(samplesheet_test2$Age, samplesheet_test2$Percent_improvement, method = "kendall")
# data:  samplesheet_test2$Age and samplesheet_test2$Percent_improvement
# z = -3.2809, p-value = 0.001035
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.06073489 
cor(samplesheet_test2$Age, samplesheet_test2$Percent_improvement) #-0.08346529



# Sex
samplesheet_test <- unique(samplesheet_OCD[, c("Resp_nr", "Sex", "Percent_improvement")])
wilcox.test(samplesheet_test$Percent_improvement ~ samplesheet_test$Sex)
# data:  samplesheet_test$Percent_improvement by samplesheet_test$Sex
# W = 91409, p-value = 0.02953
# alternative hypothesis: true location shift is not equal to 0
mean(samplesheet_test$Percent_improvement[which(samplesheet_test$Sex == "M")], na.rm = TRUE) #57.84288
mean(samplesheet_test$Percent_improvement[which(samplesheet_test$Sex == "F")], na.rm = TRUE) #61.63827

### MEDICATION ###

samplesheet_test <- unique(samplesheet_OCD[, c("Resp_nr", "Psychoactive_medicine", "Antidepressants", "Stimulants", "Antipsychotics", 
                                               "Mood_stabilizers",  "Sedatives", "Anxiolytics", "Percent_improvement")])

# psychoactive medicine
wilcox.test(samplesheet_test$Percent_improvement ~ samplesheet_test$Psychoactive_medicine)
# data:  samplesheet_test$Percent_improvement by samplesheet_test$Psychoactive_medicine
# W = 75604, p-value = 0.5572
# alternative hypothesis: true location shift is not equal to 0

# antidepressants
wilcox.test(samplesheet_test$Percent_improvement ~ samplesheet_test$Antidepressants)
# data:  samplesheet_test$Percent_improvement by samplesheet_test$Antidepressants
# W = 66123, p-value = 0.9037
# alternative hypothesis: true location shift is not equal to 0

# Stimulants
wilcox.test(samplesheet_test$Percent_improvement ~ samplesheet_test$Stimulants)
# data:  samplesheet_test$Percent_improvement by samplesheet_test$Stimulants
# W = 8788, p-value = 0.1301
# alternative hypothesis: true location shift is not equal to 0

# Antipsychotics
wilcox.test(samplesheet_test$Percent_improvement ~ samplesheet_test$Antipsychotics)
# data:  samplesheet_test$Percent_improvement by samplesheet_test$Antipsychotics
# W = 19984, p-value = 0.1213
# alternative hypothesis: true location shift is not equal to 0

# Mood_stabilizers
wilcox.test(samplesheet_test$Percent_improvement ~ samplesheet_test$Mood_stabilizers)
# data:  samplesheet_test$Percent_improvement by samplesheet_test$Mood_stabilizers
# W = 9546.5, p-value = 0.02311
# alternative hypothesis: true location shift is not equal to 0
mean(samplesheet_test$Percent_improvement[which(samplesheet_test$Mood_stabilizers == "Y")], na.rm = TRUE) #44.49731
mean(samplesheet_test$Percent_improvement[which(samplesheet_test$Mood_stabilizers == "N")], na.rm = TRUE) #60.94851

# Sedatives
wilcox.test(samplesheet_test$Percent_improvement ~ samplesheet_test$Sedatives)
# data:  samplesheet_test$Percent_improvement by samplesheet_test$Sedatives
# W = 4031, p-value = 0.619
# alternative hypothesis: true location shift is not equal to 0

# Anxiolytics
wilcox.test(samplesheet_test$Percent_improvement ~ samplesheet_test$Anxiolytics)
# data:  samplesheet_test$Percent_improvement by samplesheet_test$Anxiolytics
# W = 4739.5, p-value = 0.1355
# alternative hypothesis: true location shift is not equal to 0



### COMORDBIDITY ###

samplesheet_test <- unique(samplesheet_OCD[, c("Resp_nr", "Percent_improvement", "AnyComorbidity_CURRENT", "Depression_combined_CURRENT", 
                                                                                                 "Bipolar_combined_CURRENT", "PanicDisorder_CURRENTMONTH", "Agoraphobia_CURRENT", 
                                                                                                 "SocialPhobia_combined_CURRENTMONTH", "PTSD_CURRENTMONTH", "DrugDependence_12MONTH", 
                                                                                                 "Anorexia_CURRENT3MONTHS", "Bulimia_CURRENT3MONTHS", "GAD_CURRENT6MONTHS")])

# Comorbidity
wilcox.test(samplesheet_test$Percent_improvement ~ samplesheet_test$AnyComorbidity_CURRENT)
# data:  samplesheet_test$Percent_improvement by samplesheet_test$AnyComorbidity_CURRENT
# W = 88321, p-value = 0.207
# alternative hypothesis: true location shift is not equal to 0

# Depression
wilcox.test(samplesheet_test$Percent_improvement ~ samplesheet_test$Depression_combined_CURRENT)
# data:  samplesheet_test$Percent_improvement by samplesheet_test$Depression_combined_CURRENT
# W = 46766, p-value = 0.03126
# alternative hypothesis: true location shift is not equal to 0
mean(samplesheet_test$Percent_improvement[which(samplesheet_test$Depression_combined_CURRENT == "Y")], na.rm = TRUE) #55.38612
mean(samplesheet_test$Percent_improvement[which(samplesheet_test$Depression_combined_CURRENT == "N")], na.rm = TRUE) #61.11536

# Bipolar
wilcox.test(samplesheet_test$Percent_improvement ~ samplesheet_test$Bipolar_combined_CURRENT)
# data:  samplesheet_test$Percent_improvement by samplesheet_test$Bipolar_combined_CURRENT
# W = 2525, p-value = 0.4127
# alternative hypothesis: true location shift is not equal to 0

# Panic
wilcox.test(samplesheet_test$Percent_improvement ~ samplesheet_test$PanicDisorder_CURRENTMONTH)
# data:  samplesheet_test$Percent_improvement by samplesheet_test$PanicDisorder_CURRENTMONTH
# W = 21884, p-value = 0.02474
# alternative hypothesis: true location shift is not equal to 0
mean(samplesheet_test$Percent_improvement[which(samplesheet_test$PanicDisorder_CURRENTMONTH == "Y")], na.rm = TRUE) #67.07555
mean(samplesheet_test$Percent_improvement[which(samplesheet_test$PanicDisorder_CURRENTMONTH == "N")], na.rm = TRUE) #59.73334

# Agoraphobia
wilcox.test(samplesheet_test$Percent_improvement ~ samplesheet_test$Agoraphobia_CURRENT)
# data:  samplesheet_test$Percent_improvement by samplesheet_test$Agoraphobia_CURRENT
# W = 15718, p-value = 0.1788
# alternative hypothesis: true location shift is not equal to 0

# Social phobia
wilcox.test(samplesheet_test$Percent_improvement ~ samplesheet_test$SocialPhobia_combined_CURRENTMONTH)
# data:  samplesheet_test$Percent_improvement by samplesheet_test$SocialPhobia_combined_CURRENTMONTH
# W = 31023, p-value = 0.611
# alternative hypothesis: true location shift is not equal to 0

# PTSD
wilcox.test(samplesheet_test$Percent_improvement ~ samplesheet_test$PTSD_CURRENTMONTH)
# data:  samplesheet_test$Percent_improvement by samplesheet_test$PTSD_CURRENTMONTH
# W = 7751.5, p-value = 0.1569
# alternative hypothesis: true location shift is not equal to 0

# Anorexia
wilcox.test(samplesheet_test$Percent_improvement ~ samplesheet_test$Anorexia_CURRENT3MONTHS)
# data:  samplesheet_test$Percent_improvement by samplesheet_test$Anorexia_CURRENT3MONTHS
# W = 5025.5, p-value = 0.5502
# alternative hypothesis: true location shift is not equal to 0

# Bulimia
wilcox.test(samplesheet_test$Percent_improvement ~ samplesheet_test$Bulimia_CURRENT3MONTHS)
# data:  samplesheet_test$Percent_improvement by samplesheet_test$Bulimia_CURRENT3MONTHS
# W = 3218.5, p-value = 0.03534
# alternative hypothesis: true location shift is not equal to 0
mean(samplesheet_test$Percent_improvement[which(samplesheet_test$Bulimia_CURRENT3MONTHS == "Y")], na.rm = TRUE) #34.53465
mean(samplesheet_test$Percent_improvement[which(samplesheet_test$Bulimia_CURRENT3MONTHS == "N")], na.rm = TRUE) #60.48396


# GAD
wilcox.test(samplesheet_test$Percent_improvement ~ samplesheet_test$GAD_CURRENT6MONTHS)
# data:  samplesheet_test$Percent_improvement by samplesheet_test$GAD_CURRENT6MONTHS
# W = 53885, p-value = 0.02359
# alternative hypothesis: true location shift is not equal to 0
mean(samplesheet_test$Percent_improvement[which(samplesheet_test$GAD_CURRENT6MONTHS == "Y")], na.rm = TRUE) #55.61882
mean(samplesheet_test$Percent_improvement[which(samplesheet_test$GAD_CURRENT6MONTHS == "N")], na.rm = TRUE) #61.24969


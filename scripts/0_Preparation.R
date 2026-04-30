
################################################################################
# SETUP
###############################################################################

# PATH TO WORKING DIRECTORY:
dir_gen <- "S://Project/WP-epigenetics/06_DNAmAge/"

################################################################################
# SET UP FOLDER STRUCTURE FOR THE output FILES
################################################################################

# SET WORKING DIRECTORY
setwd(dir_gen)

# output FOLDER STRUCTURE
if (!dir.exists("output")) dir.create("output")
if (!dir.exists("output/RData")) dir.create("output/RData")
if (!dir.exists("output/Figures")) dir.create("output/Figures")
if (!dir.exists("output/Tables")) dir.create("output/Tables")

if (!dir.exists("output/Figures/CorrPlots")) dir.create("output/Figures/CorrPlots")
if (!dir.exists("output/Figures/Performance")) dir.create("output/Figures/Performance")
if (!dir.exists("output/Figures/Results")) dir.create("output/Figures/Results")
if (!dir.exists("output/Figures/SingleTimePoint")) dir.create("output/Figures/SingleTimePoint")
if (!dir.exists("output/Figures/Baseline_Phase2")) dir.create("output/Figures/Baseline_Phase2")
if (!dir.exists("output/Figures/TimeOnly")) dir.create("output/Figures/TimeOnly")
if (!dir.exists("output/Figures/ResponseTime")) dir.create("output/Figures/ResponseTime")
if (!dir.exists("output/Figures/CaseCtrl")) dir.create("output/Figures/CaseCtrl")
if (!dir.exists("output/Figures/CaseFollowUpCtrl")) dir.create("output/Figures/CaseFollowUpCtrl")
if (!dir.exists("output/Figures/BMIsmoking")) dir.create("output/Figures/BMIsmoking")

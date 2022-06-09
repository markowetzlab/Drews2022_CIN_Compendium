## Script to add noise to feature data and compare signature activities

args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(YAPSA)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(patchwork)
library(Cairo)
library(mclust)
# theme_set(theme_tufte(base_size = 12, base_family = "Arial"))
theme_set(theme_light(base_size = 12, base_family = "Arial"))

BASE="~/cnsigs2_revisions"
# BASE="/Users/drews01/data/phd/cnsigs2/cnsigs2_revisions"
OUTPATH=file.path(BASE, "MC_simulation_overfitting")
SAVEINTERMEDIATEFILES=TRUE
# For runs on the cluster
# NAMEID=paste0("_fullTCGA_1000sims_10pGaussian_10pSamplePoisson_run", args[1])
# For runs on the mac
NAMEID="_fullTCGA_1000sims_10pGaussian_10pSamplePoisson"


## Step 1
ORIECNF=file.path(BASE, "out/2_TCGA_PCAWG_ECNF.rds")
# From cluster
# ORIECNF=file.path(BASE, "../data/3_Pancancer_Signatures/1_tcga_filtered_ecnf.rds")
NUMSIMS=1e3
## Step 2
INPUTMODELS=file.path(BASE, "CINsignatures/Mixmodels_merged_components.rds")
UNINFPRIOR=TRUE
## Step 3
# SIGNATUREFILE=file.path(BASE, "CINsignatures/Signature_Compendium_v5_Cosine-0.74_Signatures.rds")
SIGNATUREFILE=file.path(BASE, "CINsignatures/Signature_Compendium_v5_Cosine-0.74_Signatures_NAMESAPRIL21.rds")

## step 4
# ACTIVITIES=file.path(BASE, "CINsignatures/Signature_Compendium_v5_Cosine-0.74_Activities.rds")
ACTIVITIES=file.path(BASE, "CINsignatures/Signature_Compendium_v5_Cosine-0.74_Activities_NAMESAPRIL21.rds")


## Step 0: Load data and functions
source(file.path(BASE, "scripts/Sig_activity_with_noise_functions.R"))
dir.create(OUTPATH, recursive = TRUE, showWarnings = FALSE)
lOriECNF = readRDS(ORIECNF)
allModels = readRDS(INPUTMODELS)
W = readRDS(SIGNATUREFILE)
dtOri = data.table(melt(readRDS(ACTIVITIES)))



## Step 1: Simulate feature distributions with noise
print("Start simulating data...")
lSimulation = addNoiseToFeatures(lOriECNF, allFeatures = c("changepoint", "segsize",
                                                            "bpchrarm", "osCN", "bp10MB"), 
                                 SDPROP = 20, FINALWIDTH = 0.1, RANGECNAS = 0.1, NUMSIMS = NUMSIMS)
# lSimulation = mergeSims(lSimulation1, lSimulation2)
## Will become too large
# if(SAVEINTERMEDIATEFILES) {
#   saveRDS(lSimulation, file.path(OUTPATH, paste0("0_Features_plus_noise_", NAMEID, ".rds")))
# }



## Step 2: Derive SxC matrices
print("Start deriving SxC matrices data...")
lMatrices = deriveSxCMatrices(lSimulation, allModels = allModels, 
                              allFeatures = names(allModels), UNINFPRIOR = TRUE)
if(SAVEINTERMEDIATEFILES) {
  saveRDS(lMatrices, file.path(OUTPATH, paste0("1_SxC_Matrices", NAMEID, ".rds")))
}



## Step 3: Calculate signature activities
print("Start calculating signature activities...")
lSignatures = calculateSignatureActivities(lMatrices, W = W)
if(SAVEINTERMEDIATEFILES) {
  saveRDS(lSignatures, file.path(OUTPATH, paste0("2_Activities", NAMEID, ".rds")))
}



# Step 4: Visualise results
# Plot boxplot for each signature and sample
lPlots = plotAllSigs(lSignatures, dtOri, orderOri = TRUE, DOTSIZE = 0.25, TITLE = "1000 runs full TCGA", 
                     FILEOUT = file.path(OUTPATH, paste0("3_Boxplots_per_sig", NAMEID, ".txt")))

CairoPNG(file.path(OUTPATH, paste0("3_Boxplots_per_sig", NAMEID, ".png")), width = 2560, height = 1440)
wrap_plots(lPlots, ncol = 5, nrow = 4); dev.off()


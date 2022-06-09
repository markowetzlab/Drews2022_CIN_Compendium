# Call 17 TCGA copy number signatures on other cohorts

rm(list=ls(all=TRUE))

library(this.path)
library(data.table)
library(GenomicRanges)
library(foreach)
library(doMC)
library(YAPSA)

## Data
BASE=dirname(this.path())
OUT=file.path(BASE, "input/Activities")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## Global parameters
# CN profiles and their correction
WIGGLE=0.1
IGNOREDELS=FALSE
CORES=6

# Extract features (trailing "/" important because of legacy code)
REMOVEQUIET=TRUE # Remove samples with less than 20 CNAs. Otherwise activities will be quite unreliable.
PREPATH=file.path(BASE, "input/Data_preparation/refgenome/")
RMNORM=TRUE

# Calculate SxC matrix
INPUTMODELS=file.path(BASE, "input/Data_preparation/Mixmodels_merged_components.rds")
UNINFPRIOR="TRUE"

# Call activities
SIGNATUREFILE=file.path(BASE, "input/Data_preparation/Signature_Compendium_v5_Cosine-0.74_Signatures_NAMESAPRIL21.rds")

## Load functions
source(file.path(BASE, "input/Data_preparation/main_functions.R"))
source(file.path(BASE, "EDF_6_Data_preparation_Signature_activities_from_CN_profiles_functions.R"))


## Go over data and produce activities and CxS matrix
## Gold standard ASCAT on SNP6
DAT=file.path(BASE, "input/TCGA_478_Samples_SNP6_GOLD.rds")
OUTFILEADDON=""
mActSNP6GOLD = calcSigsFromCNProfiles(DAT, OUT, OUTFILEADDON, WIGGLE, IGNOREDELS, CORES, REMOVEQUIET, PREPATH, RMNORM, 
                                         INPUTMODELS, UNINFPRIOR, SIGNATUREFILE, SAVEALLFILES=FALSE)

## ASCAT used with standard settings on a range of genomic technologies
PENALTY=c(35, 50, 70, 100, 140)
for(PEN in PENALTY) {
    DAT=file.path(BASE, paste0("input/ASCAT_WES_pilot_ontarget/ASCAT_TCGA_WES_478samples_penalty", PEN, ".rds"))
    OUTFILEADDON=paste0("_WES_ontarget_penalty", PEN)
    mActWESontarget = calcSigsFromCNProfiles(DAT, OUT, OUTFILEADDON, WIGGLE, IGNOREDELS, CORES, REMOVEQUIET, PREPATH, RMNORM, 
                                               INPUTMODELS, UNINFPRIOR, SIGNATUREFILE, SAVEALLFILES=FALSE) 
}
for(PEN in PENALTY) {
    DAT=file.path(BASE, paste0("input/NewASCAT_CELnoNorm/NewASCAT_478TCGA_CELnoNorm_penalty", PEN, ".rds"))
    OUTFILEADDON=paste0("_NewASCAT_CELnoNorm_penalty", PEN)
    mActCELnoNorm = calcSigsFromCNProfiles(DAT, OUT, OUTFILEADDON, WIGGLE, IGNOREDELS, CORES, REMOVEQUIET, PREPATH, RMNORM, 
                                               INPUTMODELS, UNINFPRIOR, SIGNATUREFILE, SAVEALLFILES=FALSE) 
}
for(PEN in PENALTY) {
    DAT=file.path(BASE, paste0("input/PCAWG_WGS_downSNP6/478_PCAWG_WGS_downSNP6_withNormals_penalty", PEN, ".rds"))
    OUTFILEADDON=paste0("_WGSdownSNP6_withNormals_penalty", PEN)
    mActWGSdownSNP6 = calcSigsFromCNProfiles(DAT, OUT, OUTFILEADDON, WIGGLE, IGNOREDELS, CORES, REMOVEQUIET, PREPATH, RMNORM, 
                                           INPUTMODELS, UNINFPRIOR, SIGNATUREFILE, SAVEALLFILES=FALSE) 
}

## ASCAT.sc used for settings with binned reads (e.g. off-target)
PENALTY=c(0.001, 0.01, 0.05)
for(PEN in PENALTY) {
    DAT=file.path(BASE, paste0("input/ASCATsc_WES_off_target/ASCAT.sc_WES_offtarget_30k_alpha", PEN, ".rds"))
    OUTFILEADDON=paste0("_WES_offtarget_30K_alpha", PEN)
    mActWESofftarget = calcSigsFromCNProfiles(DAT, OUT, OUTFILEADDON, WIGGLE, IGNOREDELS, CORES, REMOVEQUIET, PREPATH, RMNORM, 
                                             INPUTMODELS, UNINFPRIOR, SIGNATUREFILE, SAVEALLFILES=FALSE) 
}
for(PEN in PENALTY) {
    DAT=file.path(BASE, paste0("input/ASCATsc_WGS_downsWGS/ASCAT.sc_WGS_downsWGS_alpha", PEN, ".rds"))
    OUTFILEADDON=paste0("_WGS_downsWGS_alpha", PEN)
    mActWGSdownsWGS = calcSigsFromCNProfiles(DAT, OUT, OUTFILEADDON, WIGGLE, IGNOREDELS, CORES, REMOVEQUIET, PREPATH, RMNORM, 
                                              INPUTMODELS, UNINFPRIOR, SIGNATUREFILE, SAVEALLFILES=FALSE) 
}

## Different copy number callers (full WGS data)
METHOD=c("ABSOLUTE", "Battenberg", "Consensus", "Sclust")
for(MET in METHOD) {
    DAT=file.path(BASE, paste0("input/PCAWG_WGS_full_data/PCAWG_", MET, "_CN_profiles.rds"))
    OUTFILEADDON=paste0("_PCAWG_", MET, "_CN_profiles")
    mActPCAWGFULLWGS = calcSigsFromCNProfiles(DAT, OUT, OUTFILEADDON, WIGGLE, IGNOREDELS, CORES, REMOVEQUIET, PREPATH, RMNORM, 
                                             INPUTMODELS, UNINFPRIOR, SIGNATUREFILE, SAVEALLFILES=FALSE) 
}

METHOD=c("CNVkit_Refitted", "CNVkit_Standard", "Sequenca")
for(MET in METHOD) {
    DAT=file.path(BASE, paste0("input/ShallowWGS_Two_callers/TCGA_WES_", MET, ".rds"))
    OUTFILEADDON=paste0("_TCGA_WES_", MET)
    mActTCGAsWGS = calcSigsFromCNProfiles(DAT, OUT, OUTFILEADDON, WIGGLE, IGNOREDELS, CORES, REMOVEQUIET, PREPATH, RMNORM, 
                                              INPUTMODELS, UNINFPRIOR, SIGNATUREFILE, SAVEALLFILES=FALSE) 
}



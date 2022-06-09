## Compare TCGA signatures with OV signatures
##
## Compare exposures on TCGA OV. Compare signature definitions by mapping TCGA signatures back
## into OV signature feature space. CN feature is omited for that analysis since it doesn't exist
## in TCGA signatures. 
## Metric of comparison is cosine similarity. Significance testing is performed with a permutation
## testing from library (PharmacoGx; a bit overkill to load the whole package. 
## optimise in the future).
## Threshold for plotting is 0.9 to show only the strongest interactions.

rm(list=ls(all=TRUE))

# Standard
library(this.path)
library(data.table)
library(reshape2)
# For plotting
library(ggplot2)
library(ggthemes)
library(lemon)
# For accessing the OV mixture models
library(flexmix)
# For cosine function
library(lsa)
# For cosine permutation testing - not needed as I copied cosinePerm directly into functions file
# library(PharmacoGx)

theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))


## User input
BASE=dirname(this.path())
OUT=file.path(BASE, "output")
OUTTABLES=file.path(BASE, "output")

EXP=file.path(BASE, "input/Signature_Compendium_v5_Cosine-0.74_Activities_THRESH95_NAMESAPRIL21.rds")
FUNCTIONS=file.path(BASE, "SuppFig_25_TCGA_vs_OV_Signatures_Functions.R")

TCGA_MM=file.path(BASE, "input/2_combined_mixmodels_merged_components.rds")
TCGA_SIGS=file.path(BASE, "input/Signature_Compendium_v5_Cosine-0.74_Signatures_NAMESAPRIL21.rds")
TCGA_EXP=file.path(BASE, "input/Signature_Compendium_v5_Cosine-0.74_Activities_THRESH95_NAMESAPRIL21.rds")

OV_MM=file.path(BASE, "input/OV_Copy_number_signatures/OV_Signatures_component_parameters.rds")
OV_SIGS=file.path(BASE, "input/OV_Copy_number_signatures/feat_sig_mat.rds")
OV_EXP=file.path(BASE, "input/OV_Copy_number_signatures/Export-matrix_OV_Sigs_on_TCGA-OV_12112019.rds")

THRESHOLD=0.85
NOCN=TRUE


## Load data
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
exp = readRDS(EXP)

tcgaMM = readRDS(TCGA_MM)
tcgaSigs = readRDS(TCGA_SIGS)
tcgaExp = readRDS(TCGA_EXP)

ovMM = readRDS(OV_MM)
ovSigs = readRDS(OV_SIGS)
ovExp = readRDS(OV_EXP)

source(FUNCTIONS)

## Compare exposures => test and plot (takes a few seconds as it runs 10.000 permutation tests per signature combination)
lExp = compareExposures(EXP1 = ovExp, EXP2 = tcgaExp, NPERM = 1e4, CRIT = "MAX")
dfExp = melt(lExp$expCos)
dfExpCorr = melt(lExp$expCos - lExp$expCosSim)

## Rename OV sigs
dfExp$Var2 = factor(as.character(dfExp$Var2), levels = paste0("s", 1:7), labels = paste0("OV", 1:7))
dfExpCorr$Var2 = factor(as.character(dfExpCorr$Var2), levels = paste0("s", 1:7), labels = paste0("OV", 1:7))

## Plot cosine similarities
pExp = ggplot(dfExp, aes(x = Var1, y = Var2, fill = value)) + geom_tile() + 
  scale_fill_gradient(low = "#ffffe5", high = "#004529", breaks = c(0.1, 0.5, 0.9)) + 
  labs(x = "CIN signatures", y = "OV signatures", fill = "Cosine similarity") + 
  theme(legend.position = "bottom", legend.key.width = unit(0.5,"cm"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), aspect.ratio = 7/17) +
  coord_capped_cart(left = "both", top = "both")
pExp = addSmallLegend(pExp, textSize = 6, spaceLegend = 0.5)

ggsave(file.path(OUT, "SuppFig_25_TCGA_vs_OV_Signatures_Activities.svg"), pExp, width = 90, height = 60, units = "mm")


## Plot cosine similarities corrected by maximum value
pExpCorr = ggplot(dfExpCorr, aes(x = Var1, y = Var2, fill = value)) + geom_tile() + 
  scale_fill_gradient2(low = "#0571b0", mid = "white", high = "#ca0020",
                       breaks = c(-0.25, 0, 0.25)) +
  labs(x=  "CIN signatures", y = "OV signatures", fill = "Corrected cosine similarity") + 
  theme(legend.position = "bottom", legend.key.width = unit(0.5,"cm"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), aspect.ratio = 7/17) +
  coord_capped_cart(left = "both", top = "both")
pExpCorr = addSmallLegend(pExpCorr, textSize = 6, spaceLegend = 0.5)

ggsave(file.path(OUT, "SuppFig_25_TCGA_vs_OV_Signatures_Activities_Corrected.svg"), pExpCorr, width = 90, height = 60, units = "mm")



## Lift over signature definition
## Get conversion matrix
convMatrix = getConversionMatrix(MODELORIGIN = tcgaMM, MODELTARGET = ovMM, plotHeatmap = FALSE, uninfPrior = FALSE, noCN = NOCN)

## Liftover and compare signatures)
tcgaSigsLifted = liftOverSigs(tcgaSigs, convMatrix)

## Compare signatures => test and plot (again takes a few seconds)
lSig = compareSigs(LIFTEDSIGS = tcgaSigsLifted, QUERYSIGS = ovSigs, NPERM = 1e4, CRIT = "MAX")
dfSig = melt(lSig$sigCos)
dfSigCorr = melt(lSig$sigCos - lSig$sigCosSim)

## Rename OV sigs
dfSig$Var2 = factor(as.character(dfSig$Var2), levels = paste0("s", 1:7), labels = paste0("OV", 1:7))
dfSigCorr$Var2 = factor(as.character(dfSigCorr$Var2), levels = paste0("s", 1:7), labels = paste0("OV", 1:7))


## Plot cosine similarities
pSig = ggplot(dfSig, aes(x = Var1, y = Var2, fill = value)) + geom_tile() + 
  scale_fill_gradient(low = "#ffffe5", high = "#004529", breaks = c(0.1, 0.5, 0.9)) + 
  labs(x= "CIN signatures", y = "OV signatures", fill = "Cosine similarity") + 
  theme(legend.position = "bottom", legend.key.width = unit(0.5,"cm"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), aspect.ratio = 7/17) +
  coord_capped_cart(left = "both", top = "both")

pSig = addSmallLegend(pSig, textSize = 6, spaceLegend = 0.5)

ggsave(file.path(OUT, "SuppFig_25_TCGA_vs_OV_Signatures_Definitions.svg"), pSig, width = 90, height = 60, units = "mm")


## Plot cosine similarities corrected by maximum value
pSigCorr = ggplot(dfSigCorr, aes(x = Var1, y = Var2, fill = value)) + geom_tile() + 
  scale_fill_gradientn(colours = c("#0571b0", "white", "#ca0020"), 
                       values = scales::rescale(c(min(dfSigCorr$value), 0, max(dfSigCorr$value)))) + 
  labs(x= "CIN signatures", y = "OV signatures", fill = "Corrected cosine similarity") + 
  theme(legend.position = "bottom", legend.key.width = unit(0.5,"cm"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), aspect.ratio = 7/17) +
  coord_capped_cart(left = "both", top = "both")

pSigCorr = addSmallLegend(pSigCorr, textSize = 6, spaceLegend = 0.5)

ggsave(file.path(OUT, "SuppFig_25_TCGA_vs_OV_Signatures_Definitions_Corrected.svg"), pSigCorr, width = 90, height = 60, units = "mm")


#### Save output for web portal
lOut = list("TCGAOVActivities" = dfExp,
            "TCGAOVActivitiesCorrected" = dfExpCorr,
            "TCGAOVDefinitions" = dfSig,
            "TCGAOVDefinitionsCorrected" = dfSigCorr)
saveRDS(lOut, file.path(OUTTABLES, "Covariates_TCGAOV_sigs.rds"))

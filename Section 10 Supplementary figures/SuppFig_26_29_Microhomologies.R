#### Correlation with microhomologies at SV breakpoints

rm(list=ls(all=TRUE))

library(this.path)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(lemon)
library(patchwork)
library(scales)
library(viridisLite)
library(ggpubr)
library(rstatix)
library(stringr)
library(MASS)

theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5),
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))

## Data paths
BASE=dirname(this.path())
OUT=file.path(BASE, "output")
OUTTABLE=file.path(BASE, "output")


# If you are interested how microhomologies are obtained, please refer to script
# "SuppMat_Summarise_microhomologies_PCAWG_SVs.R" in Section 11 Miscellaneous.
DTHOM=file.path(BASE, "input/PCAWG_Microhomologies_at_SVs.rds")
PCAWGACT=file.path(BASE, "input/PCAWG_signature_activities_THRESH095_NAMESAPRIL21.rds")
META=file.path(BASE, "input/PCAWG_1900_samples_CINSig_activities_metadata_plus_HRDetect.rds")
META2=file.path(BASE, "input/pcawg_donor_clinical_August2016_v9.fixed.tsv")
META3=file.path(BASE, "input/donor.tsv")
META4=file.path(BASE, "input/TCGA_PCAWG_links_plusCancer.rds")
CCS=file.path(BASE, "input/Lundberg_2020_Cell_cycle_score.rds")
CNAS=file.path(BASE, "input/CNAs_PCAWG_1900samples.rds")

dtHom=readRDS(DTHOM)
dtPCAWG = melt(readRDS(PCAWGACT))
dtMeta = readRDS(META)
dtMeta2 = fread(META2)
dtMeta3 = fread(META3)
dtMeta4 = readRDS(META4)

dtCCS = data.table(readRDS(CCS))
dtCNAs = readRDS(CNAS)

#### Part 0: Prepare SVs
## SVs not needed anymore - remove NAs
dtHom = dtHom[ ! is.na(dtHom$HOMLEN),]

## Based on some reviews
dtHom$Pathway = cut(dtHom$HOMLEN, breaks = c(0,2,20,1e4), labels = c("NHEJ", "TMEJ", "SSA"),
                    include.lowest = TRUE, right = FALSE)
## Remove all NAs
dtHom = dtHom[ ! is.na(dtHom$Pathway),]



#### Part I: CNAs vs SVs 
#### => Mitotic catastrophe analysis and figure (SuppFig 14) in script SuppFig_14_Mitotic_Catastrophe.R
#### Here just similar figures with might be interesting to study in the context of microhomologies
## Plot number of CNAs vs number of SVs
dtStat = data.table(table(dtHom$Sample))
colnames(dtStat) = c("Sample", "SVs")

## Load number of CNAs in the PCAWG cohort
dtStat$CNAs = dtCNAs$CNAs[ match(dtStat$Sample, dtCNAs$PCAWG) ]

## Transfer IHR signatures
mPCAWG = readRDS(PCAWGACT)

dtStat$CX2 = mPCAWG[, "CX2"][ match(dtStat$Sample, rownames(mPCAWG)) ]
dtStat$CX3 = mPCAWG[, "CX3"][ match(dtStat$Sample, rownames(mPCAWG)) ]
dtStat$CX5 = mPCAWG[, "CX5"][ match(dtStat$Sample, rownames(mPCAWG)) ]

## Remove all incomplete samples
dtStatFree = dtStat[ ! is.na(dtStat$CX2),]

## Apply robust linear regression to identify samples with heavy CNA load
# Intercept doesn't make sense => remove intercept from model
modelRLM = rlm(CNAs ~ SVs + 0, data = dtStatFree)
dfModel = data.frame(Sample = dtStatFree$Sample, resid = modelRLM$residuals, weight = modelRLM$w)
dfModel = dfModel[order(dfModel$w), ]
cnaheavysamples = dfModel$Sample[ dfModel$weight < 1 & dfModel$resid > 0 ]

dtStatFree$CNAHeavy = factor(dtStatFree$Sample %in% cnaheavysamples, levels = c("FALSE", "TRUE"),
                             labels = c("Linear", "MitoticCatastrophe"))

pA = ggplot(dtStatFree, aes(x = SVs, y = CNAs, colour = CNAHeavy)) + geom_point(alpha = 0.5) +
  coord_capped_cart(left = "both", bottom = "both") + theme(aspect.ratio = 1) +
  geom_abline(intercept = 0, slope = coef(modelRLM)[1]) +
  labs(x = "Number of SVS", y = "Number of CNAs") + theme(legend.position = c(0.8, 0.8)) +
  scale_colour_manual(values = c("MitoticCatastrophe" = "#1f78b4", "Linear" = "#33a02c"))

pB = ggplot(dtStatFree, aes(x = SVs, y = CNAs, colour = CNAHeavy)) + geom_point(alpha = 0.5) +
  coord_capped_cart(left = "both", bottom = "both") + scale_y_log10() +
  scale_x_log10() + labs(x = "Number of SVS", y = "Number of CNAs") +
  theme(legend.position = "none") +
  scale_colour_manual(values = c("MitoticCatastrophe" = "#1f78b4", "Linear" = "#33a02c"))

## Save for supp figure
pCNAs = pA + pB
# ggsave(file.path(OUT, "SuppFig_14_Mitotic_Catastrophe.svg"), pCNAs, width = 180, height = 90, units = "mm")



## Check CNA heavy samples for CX2,3,5
dtStatFree[["CX2"]] = scale(dtStatFree[["CX2"]])
dtStatFree[["CX3"]] = scale(dtStatFree[["CX3"]])
dtStatFree[["CX5"]] = scale(dtStatFree[["CX5"]])
dtStatMelt = melt(dtStatFree, id.vars = c("Sample", "CNAHeavy"), measure.vars = c("CX2","CX5","CX3"))
dtStatMelt$CNAHeavy = factor(dtStatMelt$CNAHeavy, levels = c("Linear", "MitoticCatastrophe"))
pC = ggplot(dtStatMelt, aes(x = variable, y = value, fill = CNAHeavy)) + geom_boxplot() +
  coord_capped_cart(left = "both", bottom = "both", ylim = c(-1, 4)) +
  labs(y = "Scaled signature activity", x = "CIN signature", fill = "CNA heavy samples") +
  scale_fill_manual(values = c("MitoticCatastrophe" = "#1f78b4", "Linear" = "#33a02c"))

# ggsave(file.path(OUT, "SuppFig_14_Mitotic_Catastrophe_Activities.svg"), pC, width = 90, height = 90, units = "mm")



## Test significance
testCNAHeavy1 <- dtStatMelt %>%
  group_by(variable) %>%
  t_test(value ~ CNAHeavy, var.equal = FALSE, p.adjust.method = "BH")
testCNAHeavy2 <- dtStatMelt %>%
  group_by(CNAHeavy) %>%
  t_test(value ~ variable, var.equal = FALSE, p.adjust.method = "BH")





#### Cell cycle score in CNA heavy group - only possible for samples of the TCGA
dtCCS$Sample = substr(dtCCS$Sample, 1, 12)
dtCCS$PCAWG = dtMeta$samplename[ match(dtCCS$Sample, dtMeta$tcga_donor_barcode) ]


dtStatFree$CCS = dtCCS$CCS[ match(dtStatFree$Sample, dtCCS$PCAWG) ]
## According to Lundberg's script on github
dtStatFree$CCS = factor(as.numeric(dtStatFree$CCS), levels = c(1,2,3), labels = c("Low", "Medium", "High"))


## Significant enrichment from low to high for linear to CNA heavy samples
fisher.test(table(dtStatFree$CNAHeavy, dtStatFree$CCS)[,c(1,3)])

## Plot barplot
pCCS = ggplot(dtStatFree[ ! is.na(dtStatFree$CCS) & dtStatFree$CCS != "Medium", ], aes(x = CNAHeavy, fill = CCS)) + geom_bar(position = "fill") +
  scale_fill_manual(values = c("Low" = "black", "Medium" = "grey", "High" = "gold")) + scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "PCWAG / TCGA samples", y = "Proportion of samples", fill = "Cell cycle\nscore (CCS)") +
  coord_capped_cart(left = "both", bottom = "both")

# ggsave(file.path(OUT, "SuppFig_14_Mitotic_Catastrophe_Cell_cycle_score.svg"), pCCS, width = 50/25.4, height = 45/25.4)

# write.table(dtStatFree, file.path(OUTTABLE, "Mitotic_catastrophe_analysis.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



#### Part 2: Microhomology figures
## Convert microhoms to pathways
mPath = unclass(table(dtHom$Sample, dtHom$Pathway))

## Remove empty samples
mPath = mPath[ rowSums(mPath) != 0, ]
LOWERSVTHRESH=TRUE
if(LOWERSVTHRESH) {
  ## 21 is 1st quantile of the number of SV distribution. Just pick this value for now.
  mPath = mPath[ rowSums(mPath) > 21, ]
}


## Normalise or not. That's the question.
NORMALISE=TRUE
if(NORMALISE) {
  mPath = mPath / rowSums(mPath)
}
dtPath = data.table(melt(mPath))
colnames(dtPath) = c("Sample", "Pathway", "Value")

## Initial plots about repair pathways
dtPath$Pathway = factor(dtPath$Pathway, levels = c("NHEJ", "TMEJ", "SSA"))


## Not needed anymore
# pA = ggplot(dtPath, aes(x = Pathway, y = Value)) + geom_boxplot() +
#   coord_capped_cart(bottom = "both", left = "both") + labs(y = "Frequency or proportion")
#
# pB = ggplot(dtPath, aes(x = Pathway, y = Value, group = Sample)) + geom_line(alpha = 0.2) +
#   coord_capped_cart(bottom = "both", left = "both") + labs(y = "Frequency or proportion")



### Test pathways
mPCAWGScaled = scale(mPCAWG)

allPathways = levels(dtPath$Pathway)
allSigs = colnames(mPCAWGScaled)
DOMINANT=1.25
lAll = lapply(allPathways, function(thisPath) {

  dtRepair = dtPath[ dtPath$Pathway == thisPath, ]
  lSigs = lapply(allSigs, function(thisSig) {

    ## Identify dominant samples
    vStatus = cut(mPCAWGScaled[, thisSig], breaks = c(-100, 0, DOMINANT, 100),
        labels = c("Depleted", "Intermediary", "Dominant"))

    ## Transfer status and activity
    dtRepair$Signature = thisSig
    dtRepair$Status = vStatus[ match(dtRepair$Sample, rownames(mPCAWGScaled)) ]
    dtRepair$Activity = mPCAWGScaled[ match(dtRepair$Sample, rownames(mPCAWGScaled)), thisSig ]

    ## Test and plot
    tTest = t.test(dtRepair$Value[ dtRepair$Status == "Depleted" ],
           dtRepair$Value[ dtRepair$Status == "Dominant" ], var.equal = FALSE)
    out = c(thisSig, thisPath, signif(c(tTest$estimate[1], tTest$estimate[2], tTest$estimate[2] - tTest$estimate[1],
            tTest$p.value), 4))

    ## Return results
    lOut = list(dtRepair, out)
    return(lOut)
  })

  ## Simplify results
  lRepair = lapply(lSigs, function(x) x[[1]])
  dtCollectRepair = data.table(do.call(rbind, lRepair))

  lStats = lapply(lSigs, function(x) x[[2]])
  dtStats = data.table(do.call(rbind, lStats))

  colnames(dtStats) = c("Signature", "Pathway", "MeanDepleted", "MeanDominant",
                        "Differnce", "pVal")
  dtStats$MeanDepleted = as.numeric(dtStats$MeanDepleted)
  dtStats$MeanDominant = as.numeric(dtStats$MeanDominant)
  dtStats$Differnce = as.numeric(dtStats$Differnce)
  dtStats$pVal = as.numeric(dtStats$pVal)
  dtStats$pAdj = p.adjust(dtStats$pVal, method = "BH")

  ## Mark samples and signatures as significant
  dtCollectRepair$Significant = factor(dtCollectRepair$Signature %in% dtStats$Signature[ dtStats$pAdj < 0.05 ],
                                       levels = c(TRUE, FALSE), labels = c("Sign", "NS"))
  dtCollectRepair$Plot = paste0(dtCollectRepair$Status, "_", dtCollectRepair$Significant)

  ## Collect for output
  lOut = list(dtCollectRepair, dtStats)
  return(lOut)

})

## Split results again
lRepair = lapply(lAll, function(x) x[[1]])
dtCollectRepair = data.table(do.call(rbind, lRepair))

lStats = lapply(lAll, function(x) x[[2]])
dtStats = data.table(do.call(rbind, lStats))


## Make big plot for supplementary figure
dtRepairFilt = dtCollectRepair[ ! grepl("NA", dtCollectRepair$Plot), ]
dtRepairFilt$Plot = factor(dtRepairFilt$Plot, levels = c("Depleted_Sign", "Intermediary_Sign",
                                                         "Dominant_Sign", "Depleted_NS",
                                                         "Intermediary_NS", "Dominant_NS"))
dtRepairFilt$Signature = factor(dtRepairFilt$Signature, levels = paste0("CX", 1:17))
pAllA = ggplot(dtRepairFilt[ dtRepairFilt$Pathway == "NHEJ", ], aes(x = Signature, y = Value, fill = Plot)) +
  geom_boxplot(outlier.size = 0.1) +
  facet_grid(Pathway ~ ., scales = "free_y") + labs(y = "Proportion of SVs") +
  scale_fill_manual(values = c("Depleted_Sign" = "#edf8b1", "Intermediary_Sign" = "#7fcdbb",
                               "Dominant_Sign" = "#2c7fb8", "Depleted_NS" = "white",
                               "Intermediary_NS" = "white", "Dominant_NS" = "white")) +
  coord_capped_cart(left = "both", bottom = "both", ylim = c(0.35, 0.8)) +
  scale_y_continuous(labels = scales::percent_format())
ggsave(file.path(OUT, "SuppFig_29_A_Microhomologies_all_Signatures.svg"), pAllA,
       width = 180, height = 50, units = "mm")

pAllB = ggplot(dtRepairFilt[ dtRepairFilt$Pathway == "TMEJ", ], aes(x = Signature, y = Value, fill = Plot)) +
  geom_boxplot(outlier.size = 0.1) +
  facet_grid(Pathway ~ ., scales = "free_y") + labs(y = "Proportion of SVs") +
  scale_fill_manual(values = c("Depleted_Sign" = "#edf8b1", "Intermediary_Sign" = "#7fcdbb",
                               "Dominant_Sign" = "#2c7fb8", "Depleted_NS" = "white",
                               "Intermediary_NS" = "white", "Dominant_NS" = "white")) +
  coord_capped_cart(left = "both", bottom = "both", ylim = c(0.05, 0.5)) +
  scale_y_continuous(labels = scales::percent_format())
ggsave(file.path(OUT, "SuppFig_29_B_Microhomologies_all_Signatures.svg"), pAllB,
       width = 180, height = 50, units = "mm")

pAllC = ggplot(dtRepairFilt[ dtRepairFilt$Pathway == "SSA", ], aes(x = Signature, y = Value, fill = Plot)) +
  geom_boxplot(outlier.size = 0.1) +
  facet_grid(Pathway ~ ., scales = "free_y") + labs(y = "Proportion of SVs") +
  scale_fill_manual(values = c("Depleted_Sign" = "#edf8b1", "Intermediary_Sign" = "#7fcdbb",
                               "Dominant_Sign" = "#2c7fb8", "Depleted_NS" = "white",
                               "Intermediary_NS" = "white", "Dominant_NS" = "white")) +
  coord_capped_cart(left = "both", bottom = "both", ylim = c(0, 0.01)) +
  scale_y_continuous(labels = scales::percent_format())
ggsave(file.path(OUT, "SuppFig_29_C_Microhomologies_all_Signatures.svg"), pAllC,
       width = 180, height = 50, units = "mm")


write.table(dtRepairFilt, file.path(OUTTABLE, "Supp_Table_Microhomologies.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
saveRDS(dtRepairFilt, file.path(OUTTABLE, "Supp_Table_Microhomologies.rds"))




#### Part 3: Repair pathways for CNA heavy samples
dfPath = data.frame(melt(mPath))
colnames(dfPath) = c("Sample", "Pathway", "SVs")
dfPath$CNAHeavy = factor(dfPath$Sample %in% cnaheavysamples, levels = c("TRUE", "FALSE"),
                             labels = c("MitoticCatastrophe", "Linear"))

## Dominant signature classification (either CX2, CX5 or CX3)
mPCAWGScaled = scale(mPCAWG)
dfPS = data.frame(mPCAWGScaled[,c("CX2", "CX3", "CX5")])
dfPS$Class = "Neither"
dfPS$Class[ dfPS$CX2 > 1.25 & dfPS$CX2 > dfPS$CX3 & dfPS$CX2 > dfPS$CX5 ] = "CX2"
dfPS$Class[ dfPS$CX3 > 1.25 & dfPS$CX3 > dfPS$CX2 & dfPS$CX3 > dfPS$CX5 ] = "CX3"
dfPS$Class[ dfPS$CX5 > 1.25 & dfPS$CX5 > dfPS$CX2 & dfPS$CX5 > dfPS$CX3 ] = "CX5"

# NA samples are samples without arrays or dCIN
dfPath$SigClass = factor(dfPS$Class[ match(dfPath$Sample, rownames(dfPS)) ],
                         levels = c("Neither", "CX2", "CX5", "CX3"))
dfPathNoNA = dfPath[ ! is.na(dfPath$SigClass), ]
p7a = ggplot(dfPathNoNA, aes(x = Pathway, y = SVs, fill = SigClass)) +
  # geom_jitter(width = 0.1, height = 0, alpha = 0.5) +
  geom_boxplot(outlier.size = 0.25) +
  coord_capped_cart(left = "both", bottom = "both") +
  labs(y = "Frequency of SVs with microhomologies", x = "DNA repair pathway",
       fill = "Dominant IHR\nsignature")
## Zoom for SSA
p7b = ggplot(dfPathNoNA[ dfPathNoNA$Pathway=="SSA", ], aes(x = Pathway, y = SVs, fill = SigClass)) +
  # geom_jitter(width = 0.1, height = 0, alpha = 0.5) +
  geom_boxplot(outlier.size = 0.25) + theme(legend.position = "none") +
  coord_capped_cart(left = "both", bottom = "both", ylim = c(0,0.035)) +
  labs(y = "Frequency of SVs with microhomologies", x = "DNA repair pathway",
       fill = "Dominant IHR\nsignature")
p7 = p7a + (plot_spacer() / p7b)

## This is the plot to show in the supps
cairo_pdf(file.path(OUT, "SuppFigure_26_Microhomologies.pdf"), width = 120/25.4, height = 90/25.4)
print(p7); dev.off()
ggsave(file.path(OUT, "SuppFigure_26_Microhomologies.svg"), p7, width = 120/25.4, height = 90/25.4)

# Test significance
dfPath %>% group_by(Pathway) %>% t_test(SVs ~ SigClass, var.equal = FALSE, p.adjust.method = "BH")

# Neither group is significantly higher than CX2. But this relationship is driven by outliers.
# Therefore we perform a proportion test at zero to test how many patients have zero evidence for SSA and
# which have encountered at least 1 instance of SSA-mediated homology.
# Extract samples and identify samples with zero and non-zero SSA activity:
neither=dfPath[dfPath$Pathway=="SSA" & dfPath$SigClass=="Neither",]
neither$zero = neither$SVs == 0
cx2=dfPath[dfPath$Pathway=="SSA" & dfPath$SigClass=="CX2",]
cx2$zero = cx2$SVs == 0

# CX2-dominant samples have 7% more samples above zero and this is significantly more.
prop.test(matrix(c(table(neither$zero), table(cx2$zero)), nrow = 2, byrow = TRUE))

# CX5-dominant samples have 20% more samples above zero.
cx5=dfPath[dfPath$Pathway=="SSA" & dfPath$SigClass=="CX5",]
cx5$zero = cx5$SVs == 0
prop.test(matrix(c(table(neither$zero), table(cx5$zero)), nrow = 2, byrow = TRUE))

# CX3-dominant samples have 26% more samples above zero.
cx3=dfPath[dfPath$Pathway=="SSA" & dfPath$SigClass=="CX3",]
cx3$zero = cx3$SVs == 0
prop.test(matrix(c(table(neither$zero), table(cx3$zero)), nrow = 2, byrow = TRUE))


#### For the interested, here is a bit of code with microhoms split by cancer type and other covariates
# #### Other ways of carving up the data!
# ## Binarisation
# dfPB = data.frame(mPCAWG[,c("CX2", "CX3", "CX5")] > 0)
# dfPB$Class = NA
# dfPB$Class[ dfPB$CX2 == FALSE & dfPB$CX5 == FALSE & dfPB$CX3 == FALSE ] = "Background"
# dfPB$Class[ dfPB$CX2 == TRUE & dfPB$CX5 == FALSE & dfPB$CX3 == FALSE ] = "CX2 only"
# dfPB$Class[ dfPB$CX2 == FALSE & dfPB$CX5 == TRUE & dfPB$CX3 == FALSE ] = "CX5 only"
# dfPB$Class[ dfPB$CX2 == FALSE & dfPB$CX5 == FALSE & dfPB$CX3 == TRUE ] = "CX3 only"
# dfPB$Class[ dfPB$CX2 == TRUE & dfPB$CX5 == TRUE & dfPB$CX3 == FALSE ] = "CX2+CX5"
# dfPB$Class[ dfPB$CX2 == FALSE & dfPB$CX5 == TRUE & dfPB$CX3 == TRUE ] = "CX5+CX3"
# dfPB$Class[ dfPB$CX2 == TRUE & dfPB$CX5 == FALSE & dfPB$CX3 == TRUE ] = "CX2+CX3"
# dfPB$Class[ dfPB$CX2 == TRUE & dfPB$CX5 == TRUE & dfPB$CX3 == TRUE ] = "All present"
#
# dfPath$SigClass2 = factor(dfPB$Class[ match(dfPath$Sample, rownames(dfPS)) ],
#                          levels = c("Background", "CX2 only", "CX5 only", "CX3 only",
#                                     "CX2+CX5", "CX2+CX3", "CX5+CX3", "All present"))
# dfPathNoNA = dfPath[ ! is.na(dfPath$SigClass2), ]
# p8 = ggplot(dfPathNoNA, aes(x = Pathway, y = SVs, fill = SigClass2)) +
#   # geom_jitter(width = 0.1, height = 0, alpha = 0.5) +
#   geom_boxplot(outlier.size = 0.25) +
#   coord_capped_cart(left = "both", bottom = "both") +
#   labs(y = "Signature activity", x = "CIN signature", fill = "Presence of\nIHR signature")
#
# # Not useful as we don't want all combinationes tested
# # dfPath %>% group_by(Pathway) %>% t_test(SVs ~ SigClass2, var.equal = FALSE, p.adjust.method = "BH")
#
#
# ## Specific cancer types
# ## Somehow hard to find all the metadata for those samples
# dfPath$Cancer = dtMeta$projectcode[ match(dfPath$Sample, dtMeta$samplename) ]
# ## sum(is.na(dfPath$Cancer)) => 921
#
# dfPath$Cancer[ is.na(dfPath$Cancer) ] = dtMeta2$project_code[ match(dfPath$Sample[ is.na(dfPath$Cancer) ],
#                                                               dtMeta2$tcga_donor_uuid) ]
# ## sum(is.na(dfPath$Cancer)) => 921
#
# dfPath$Cancer[ is.na(dfPath$Cancer) ] = dtMeta3$project_code[ match(dfPath$Sample[ is.na(dfPath$Cancer) ],
#                                                                     dtMeta3$samplename) ]
# ## sum(is.na(dfPath$Cancer)) => 921
#
# dfPath$Cancer[ is.na(dfPath$Cancer) ] = dtMeta4$Cancer[ match(dfPath$Sample[ is.na(dfPath$Cancer) ],
#                                                               dtMeta4$PCAWG) ]
# ## sum(is.na(dfPath$Cancer)) => 795
# ## 795/3 = 265 samples without cancer information
# ## So weird. Why do I have no cancer information for those 265 samples?
#
# ## Add order
# dfOrder = aggregate(SVs ~ Cancer, data = dfPath[ dfPath$Pathway == "NHEJ", ], mean)
# dfPath$Cancer = factor(dfPath$Cancer, levels = dfOrder$Cancer[ order(dfOrder$SVs, decreasing = TRUE) ])
#
# ## Plot 1: Show activity of repair pathway by tumour type
# p9 = ggplot(dfPath, aes(x = Cancer, y = SVs, colour = Pathway)) + geom_boxplot() +
#   scale_x_discrete(guide = guide_axis(angle = 90))
#
#
# ## Plot 2: Add dominant IHR signature
# dfPathNoNA = dfPath[ ! is.na(dfPath$SigClass),  ]
# p10 = ggplot(dfPathNoNA, aes(x = Pathway, y = SVs, fill = SigClass)) +
#   facet_wrap(. ~ Cancer) +
#   geom_boxplot(outlier.size = 0.25) +
#   coord_capped_cart(left = "both", bottom = "both") +
#   labs(y = "Signature activity", x = "DNA repair pathway", fill = "Presence of\nIHR signature")
#
#
# ## Select a few cancer types
# dfPathFilt = dfPathNoNA[ grepl("BRCA|GBM|OV-|UCEC", dfPathNoNA$Cancer), ]
# dfPathFilt$Cancer = sapply(str_split(dfPathFilt$Cancer, pattern = "-"), function(x) return(x[1]))
# p11 = ggplot(dfPathFilt, aes(x = Pathway, y = SVs, fill = SigClass)) +
#   facet_wrap(. ~ Cancer) +
#   geom_boxplot(outlier.size = 0.25) +
#   coord_capped_cart(left = "both", bottom = "both") +
#   labs(y = "Signature activity", x = "DNA repair pathway", fill = "Presence of\nIHR signature")
#
#
#
# ## By repair pathway
# TOPX = 0.9
# # quantile(dfPath$SVs[ dfPath$Pathway == "NHEJ" ], probs = TOPX)
# dfTMEJ = dfPath[ dfPath$Pathway == "TMEJ", ]
# threshTMEJ = quantile(dfTMEJ$SVs, probs = TOPX)
# dfTMEJ$ClassPathway = cut(dfTMEJ$SVs,breaks = c(0, threshTMEJ, 1), labels = c("Low", "High"),
#                           include.lowest = TRUE)
# p12 = ggplot(dfTMEJ, aes(x = ClassPathway, fill = SigClass)) + geom_bar(position = "fill") +
#   labs(x = "Top 10% quantile of TMEJ samples", y = "Proportion", fill = "Dominant IHR\nsignature")
#
#
#
#
#
# #### Correlation with signatures
# dfMPath = data.frame(mPath)
# allSigs = levels(dtPCAWG$Var2)
# lCors= lapply(allSigs, function(thisSig) {
#
#   dtSig = dtPCAWG[ dtPCAWG$Var2 == thisSig, ]
#   dfMPath$Activity = dtSig$value[ match(rownames(dfMPath), dtSig$Var1) ]
#
#   ## Do correlations manually - boooh!
#   corNHEJ = cor.test(dfMPath$NHEJ, dfMPath$Activity, method = 'spearman', exact = FALSE)
#   corTMEJ = cor.test(dfMPath$TMEJ, dfMPath$Activity, method = 'spearman', exact = FALSE)
#   corSSA = cor.test(dfMPath$SSA, dfMPath$Activity, method = 'spearman', exact = FALSE)
#
#   ## Create output
#   dtOut = data.table("Sig" = thisSig, "Pathway" = c("NHEJ", "TMEJ", "SSA"),
#                      "Rho" = c(signif(corNHEJ$estimate, 4), signif(corTMEJ$estimate, 4),
#                                signif(corSSA$estimate, 4)),
#                      "pVal" = c(signif(corNHEJ$p.value, 4), signif(corTMEJ$p.value, 4),
#                                 signif(corSSA$p.value, 4)))
#   return(dtOut)
#
# })
#
# dtCors = rbindlist(lCors)
# dtCors$pAdj = p.adjust(dtCors$pVal, method = "BH")
# dtCors$Plot = factor(dtCors$pAdj < 0.05, levels = c("TRUE", "FALSE"))
#
# ## Plot
# dtCors$Sig = factor(dtCors$Sig, levels = paste0("CX", 1:17))
# dtCors$Pathway = factor(dtCors$Pathway, levels = c("NHEJ", "TMEJ", "SSA"))
# pC = ggplot(dtCors, aes(x = Sig, y = Rho, group = Pathway, colour = Pathway)) + geom_line() +
#   geom_point(size = 2, aes(shape = Plot)) + coord_capped_cart(bottom = "both", left = "both")
#
# ## Just IHR sigs
# dtIHR = dtCors[ dtCors$Sig %in% c("CX2", "CX3", "CX5"), ]
# dtIHR$Sig = factor(as.character(dtIHR$Sig), levels = c("CX2", "CX5", "CX3"))
# pD = ggplot(dtIHR, aes(x = Sig, y = Rho, group = Pathway, colour = Pathway)) + geom_line() +
#   geom_point(size = 2) + coord_capped_cart(bottom = "both", left = "both")
#
#
# # (pA+pB)/pC
#
#
# #### Classify samples by SSA and see how sigs differ
# ## SSA happens at such low frequencies that the mere presence of it should be a telltale sign
# dtSSA = data.table("Sample" = rownames(mPath),
#                    "SSA" = factor(mPath[,"SSA"] > 0, levels = c("TRUE", "FALSE"),
#                                   labels = c("Active", "Inactive")))
#
# dtPCAWGScaled = melt(scale(readRDS(PCAWGACT)))
# dtPCAWGScaled$SSA = dtSSA$SSA[ match(dtPCAWG$Var1, dtSSA$Sample) ]
# dtPCAWGScaled = dtPCAWGScaled[ ! is.na(dtPCAWGScaled$SSA), ]
#
# pE = ggplot(dtPCAWGScaled, aes(x = Var2, y = value, colour = SSA)) +
#   geom_boxplot() + coord_cartesian(ylim = c(-2.2, 4))
#
# ## IHR sigs only
# dtIHRScaled = dtPCAWGScaled[ dtPCAWGScaled$Var2 %in% c("CX2", "CX3", "CX5"), ]
# dtIHRScaled$Var2 = factor(as.character(dtIHRScaled$Var2), levels = c("CX2", "CX5", "CX3"))
# pF = ggplot(dtIHRScaled, aes(x = SSA, y = value, colour = Var2)) +
#   geom_boxplot() + coord_cartesian(ylim = c(-2.2, 4))
#
# ## Check p-values
# dtPVals = dtIHRScaled[dtIHRScaled$SSA=="Active",]
# pairwise.t.test(dtPVals$value, dtPVals$Var2, pool.sd = FALSE, p.adjust.method = "BH")
#
#
#
# #### Classify samples by their dominant signature and then compare sigs
# dtDom = data.table("Sample" = rownames(mPath),
#                    "DominantRepair" = colnames(mPath)[ apply(mPath, 1, which.max) ] )
# # dtDom$DominantRepair[ is.na(dtDom$DominantRepair) ] = "NHEJ-TMEJ"
#
# dtPCAWGScaled = melt(scale(readRDS(PCAWGACT)))
# dtPCAWGScaled$DominantRepair = dtDom$DominantRepair[ match(dtPCAWG$Var1, dtDom$Sample) ]
# dtPCAWGScaled = dtPCAWGScaled[ ! is.na(dtPCAWGScaled$DominantRepair), ]
#
# pG = ggplot(dtPCAWGScaled, aes(x = Var2, y = value, colour = DominantRepair)) +
#   geom_boxplot() + coord_cartesian(ylim = c(-2.2, 4))
#
# ## Just IHR sigs
# dtIHRScaled = dtPCAWGScaled[ dtPCAWGScaled$Var2 %in% c("CX2", "CX3", "CX5"), ]
# dtIHRScaled$Var2 = factor(as.character(dtIHRScaled$Var2), levels = c("CX2", "CX5", "CX3"))
# pH = ggplot(dtIHRScaled, aes(x = Var2, y = value, colour = DominantRepair)) +
#   geom_boxplot() + coord_cartesian(ylim = c(-2.2, 4))
#
#
# stat.test <- dtIHRScaled %>%
#   group_by(DominantRepair) %>%
#   t_test(value ~ Var2, var.equal = FALSE, p.adjust.method = "BH")

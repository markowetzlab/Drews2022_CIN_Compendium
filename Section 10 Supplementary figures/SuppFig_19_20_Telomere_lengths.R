#### Test whole chromosomes and CIN signatures
rm(list=ls(all=TRUE))

library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(lemon)
library(patchwork)
# Robust linear regression
library(MASS)
# F (Wald) test
library(sfsmisc)
library(this.path)

theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))

BASE=dirname(this.path())
OUT=file.path(BASE, "output")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

ACT=file.path(BASE, "input/Signature_Compendium_v5_Cosine-0.74_Activities_THRESH95_NAMESAPRIL21.rds")
TL=file.path(BASE, "input/Barthel2017_SuppTab1_TelomereLength_TCGA.csv")


## Load data
mAct = readRDS(ACT)
dtAct = data.table(melt(scale(mAct, center = FALSE)))
colnames(dtAct) = c("Sample", "Signature", "Activity")

tlFull = fread(TL)


## Functions
# Source: https://stackoverflow.com/questions/52297978/decrease-overal-legend-size-elements-and-text
addSmallLegend = function(myPlot, pointSize = 0.5, textSize = 3, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize),
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}






# Remove sWGS and WXS since I'm not sure about how good TelSeq estimates TL lengths from off-target reads.
# Cannot access supp material of Ding et al., NAR (2014). Sent email to Richard Durbin on 29/10/19.
# Supp Figure 2, Panel A and C in Barthel et al., Nat Gen (2019) shows quite strong differences between sequencing technologies.
# So for now focussing on WGS-derived TL values only.
tl = tlFull[ tlFull$LibraryType == "WGS", ]


## Seven comparisons:
### - telomere length,
### - Telomerase Signature Score,
### - TERT and TERC amp, 
### - TERT expression, 
### - TERC expression, 
### - TERRA expression, 
### - ATRX/DAXX amplification

## Transfer telomere length
dtAct$tTL = tl$tTL[ match(substr(dtAct$Sample,1,12), substr(tl$PatientID,1,12)) ]
dtAct$nTL = tl$nTL[ match(substr(dtAct$Sample,1,12), substr(tl$PatientID,1,12)) ]
dtAct$tlLength = dtAct$tTL - dtAct$nTL

## Telomerase Signature Score - TSS
dtAct$TSS = tl$TelomeraseSignatureScore[ match(substr(dtAct$Sample,1,12), substr(tl$PatientID, 1,12)) ]

## TERT and TERC amplification
dtAct$TERTamp = tl$TERTamp[ match(substr(dtAct$Sample,1,12), substr(tl$PatientID, 1,12)) ]
dtAct$TERCamp = tl$TERCamp[ match(substr(dtAct$Sample,1,12), substr(tl$PatientID, 1,12)) ]
dtAct$TERTamp[ dtAct$TERTamp == "" ] = NA
dtAct$TERCamp[ dtAct$TERCamp == "" ] = NA

dtAct$TERTTERC = NA
dtAct$TERTTERC[ is.na(dtAct$TERTTERC) & dtAct$TERTamp == "Amplification" & dtAct$TERCamp == "Amplification" ] = "Both amplified"
dtAct$TERTTERC[ is.na(dtAct$TERTTERC) & dtAct$TERTamp == "Amplification" ] = "TERT amplified"
dtAct$TERTTERC[ is.na(dtAct$TERTTERC) & dtAct$TERCamp == "Amplification" ] = "TERC amplified"
dtAct$TERTTERC[ is.na(dtAct$TERTTERC) & dtAct$TERTamp == "Neutral/Deletion" & dtAct$TERCamp == "Neutral/Deletion" ] = "Both WT"

dtAct$TERTTERC = factor(dtAct$TERTTERC, levels = c("Both amplified", "TERT amplified", "TERC amplified", "Both WT"))



## TERT and TERC expression
dtAct$TERTexpr = tl$TERTexpr[ match(substr(dtAct$Sample,1,12), substr(tl$PatientID, 1,12)) ]
dtAct$TERCexpr = tl$TERCexpr[ match(substr(dtAct$Sample,1,12), substr(tl$PatientID, 1,12)) ]


## TERRA expression
dtAct$TERRA = tl$TERRAexpr[ match(substr(dtAct$Sample,1,12), substr(tl$PatientID,1,12)) ]
dtAct$TERRA[ dtAct$TERRA == "" ] = NA
dtAct$TERRA = factor(dtAct$TERRA, levels = c("TERRA expr", "TERRA non-expr"), labels = c("Expressed", "Non-expressed"))


## ATRX/DAXX expression
dtAct$ATRXDAXXstatus = tl$ATRXDAXXstatus[ match(substr(dtAct$Sample,1,12), substr(tl$PatientID,1,12)) ]
dtAct$ATRXDAXXstatus[ dtAct$ATRXDAXXstatus == "" ] = NA
dtAct$ATRXDAXXstatus = factor(dtAct$ATRXDAXXstatus, levels = c("ATRX/DAXX alt", "ATRX/DAXX wt"), labels = c("Mutated", "WT"))


#### Do all the statistical tests now
allSigs = levels(dtAct$Signature)
lLength = lapply(allSigs, function(thisSig) {
  
  ## Test 1: telomere length
  thisFilt1 = dtAct[ dtAct$Signature == thisSig & ! is.na(dtAct$tlLength), ]
  rlmModel1 = rlm(tlLength ~ Activity, data = thisFilt1)
  fTestModel1 = f.robftest(rlmModel1, var = -1)

  outNums1 = signif(c(coefficients(rlmModel1)[1], coefficients(rlmModel1)[2], fTestModel1$statistic, fTestModel1$p.value), 4)
  
  
  ## Test 2: Telomerase Signature Score
  thisFilt2 = dtAct[ dtAct$Signature == thisSig & ! is.na(dtAct$TSS), ]
  corTest2 = cor.test(thisFilt2$Activity, thisFilt2$TSS, method = "spearman", exact = FALSE, use = "pairwise")
  
  outNums2 = signif(c(corTest2$estimate, corTest2$p.value), 4)
  
  
  ## Test 3: TERT / TERC amplification
  dtFilt3 = dtAct[ dtAct$Signature == thisSig & ! is.na(dtAct$TERTTERC), ]
  # Three tests needed: WT-TERC, WT-TERT, WT-BOTH
  test31 = t.test(dtFilt3$Activity[ dtFilt3$TERTTERC == "Both WT" ], dtFilt3$Activity[ dtFilt3$TERTTERC == "TERC amplified" ], var.equal = FALSE)
  change31 = test31$estimate[2] - test31$estimate[1]
  
  test32 = t.test(dtFilt3$Activity[ dtFilt3$TERTTERC == "Both WT" ], dtFilt3$Activity[ dtFilt3$TERTTERC == "TERT amplified" ], var.equal = FALSE)
  change32 = test32$estimate[2] - test32$estimate[1]
  
  test33 = t.test(dtFilt3$Activity[ dtFilt3$TERTTERC == "Both WT" ], dtFilt3$Activity[ dtFilt3$TERTTERC == "Both amplified" ], var.equal = FALSE)
  change33 = test33$estimate[2] - test33$estimate[1]
  
  outNums3 = signif(c(change31, test31$p.value, change32, test32$p.value, change33, test33$p.value), 4)
  
  
  ## Test 4: TERC expression
  thisFilt4 = dtAct[ dtAct$Signature == thisSig & ! is.na(dtAct$TERCexpr), ]
  corTest4 = cor.test(thisFilt4$Activity, thisFilt4$TERCexpr, method = "spearman", exact = FALSE, use = "pairwise")
  
  outNums4 = signif(c(corTest4$estimate, corTest4$p.value), 4)
  
  
  ## Test 5: TERT expression
  thisFilt5 = dtAct[ dtAct$Signature == thisSig & ! is.na(dtAct$TERTexpr), ]
  corTest5 = cor.test(thisFilt5$Activity, thisFilt5$TERTexpr, method = "spearman", exact = FALSE, use = "pairwise")
  
  outNums5 = signif(c(corTest5$estimate, corTest5$p.value), 4)
  
  
  ## Test 6: TERRA expression
  thisFilt6 = dtAct[ dtAct$Signature == thisSig & ! is.na(dtAct$TERRA), ]
  test6 = t.test(thisFilt6$Activity[ thisFilt6$TERRA == "Non-expressed" ], thisFilt6$Activity[ thisFilt6$TERRA == "Expressed" ], var.equal = FALSE)
  change6 = test6$estimate[2] - test6$estimate[1]
  
  outNums6 = signif(c(change6, test6$p.value), 4)
  
  
  ## Test 7: ATRX/DAXX amplification
  dtFilt7 = dtAct[ dtAct$Signature == thisSig & ! is.na(dtAct$ATRXDAXXstatus), ]
  test7 = t.test(dtFilt7$Activity[ dtFilt7$ATRXDAXXstatus == "WT" ], dtFilt7$Activity[ dtFilt7$ATRXDAXXstatus == "Mutated" ], var.equal = FALSE)
  change7 = test7$estimate[2] - test7$estimate[1]
  
  outNums7 = signif(c(change7, test7$p.value), 4)
  
  out = c(thisSig, outNums1, outNums2, outNums3, outNums4, outNums5, outNums6, outNums7)
  return(out)
  
})


## Combine results
dtLength = data.table(do.call(rbind, lLength))
colnames(dtLength) = c("Signature", "Intercept", "Slope", "F", "pValRLM", "RhoTSS", "pValTSS", "DiffWTTERC", "pValWTTERC", "DiffWTTERT", "pValWTTERT", "DiffWTBoth", "pValWTBoth",
                       "RhoTERC", "pValTERC", "RhoTERT", "pValTERT", "DiffTERRA", "pValTERRA", "DiffATL", "pValATL")
dtLength$Intercept = as.numeric(dtLength$Intercept)
dtLength$Slope = as.numeric(dtLength$Slope)
dtLength$F = as.numeric(dtLength$F)

dtLength$RhoTSS = as.numeric(dtLength$RhoTSS)
dtLength$pValTSS = as.numeric(dtLength$pValTSS)

dtLength$DiffWTTERC = as.numeric(dtLength$DiffWTTERC)
dtLength$pValWTTERC = as.numeric(dtLength$pValWTTERC)
dtLength$DiffWTTERT = as.numeric(dtLength$DiffWTTERT)
dtLength$pValWTTERT = as.numeric(dtLength$pValWTTERT)
dtLength$DiffWTBoth = as.numeric(dtLength$DiffWTBoth)
dtLength$pValWTBoth = as.numeric(dtLength$pValWTBoth)

dtLength$RhoTERC = as.numeric(dtLength$RhoTERC)
dtLength$pValTERC = as.numeric(dtLength$pValTERC)
dtLength$RhoTERT = as.numeric(dtLength$RhoTERT)
dtLength$pValTERT = as.numeric(dtLength$pValTERT)

dtLength$DiffTERRA = as.numeric(dtLength$DiffTERRA)
dtLength$pValTERRA = as.numeric(dtLength$pValTERRA)

dtLength$DiffATL = as.numeric(dtLength$DiffATL)
dtLength$pValATL = as.numeric(dtLength$pValATL)


## Correct p-values
dtLength$pAdjRLM = p.adjust(dtLength$pValRLM, method = "BH")
dtLength$pAdjTSS = p.adjust(dtLength$pValTSS, method = "BH")
dtLength$pAdjWTTERC = p.adjust(dtLength$pValWTTERC, method = "BH")
dtLength$pAdjWTTERT = p.adjust(dtLength$pValWTTERT, method = "BH")
dtLength$pAdjWTBoth = p.adjust(dtLength$pValWTBoth, method = "BH")
dtLength$pAdjTERC = p.adjust(dtLength$pValTERC, method = "BH")
dtLength$pAdjTERT = p.adjust(dtLength$pValTERT, method = "BH")
dtLength$pAdjTERRA = p.adjust(dtLength$pValTERRA, method = "BH")
dtLength$pAdjATL = p.adjust(dtLength$pValATL, method = "BH")

write.table(dtLength, file.path(OUT, "Covariates_telomeres.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


#### Plotting!
## Panel A: Robust linear model of telomere length
pA = ggplot(dtLength, aes(x = Slope, y = -log10(pAdjRLM), colour = ifelse(dtLength$pAdjRLM < 0.05, "Significant", "NS"))) + 
  geom_point() + geom_text_repel(aes(label = Signature), force = 100) + 
  coord_capped_cart(left = "both", bottom = "both") +
  labs(x = "Slope of robust linear model", y = "-log(q-value)", title = "Telomere length [kB]", 
       colour = "q-Value") + 
  geom_vline(xintercept = 0, colour = "grey30", linetype = "dashed") +
  theme(legend.position = c(0.4, 0.8)) +
  scale_colour_manual(values = c("grey80", "black"))

pB = ggplot(dtLength, aes(x = RhoTSS, y = -log10(pAdjTSS), colour = ifelse(dtLength$pAdjTSS < 0.05, "Significant", "NS"))) + 
  geom_point() + geom_text_repel(aes(label = Signature), force = 100) + 
  coord_capped_cart(left = "both", bottom = "both") +
  labs(x = "Spearman's Rho", y = "-log(q-value)", title = "Telomerase Signature Score", 
       colour = "q-Value") + 
  geom_vline(xintercept = 0, colour = "grey30", linetype = "dashed") +
  theme(legend.position = c(0.4, 0.8)) +
  scale_colour_manual(values = c("grey80", "black"))


## Prepare panel C
dtLength$Significant = factor(dtLength$pAdjWTBoth < 0.05 | dtLength$pAdjWTTERC < 0.05 | dtLength$pAdjWTTERT < 0.05, levels = c(TRUE, FALSE), labels = c("Significant", "NS"))
dtAct$Label = dtLength$Significant[ match(dtAct$Signature, dtLength$Signature) ]
dtAct$Label = paste0(dtAct$TERTTERC, "_", dtAct$Label)
dtAct$Label[ is.na(dtAct$TERTTERC) ] = NA

dtAct$Label = factor(dtAct$Label, levels = c("Both amplified_NS", "TERC amplified_NS", "TERT amplified_NS", "Both WT_NS",
                                             "Both amplified_Significant", "TERC amplified_Significant", "TERT amplified_Significant", "Both WT_Significant"))

pC = ggplot(dtAct[ ! is.na(dtAct$TERTTERC),], aes(x = Signature, y = Activity, fill = Label)) + geom_boxplot(outlier.size = 0.1) +
  coord_capped_cart(left = "both", bottom = "both", ylim = c(-0.1, 3)) + labs(y = "Scaled activity") + theme(legend.position = c(0.8, 0.8)) +
  scale_fill_manual(values = c("Both amplified_NS" = "white", "TERC amplified_NS" = "white", "TERT amplified_NS" = "white", "Both WT_NS" = "white",
                               "Both amplified_Significant" = "#fb8072", "TERC amplified_Significant" = "#ffffb3", "TERT amplified_Significant" = "#bebada", "Both WT_Significant" = "#8dd3c7"))
pC = addSmallLegend(myPlot = pC, textSize = 6, spaceLegend = 0.5)


pD = ggplot(dtLength, aes(x = RhoTERC, y = -log10(pAdjTERC), colour = ifelse(dtLength$pAdjTERC < 0.05, "Significant", "NS"))) + 
  geom_point() + geom_text_repel(aes(label = Signature), force = 100) + 
  coord_capped_cart(left = "both", bottom = "both") +
  labs(x = "Spearman's Rho", y = "-log(q-value)", title = "TERC expression", 
       colour = "q-Value") + 
  geom_vline(xintercept = 0, colour = "grey30", linetype = "dashed") +
  theme(legend.position = c(0.4, 0.8)) +
  scale_colour_manual(values = c("grey80", "black"))

pE = ggplot(dtLength, aes(x = RhoTERT, y = -log10(pAdjTERT), colour = ifelse(dtLength$pAdjTERT < 0.05, "Significant", "NS"))) + 
  geom_point() + geom_text_repel(aes(label = Signature), force = 100) + 
  coord_capped_cart(left = "both", bottom = "both") +
  labs(x = "Spearman's Rho", y = "-log(q-value)", title = "TERT expression", 
       colour = "q-Value") + 
  geom_vline(xintercept = 0, colour = "grey30", linetype = "dashed") +
  theme(legend.position = c(0.4, 0.8)) +
  scale_colour_manual(values = c("grey80", "black"))



## Prepare panel F
dtLength$SignificantTERRA = factor(dtLength$pAdjTERRA < 0.05, levels = c(TRUE, FALSE), labels = c("Significant", "NS"))
dtAct$LabelTERRA = dtLength$SignificantTERRA[ match(dtAct$Signature, dtLength$Signature) ]
dtAct$LabelTERRA = paste0(dtAct$TERRA, "_", dtAct$LabelTERRA)
dtAct$LabelTERRA[ is.na(dtAct$TERRA) ] = NA

dtAct$LabelTERRA = factor(dtAct$LabelTERRA, levels = c("Expressed_NS", "Non-expressed_NS" , "Expressed_Significant", "Non-expressed_Significant"))

pF = ggplot(dtAct[ ! is.na(dtAct$LabelTERRA), ], aes(x = Signature, y = Activity, fill = LabelTERRA)) + geom_boxplot(outlier.size = 0.1) +
  coord_capped_cart(left = "both", bottom = "both", ylim = c(-0.1, 3)) + labs(y = "Scaled activity") + theme(legend.position = c(0.8, 0.8)) +
  scale_fill_manual(values = c("Expressed_NS" = "white", "Non-expressed_NS" = "white" , "Expressed_Significant" = "red", "Non-expressed_Significant" = "blue"))
pF = addSmallLegend(myPlot = pF, textSize = 6, spaceLegend = 0.5)

## Prepare panel G
dtLength$SignificantATL = factor(dtLength$pAdjATL < 0.05, levels = c(TRUE, FALSE), labels = c("Significant", "NS"))
dtAct$LabelATL = dtLength$SignificantATL[ match(dtAct$Signature, dtLength$Signature) ]
dtAct$LabelATL = paste0(dtAct$ATRXDAXXstatus, "_", dtAct$LabelATL)
dtAct$LabelATL[ is.na(dtAct$ATRXDAXXstatus) ] = NA

dtAct$LabelATL = factor(dtAct$LabelATL, levels = c("Mutated_NS", "WT_NS" , "Mutated_Significant", "Non-WT_Significant"))

pG = ggplot(dtAct[ ! is.na(dtAct$LabelATL), ], aes(x = Signature, y = Activity, fill = LabelATL)) + geom_boxplot(outlier.size = 0.1) +
  coord_capped_cart(left = "both", bottom = "both", ylim = c(-0.1, 3)) + labs(y = "Scaled activity") + theme(legend.position = c(0.8, 0.8)) +
  scale_fill_manual(values = c("Mutated_NS" = "white", "WT_NS" = "white" , "Mutated_Significant" = "yellow", "Non-WT_Significant" = "black"), drop = FALSE)
pG = addSmallLegend(myPlot = pG, textSize = 6, spaceLegend = 0.5)


## Save output
ggsave(file.path(OUT, "SuppFig_19_Telomere_lengths_Panel_A.svg"), pA, width = 60, height = 45, units = "mm")
ggsave(file.path(OUT, "SuppFig_19_Telomere_lengths_Panel_B.svg"), pB, width = 60, height = 45, units = "mm")
ggsave(file.path(OUT, "SuppFig_19_Telomere_lengths_Panel_C.svg"), pC, width = 120, height = 60, units = "mm")
ggsave(file.path(OUT, "SuppFig_19_Telomere_lengths_Panel_D.svg"), pD, width = 60, height = 45, units = "mm")
ggsave(file.path(OUT, "SuppFig_19_Telomere_lengths_Panel_E.svg"), pE, width = 60, height = 45, units = "mm")
ggsave(file.path(OUT, "SuppFig_20_Telomere_lengths_Panel_A.svg"), pF, width = 120, height = 60, units = "mm")
ggsave(file.path(OUT, "SuppFig_20_Telomere_lengths_Panel_B.svg"), pG, width = 120, height = 60, units = "mm")

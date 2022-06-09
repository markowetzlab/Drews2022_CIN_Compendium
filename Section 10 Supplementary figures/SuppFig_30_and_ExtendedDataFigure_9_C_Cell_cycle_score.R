## Correlating cell cycle scores with CX signature activities

rm(list=ls(all=TRUE))

library(this.path)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(lemon)
library(rstatix)

theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))

BASE=dirname(this.path())
OUT=file.path(BASE, "output")
OUTTABLE=file.path(BASE, "output")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTTABLE, showWarnings = FALSE, recursive = TRUE)

MEDIANSCALE=FALSE
EXP=file.path(BASE, "input/Signature_Compendium_v5_Cosine-0.74_Activities_THRESH95_NAMESAPRIL21.rds")
CCS=file.path(BASE, "input/Lundberg_2020_Cell_cycle_score.rds")


## Load data and merge cell cycle scores with signature activities
## For better visualisation, scaling is done with robust scaling using median and sd
if(MEDIANSCALE) {
  mExp = readRDS(EXP)
  scExp = apply(mExp, 2, function(x) (x-median(x))/sd(x))
  dtExp = data.table(melt(scExp))
} else {
  dtExp =  data.table(melt(scale(readRDS(EXP))))
}

dtCCS = data.table(readRDS(CCS))
dtCCS$Sample = substr(dtCCS$Sample, 1, 12)

dtExp$CCS = dtCCS$CCS[ match(substr(dtExp$Var1,1,12), dtCCS$Sample) ]
## According to Lundberg's script on github
dtExp$CCS = factor(as.numeric(dtExp$CCS), levels = c(1,2,3), labels = c("Low", "Medium", "High"))

## Not really needed
dtExp$CCS_ct = dtCCS$CCS_ct[ match(substr(dtExp$Var1,1,12), dtCCS$Sample) ]
dtExp$CCS_ct = (dtExp$CCS_ct-min(dtExp$CCS_ct, na.rm = TRUE))/
  (max(dtExp$CCS_ct, na.rm = TRUE)-min(dtExp$CCS_ct, na.rm = TRUE))

## Identify signatures where high is significantly different to med and low
allSigs = levels(dtExp$Var2)
lTTest = lapply(allSigs, function(thisSig) {
  
  print(thisSig)
  dtSig = dtExp[ dtExp$Var2 == thisSig, ]
  
  ## Only interested when high is significant to medium and low
  tHighVMed = t.test(dtSig$value[ dtSig$CCS == "High" ], 
                      dtSig$value[ dtSig$CCS == "Medium" ], var.equal = FALSE)
  tHighVLow = t.test(dtSig$value[ dtSig$CCS == "High" ], 
                     dtSig$value[ dtSig$CCS == "Low" ], var.equal = FALSE)
  out = c(thisSig, tHighVMed$statistic, tHighVMed$p.value, 
      tHighVLow$statistic, tHighVLow$p.value)
  return(out)
  
})

dtCor = data.table(do.call(rbind, lTTest))
colnames(dtCor) = c("Sig", "HighVMed_t", "HighVMed_pVal", "HighVLow_t","HighVLow_pVal")
dtCor$HighVMed_t = signif(as.numeric(dtCor$HighVMed_t), 4)
dtCor$HighVMed_pVal = signif(as.numeric(dtCor$HighVMed_pVal), 4)
dtCor$HighVLow_t = signif(as.numeric(dtCor$HighVLow_t), 4)
dtCor$HighVLow_pVal = signif(as.numeric(dtCor$HighVLow_pVal), 4)

## Correct p-values and identify signatures where both tests are significant
dtCor$HighVMed_pAdj = p.adjust(dtCor$HighVMed_pVal, method = "BH")
dtCor$HighVLow_pAdj = p.adjust(dtCor$HighVLow_pVal, method = "BH")
dtCor$Sign = dtCor$HighVMed_pAdj < 0.05 & dtCor$HighVLow_pAdj < 0.05

## Plot
dtExp$SigSign = dtCor$Sign[ match(dtExp$Var2, dtCor$Sig) ]

## Manual curation of results (plot and then curate)
## Sigs where directionality changes going from low to medium to high: CS14
## Sigs with negative directionality are also ignored: CS1, CS11
dtExp$SigSign[ dtExp$Var2 %in% c("CX1", "CX6", "CX14") ] = FALSE
dtExp$Col = paste(dtExp$CCS, dtExp$SigSign, sep = "_")

## Prepare plot
dtPlot = dtExp[ ! is.na(dtExp$CCS), ]
dtPlot$Col = factor(dtPlot$Col, levels = c("Low_FALSE", "Medium_FALSE", "High_FALSE",
                                           "Low_TRUE", "Medium_TRUE", "High_TRUE"))

pOut = ggplot(dtPlot, aes(x = Var2, y = value, fill = Col)) + 
  geom_hline(yintercept = 0, colour = "grey20", linetype = "dashed") +
  geom_boxplot(outlier.size = 0.15) +
  labs(x = "Signature", y = "Scaled activity", fill = "Cell cycle\nscore") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  # coord_cartesian(ylim = c(-3, 3)) +
  coord_capped_cart(left = "both", bottom = "both", ylim = c(-2,2)) +
  scale_fill_manual(values = c("Low_FALSE" = "white", "Medium_FALSE" = "white",  "High_FALSE" = "white", 
                               "Low_TRUE" = "black", "Medium_TRUE" = "grey", "High_TRUE" = "gold"))

cairo_pdf(file.path(OUT, "SuppFig_30_Cell_cycle_score.pdf"), width = 180/25.4, height =  120/25.4)
print(pOut); dev.off()

ggsave(file.path(OUT, "SuppFig_30_Cell_cycle_score.svg"), pOut, width = 180, height = 120, units = "mm")

write.table(dtPlot, file.path(OUTTABLE, "Supp_Table_Cell_cycle_score.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
saveRDS(dtPlot, file.path(OUTTABLE, "Supp_Table_Cell_cycle_score.rds"))
write.table(dtCor, file.path(OUTTABLE, "Supp_Table_Cell_cycle_score_pVals.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


### Just three IHR signatures for Extended Data Figure 4C
dtIHR = dtPlot[ dtPlot$Var2 %in% c("CX2", "CX5", "CX3"), ]
dtIHR$Var2 = factor(dtIHR$Var2, levels = c("CX2", "CX5", "CX3"))

pOut2 = ggplot(dtIHR, aes(x = Var2, y = value, fill = CCS)) + 
  geom_hline(yintercept = 0, colour = "grey20", linetype = "dashed") +
  geom_boxplot(outlier.size = 0.15) + theme(legend.position = c(0.7, 0.7)) +
  labs(x = "Signature", y = "Scaled activity", fill = "Cell cycle\nscore") +
  coord_capped_cart(left = "both", bottom = "both", ylim = c(-1.2,2)) +
  scale_fill_manual(values = c("black", "grey", "gold"))

ggsave(file.path(OUT, "ExtendedDataFigure_4_C_Cell_cycle_score.svg"), pOut2, width = 60, height = 45, units = "mm")

## T-test
statTest <- dtIHR %>%
  group_by(Var2) %>%
  t_test(value ~ CCS, var.equal = FALSE, p.adjust.method = "BH")


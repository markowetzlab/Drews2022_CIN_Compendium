#### Correlation of SV signatures with CIN signatures
rm(list=ls(all=TRUE))

library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(lemon)
library(patchwork)
library(scales)
library(viridisLite)
library(this.path)

theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6))

## 560 Breast cancers SV exposures
BASE=dirname(this.path())
OUT=file.path(BASE, "output")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

SVSIGS = file.path(BASE, "input/Nik-Zainal2018_SMTable_21.txt")
BC560 = file.path(BASE, "input/ICGC_560BCs_signature_activities_THRESH095_NAMESAPRIL21.rds")

NORMALISESV = TRUE
THRESHPLOT = 0.2

mSV = fread(SVSIGS)
if(NORMALISESV) {
  mSV[,2:ncol(mSV)] = mSV[,2:ncol(mSV)] / rowSums(mSV[,2:ncol(mSV)])
}
dtSV = melt(mSV)
dtBC = melt(readRDS(BC560))

## Fit column names
dtBC$Var1 = as.character(dtBC$Var1)
dtBC$Var1 = substr(dtBC$Var1, 1, nchar(dtBC$Var1)-1)

## Correlate
allCINSigs = levels(dtBC$Var2)
allSVSigs = levels(dtSV$variable)
lCor = lapply(allCINSigs, function(thisSig) {
  
  dtSig = dtBC[ dtBC$Var2 == thisSig, ]
  dtSV$Activity = dtSig$value[ match(dtSV$Sample, dtSig$Var1) ]
  
  ## Correlate SV sigs with each CIN sigs
  lSV = lapply(allSVSigs, function(thisSVSig) {
    
    dtSVSig = dtSV[ dtSV$variable == thisSVSig, ]
    spear = cor.test(dtSVSig$value, dtSVSig$Activity, use = "complete.obs", method = "spearman")
    
    out = c(thisSig, thisSVSig, signif(spear$estimate, 4), signif(spear$p.value, 4))
    return(out)
  })
  
  dtCorSV = data.table(do.call(rbind, lSV))
  colnames(dtCorSV) = c("CINSignature", "SVSignature", "Spearman", "pVal")
  dtCorSV$Spearman = as.numeric(dtCorSV$Spearman)
  dtCorSV$pVal = as.numeric(dtCorSV$pVal)
  dtCorSV$pAdj = p.adjust(dtCorSV$pVal, method = "BH")
  
  return(dtCorSV)
  
})

dtCor = rbindlist(lCor)
dtCor$CINSignature = factor(dtCor$CINSignature, levels = allCINSigs)
dtCor$SVSignature = factor(dtCor$SVSignature, levels = allSVSigs)


## Filter for significance
dtCorSig = dtCor[ dtCor$pAdj < 0.05 & dtCor$Spearman > 0, ]
dtCorSig$Plot = -log(dtCorSig$pAdj)

if(THRESHPLOT) {
  dtCorSig = dtCorSig[ dtCorSig$Spearman > THRESHPLOT, ]
}

## Plot
pOut = ggplot(dtCorSig, aes(x = CINSignature, y = SVSignature, fill = Spearman)) + geom_tile(aes(width = 0.94, height = 0.94)) +
  scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE) + scale_fill_viridis_c() +
  labs(x = "CIN signatures", y = "SV signatures", fill = "Spearman\ncorrelation") +
  theme(aspect.ratio = 6/17, axis.line = element_line(size = 0.5), axis.ticks = element_line(size = 0.5),
        axis.ticks.length = unit(.1, "cm")) +
  coord_capped_cart(bottom = "both", left = "both")

## Save plots
# cairo_pdf(file.path(OUT, "SuppFig_23_SV_Signatures.pdf"), width = 160/25.4, height = 50/25.4)
# print(pOut); dev.off()

ggsave(file.path(OUT, "SuppFig_23_SV_Signatures.svg"), pOut, width = 160/25.4, height = 50/25.4)

write.table(dtCorSig, file.path(OUT, "Supp_Table_SV_Signatures_plot.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
#saveRDS(dtCorSig, file.path(OUTTABLE, "Supp_Table_SV_Signatures_plot.rds"))
write.table(dtCor, file.path(OUT, "Supp_Table_SV_Signatures.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




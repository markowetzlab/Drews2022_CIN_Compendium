## Correlate ploidy values with signature activities

rm(list=ls(all=TRUE))

library(this.path)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(lemon)
library(RColorBrewer)

theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))

## Basics
## Windows
BASE=dirname(this.path())
OUTPLOTS=file.path(BASE, "output")
OUTRESULTS=file.path(BASE, "output")
dir.create(OUTPLOTS, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTRESULTS, showWarnings = FALSE, recursive = TRUE)

## Signature files
EXP=file.path(BASE, "input/Signature_Compendium_v5_Cosine-0.74_Activities_THRESH95_NAMESAPRIL21.rds")
META=file.path(BASE, "input/Metadata_TCGA_ASCAT_penalty70.rds")
WGD=file.path(BASE, "input/Haase2019_TCGA.giScores.wgd.txt")

## Mode of test: Simple => Correlation, Complex => linear regression correcting for CX4 (THE WGD signature)
MODE="COMPLEX"


## Load data
exp = readRDS(EXP)
meta = readRDS(META)
wgd = fread(WGD)

## Prepare data
dtScore = data.table(melt(scale(exp)))
colnames(dtScore) = c("Name", "Sig", "Activity")
dtScore$Ploidy = meta$ploidy[ match(dtScore$Name, meta$name) ]
dtScore$WGD = wgd$wgd[ match(substr(dtScore$Name,1,12), wgd$patient) ]

## Add cX4 so that we can correct for it during the complex testing regime
sExp = scale(exp)
dtScore$CX4 = sExp[match(dtScore$Name, rownames(sExp)),"CX4"]


## Correlate activity with ploidy
allSigs = levels(dtScore$Sig)
lTest = lapply(allSigs, function(thisSig) {
  
  dtSig = dtScore[ dtScore$Sig == thisSig, ]
  
  if(MODE == "SIMPLE") {
    ## Correlation
    corTest = cor.test(dtSig$Activity, dtSig$Ploidy, method = "spearman")
    out = c(thisSig, signif(corTest$estimate, 4), signif(corTest$p.value, 4))
  } else {
    ## Linear regression correcting for CX4
    lmWGD = lm(Activity ~ Ploidy + CX4 + WGD*CX4, dtSig)
    out = c(thisSig, signif(summary(lmWGD)$coefficients[, "t value"]["Ploidy"], 4),
            signif(summary(lmWGD)$coefficients[, "Pr(>|t|)"]["Ploidy"], 4))
  }
  
  return(out)
  
})

dtTest = data.table(do.call(rbind, lTest))
colnames(dtTest) = c("Sig", "Statistic", "pVal")
dtTest$Statistic = as.numeric(dtTest$Statistic)
dtTest$pVal = as.numeric(dtTest$pVal)
dtTest$pAdj = p.adjust(dtTest$pVal)
dtTest$logPAdj = -log(dtTest$pAdj)
dtTest$Significant = dtTest$pAdj < 0.05

## Plot figure
## Remove CX4 if mode is not simple because we correct for CX4 this result is obviously crap
if(MODE != "SIMPLE") {
  dtPlot = dtTest[ dtTest$Sig != "CX4", ]
  p1 = ggplot(dtPlot, aes(x = Statistic, y = logPAdj, colour = Significant)) + 
    geom_point() + geom_text_repel(aes(label = Sig), force = 50, max.overlaps = 15) + 
    coord_capped_cart(left = "both", bottom = "both") +
    labs(x = "t Test Statistic", y = "-log(q-value)", title = "Sample ploidy") + 
    geom_vline(xintercept = 0, colour = "grey30", linetype = "dashed") +
    theme(legend.position = c(0.3, 0.8)) +
    scale_colour_manual(values = c("grey80", "black"))
} else {
  p1 = ggplot(dtTest, aes(x = Statistic, y = logPAdj, colour = Significant)) + 
    geom_point() + geom_text_repel(aes(label = Sig), force = 50, max.overlaps = 15) + 
    coord_capped_cart(left = "both", bottom = "both") +
    labs(x = "Spearman's Rho", y = "-log(q-value)", title = "Sample ploidy") + 
    geom_vline(xintercept = 0, colour = "grey30", linetype = "dashed") +
    theme(legend.position = c(0.3, 0.8)) +
    scale_colour_manual(values = c("grey80", "black"))
}


cairo_pdf(file.path(OUTPLOTS, "SuppFig_52_Ploidy_v_activity_ComplexLM.pdf"), width = 90/25.4, height = 90/25.4)
print(p1); dev.off()
ggsave(file.path(OUTPLOTS, "SuppFig_52_Ploidy_v_activity_ComplexLM.svg"), p1, width = 90/25.4, height = 90/25.4)
ggsave(file.path(OUTPLOTS, "SuppFig_52_Ploidy_v_activity_ComplexLM.png"), p1, width = 90/25.4, height = 90/25.4)

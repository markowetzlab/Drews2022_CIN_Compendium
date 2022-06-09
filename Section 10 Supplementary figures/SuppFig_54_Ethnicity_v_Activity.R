## Checking ethnicity distribution in the TCGA

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

BASE=dirname(this.path())
OUTFIGURES=file.path(BASE, "output")

EXP=file.path(BASE, "input/Signature_Compendium_v5_Cosine-0.74_Activities_THRESH95_NAMESAPRIL21.rds")
SURV=file.path(BASE, "input/Survival_and_BRCA_Status_fullTCGA.rds")

## Load and prepare activities
exp = readRDS(EXP)
surv = readRDS(SURV)

dtExp = data.table(melt(scale(exp)))
dtExp$Race = factor(surv$Race[ match(dtExp$Var1, surv$Name) ])
dtExp$TS = factor(surv$TS[ match(dtExp$Var1, surv$Name) ])
dtExp$Cancer = factor(surv$Cancer[ match(dtExp$Var1, surv$Name) ])
dtExp$AgeY = surv$AgeY[ match(dtExp$Var1, surv$Name) ]


## Linear model showing that race does have an influence
lmRace = lm(value ~ Var2 * Race + AgeY + TS + Cancer, dtExp)


### Simple version
## Do t-tests 

## Identify signatures where high is significantly different to med and low
allSigs = levels(dtExp$Var2)
lTTest = lapply(allSigs, function(thisSig) {
  
  print(thisSig)
  dtSig = dtExp[ dtExp$Var2 == thisSig, ]
  
  ## Only interested when high is significant to medium and low
  tWvB = t.test(dtSig$value[ dtSig$Race == "White" ], 
                     dtSig$value[ dtSig$Race == "Black" ], var.equal = FALSE)
  tWvA = t.test(dtSig$value[ dtSig$Race == "White" ], 
                dtSig$value[ dtSig$Race == "Asian" ], var.equal = FALSE)
  tWvO = t.test(dtSig$value[ dtSig$Race == "White" ], 
                dtSig$value[ dtSig$Race == "Other" ], var.equal = FALSE)
  out = c(thisSig, tWvB$statistic, tWvB$p.value, 
          tWvA$statistic, tWvA$p.value,
          tWvO$statistic, tWvO$p.value)
  return(out)
  
})

dtCor = data.table(do.call(rbind, lTTest))
colnames(dtCor) = c("Sig", "WhiteVBlack_t", "WhiteVBlack_pVal", 
                    "WhiteVAsian_t","WhiteVAsian_pVal",
                    "WhiteVOther_t","WhiteVOther_pVal")
dtCor$WhiteVBlack_t = signif(as.numeric(dtCor$WhiteVBlack_t), 4)
dtCor$WhiteVBlack_pVal = signif(as.numeric(dtCor$WhiteVBlack_pVal), 4)
dtCor$WhiteVAsian_t = signif(as.numeric(dtCor$WhiteVAsian_t), 4)
dtCor$WhiteVAsian_pVal = signif(as.numeric(dtCor$WhiteVAsian_pVal), 4)
dtCor$WhiteVOther_t = signif(as.numeric(dtCor$WhiteVOther_t), 4)
dtCor$WhiteVOther_pVal = signif(as.numeric(dtCor$WhiteVOther_pVal), 4)

## Correct p-values and identify signatures where both tests are significant
dtCor$WhiteVBlack_pAdj = p.adjust(dtCor$WhiteVBlack_pVal, method = "BH")
dtCor$WhiteVAsian_pAdj = p.adjust(dtCor$WhiteVAsian_pVal, method = "BH")
dtCor$WhiteVOther_pAdj = p.adjust(dtCor$WhiteVOther_pVal, method = "BH")
dtCor$Sign = dtCor$WhiteVBlack_pAdj < 0.001 | dtCor$WhiteVAsian_pAdj < 0.001 | dtCor$WhiteVOther_pAdj < 0.001


## Plot
dtExp$SigSign = dtCor$Sign[ match(dtExp$Var2, dtCor$Sig) ]


## Manual curation of results (plot and then curate)
## Sigs where directionality changes going from low to medium to high: CS14
## Sigs with negative directionality are also ignored: CS1, CS11
# dtExp$SigSign[ dtExp$Var2 %in% c("CX1", "CX6", "CX14") ] = FALSE
dtExp$Col = paste(dtExp$Race, dtExp$SigSign, sep = "_")

## Prepare plot
dtPlot = dtExp[ ! is.na(dtExp$Race), ]
dtPlot$Col = factor(dtPlot$Col, levels = c("Asian_FALSE", "Black_FALSE", "Other_FALSE", "White_FALSE",
                                           "Asian_TRUE", "Black_TRUE", "Other_TRUE", "White_TRUE"))

pOut = ggplot(dtPlot, aes(x = Var2, y = value, fill = Col)) + 
  geom_hline(yintercept = 0, colour = "grey20", linetype = "dashed") +
  geom_boxplot(outlier.size = 0.15) +
  labs(x = "Signature", y = "Scaled activity", fill = "Cell cycle\nscore") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  # coord_cartesian(ylim = c(-3, 3)) +
  coord_capped_cart(left = "both", bottom = "both", ylim = c(-2,2)) +
  scale_fill_manual(values = c("Asian_FALSE" = "grey", "Black_FALSE" = "grey", "Other_FALSE" = "grey", "White_FALSE" = "grey",
                               "Asian_TRUE" = "#fff7bc", "Black_TRUE" = "#ec7014", "Other_TRUE" = "#41ab5d", "White_TRUE" = "#f6eff7"))


cairo_pdf(file.path(OUTFIGURES, "SuppFig_54_Ethnicity_v_Activity.pdf"), width = 180/25.4, height =  120/25.4)
print(pOut); dev.off()

ggsave(file.path(OUTFIGURES, "SuppFig_54_Ethnicity_v_Activity.svg"), pOut, width = 180, height = 120, units = "mm")

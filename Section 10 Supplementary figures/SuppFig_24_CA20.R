## Correlating CA20 scores with CS12 activities

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
EXP=file.path(BASE, "input/Signature_Compendium_v5_Cosine-0.74_Activities_THRESH95_NAMESAPRIL21.rds")
CA20=file.path(BASE, "input/Almeida_2019_S1_CA20_Scores_TCGA.txt")
OUT=file.path(BASE, "output")


## Load data and merge CA20 with signature activities
dtExp = data.table(melt(scale(readRDS(EXP))))
dtCa = fread(CA20)
dtCa$Sample = gsub("\\.", "-", substr(dtCa$`Sample ID`, 1,12))
dtExp$CA20 = dtCa$CA20[ match(substr(dtExp$Var1,1,12), dtCa$Sample) ]


## Correlate data and plot results
allSigs = levels(dtExp$Var2)
lCor = lapply(allSigs, function(thisSig) {
  
  print(thisSig)
  dtSig = dtExp[ dtExp$Var2 == thisSig, ]
  out = c(thisSig, cor(dtSig$value, dtSig$CA20, use = "pairwise.complete.obs", method = "kendall"),
          cor.test(dtSig$value, dtSig$CA20, use = "pairwise.complete.obs", method = "kendall")$p.value,
          cor(dtSig$value, dtSig$CA20, use = "pairwise.complete.obs", method = "pearson"),
          cor.test(dtSig$value, dtSig$CA20, use = "pairwise.complete.obs", method = "pearson")$p.value,
          cor(dtSig$value, dtSig$CA20, use = "pairwise.complete.obs", method = "spearman"),
          cor.test(dtSig$value, dtSig$CA20, use = "pairwise.complete.obs", method = "spearman")$p.value)
  return(out)
  
})

dtCor = data.table(do.call(rbind, lCor))
colnames(dtCor) = c("Sig", "KendallsTau", "KendallPVal", "PearsonsR","PearsonPVal", 
                    "SpearmansRho", "SpearmanPVal")
dtCor$KendallsTau = signif(as.numeric(dtCor$KendallsTau), 4)
dtCor$KendallPVal = signif(as.numeric(dtCor$KendallPVal), 4)
dtCor$PearsonsR = signif(as.numeric(dtCor$PearsonsR), 4)
dtCor$PearsonPVal = signif(as.numeric(dtCor$PearsonPVal), 4)
dtCor$SpearmansRho = signif(as.numeric(dtCor$SpearmansRho), 4)
dtCor$SpearmanPVal = signif(as.numeric(dtCor$SpearmanPVal), 4)

dtMCor = melt(dtCor, id.vars = "Sig", measure.vars = c("KendallsTau", "PearsonsR", "SpearmansRho"))
dtMCorSig = melt(dtCor, id.vars = "Sig", measure.vars = c("KendallPVal", "PearsonPVal", "SpearmanPVal"))
dtMCor$PVal = dtMCorSig$value

dtMCor$PAdj = p.adjust(dtMCor$PVal, method = "BH")


## Plot results and output
newOrder = allSigs[ order(dtMCor$value[ dtMCor$variable == "SpearmansRho" ], decreasing = TRUE) ]
dtMCor$Sig = factor(dtMCor$Sig, levels = newOrder)
dtMCor$variable = factor(dtMCor$variable, labels = c("Kendall's Tau", "Pearson's r", "Spearman's Rho"))

hLine = max(abs(dtMCor$value[ dtMCor$PAdj > 0.05 ]))

pOut = ggplot(dtMCor, aes(x = Sig, y = value, group = variable, colour = variable)) + geom_point() +
  labs(x = "Signature", y = "Correlation coefficient", colour = "CA20\ncorrelation") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  geom_hline(yintercept = 0, colour = "grey20") +
  geom_hline(yintercept = hLine, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -hLine, linetype = "dashed", colour = "grey40")

cairo_pdf(file.path(OUT, "SuppFig_24_CA20_correlation.pdf"), width = 6, height = 5)
print(pOut); dev.off()


## Simplified plot with only Kendall's tau
dtMCorSimple = dtMCor[ dtMCor$variable == "Kendall's Tau",]
pOut2 = ggplot(dtMCorSimple, aes(x = Sig, y = value)) + geom_point() +
  labs(x = "CIN Signature", y = "Kendall's Tau") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  geom_hline(yintercept = 0, colour = "grey20") +
  geom_hline(yintercept = hLine, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -hLine, linetype = "dashed", colour = "grey40") +
  coord_capped_cart(left = "both", bottom = "both")

cairo_pdf(file.path(OUT, "SuppFig_24_CA20_correlation_simple.pdf"), width = 90/25.4, height = 90/25.4)
print(pOut2); dev.off()
ggsave(file.path(OUT, "SuppFig_24_CA20_correlation_simple.svg"), pOut2, width = 90/25.4, height = 90/25.4)


write.table(dtMCorSimple, file.path(OUT, "Supp_Table_CA20.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
saveRDS(dtMCorSimple, file.path(OUT, "Supp_Table_CA20.rds"))


#### Direct comparison of Centrosome amplification and impaired HR signatures
## Find dominant IHR signature
mExp = scale(readRDS(EXP))
dfIHR = data.frame(mExp[,c("CX2", "CX3", "CX5")])
dfIHR$Class = "Neither"
dfIHR$Class[ dfIHR$CX2 > 0.5 & dfIHR$CX2 > dfIHR$CX3 & dfIHR$CX2 > dfIHR$CX5 ] = "CX2"
dfIHR$Class[ dfIHR$CX3 > 0.5 & dfIHR$CX3 > dfIHR$CX2 & dfIHR$CX3 > dfIHR$CX5 ] = "CX3"
dfIHR$Class[ dfIHR$CX5 > 0.5 & dfIHR$CX5 > dfIHR$CX2 & dfIHR$CX5 > dfIHR$CX3 ] = "CX5"
dfIHR$Class = factor(dfIHR$Class, levels = c("Neither", "CX2", "CX5", "CX3"))

## Transfer and extract signature activities
dfIHR$CA20 = dtCa$CA20[ match(substr(rownames(dfIHR),1,12), dtCa$Sample) ]


dtExp$IHR = dfIHR$Class[ match(dtExp$Var1, rownames(dfIHR)) ]
dtIHR = dtExp[ dtExp$Var2 %in% c("CX2", "CX3", "CX5"), ]
dtIHR$Var2 = factor(as.character(dtIHR$Var2), levels = c("CX2", "CX5", "CX3"))

## Plot
pOut3 = ggplot(dfIHR, aes(x = Class, y = CA20)) + geom_jitter(width = 0.2, height = 0, alpha = 0.2) +
  geom_boxplot(alpha = 0.5, outlier.colour = NA) + coord_capped_cart(left = "both", bottom = "both") + 
  labs(x = "Dominant IHR signature", y = "Centrosome amplification")

# Test significance
# dfPath %>% group_by(Pathway) %>% t_test(SVs ~ SigClass, var.equal = FALSE, p.adjust.method = "BH")
dfIHR %>% t_test(CA20 ~ Class, var.equal = FALSE, p.adjust.method = "BH")

cairo_pdf(file.path(OUT, "SuppFig_24_CA20_vs_dominant_IHR_sigs.pdf"), width = 90/25.4, height = 90/25.4)
print(pOut3); dev.off()
ggsave(file.path(OUT, "SuppFig_24_CA20_vs_dominant_IHR_sigs.svg"), pOut3, width = 90/25.4, height = 90/25.4)



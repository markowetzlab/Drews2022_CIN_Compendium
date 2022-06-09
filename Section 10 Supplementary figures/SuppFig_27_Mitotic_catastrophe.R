#### Get excess number of CNAs and their likely origin

rm(list=ls(all=TRUE))

library(this.path)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(lemon)
library(patchwork)
library(ggpubr)
library(rstatix)
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
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

PCAWGACT=file.path(BASE, "input/PCAWG_signature_activities_THRESH095_NAMESAPRIL21.rds")
META=file.path(BASE, "input/PCAWG_summary_table_combined_annotations_v4.txt")
PCAWGRAW=file.path(BASE, "input/PCAWG_1900_samples_raw_CINSig_activities_THRESH095_NAMESAPRIL21.rds")
CNAMODEL=file.path(BASE, "input/Estimated_CNAs_per_signature_and_sample_MODEL.rds")

dtPCAWG = melt(readRDS(PCAWGACT))
dtMeta = fread(META)
mRaw = readRDS(PCAWGRAW)
lmModel = readRDS(CNAMODEL)


#### Part I: CNAs vs SVs
## Plot number of CNAs vs number of SVs
#### This code is currently hacked as it needs dtHom which is generated below
dtHom = readRDS(file.path(BASE, "input/PCAWG_Microhomologies_at_SVs.rds"))
dtStat = data.table(table(dtHom$Sample))
colnames(dtStat) = c("Sample", "SVs")

## Load number of CNAs in the PCAWG cohort
dtCNAs = readRDS(file.path(BASE, "input/CNAs_PCAWG_1900samples.rds"))
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
pAB = pA + pB
ggsave(file.path(OUT, "SuppFig_27_AB_Mitotic_Catastrophe.svg"), pAB, width = 180, height = 90, units = "mm")



#### Panel C: Signature activities for each sample group
## Check CNA heavy samples for CX2,3,5
dtStatFree[["CX2"]] = scale(dtStatFree[["CX2"]])
dtStatFree[["CX3"]] = scale(dtStatFree[["CX3"]])
dtStatFree[["CX5"]] = scale(dtStatFree[["CX5"]])
dtStatMelt = melt(dtStatFree, id.vars = c("Sample", "CNAHeavy"), measure.vars = c("CX2","CX5","CX3"))
dtStatMelt$CNAHeavy = factor(dtStatMelt$CNAHeavy, levels = c("Linear", "MitoticCatastrophe"))
pC = ggplot(dtStatMelt, aes(x = variable, y = value, fill = CNAHeavy)) + geom_boxplot() +
  coord_capped_cart(left = "both", bottom = "both", ylim = c(-1, 4)) + theme(legend.position = "none") +
  labs(y = "Scaled signature activity", x = "CIN signature", fill = "CNA heavy samples") +
  scale_fill_manual(values = c("MitoticCatastrophe" = "#1f78b4", "Linear" = "#33a02c"))

ggsave(file.path(OUT, "SuppFig_27_C_Mitotic_Catastrophe.svg"), pC, width = 90, height = 90, units = "mm")





#### Calculate number of CNAs per sample
## Assume distribution of CNAs follow robust linear model
dtStatFree$ExpCNAs = predict(modelRLM)
dtStatFree$ExcessCNAs = dtStatFree$CNAs - dtStatFree$ExpCNAs

## Transfer estimated number of CNAs for each of the IHR signatures
est = mRaw %*% diag(coefficients(lmModel))
colnames(est) = colnames(mRaw)

vSumCNA = round(rowSums(est))
dtStatFree$SumCNA = vSumCNA[ match(dtStatFree$Sample, names(vSumCNA)) ]

## Transfer sig-specific estimates
dtEst = data.table(melt(est))

dtCX2 = dtEst[ dtEst$Var2 == "CX2", ]
dtCX5 = dtEst[ dtEst$Var2 == "CX5", ]
dtCX3 = dtEst[ dtEst$Var2 == "CX3", ]
dtStatFree$CX2est = dtCX2$value[ match(dtStatFree$Sample, dtCX2$Var1) ]
dtStatFree$CX5est = dtCX5$value[ match(dtStatFree$Sample, dtCX5$Var1) ]
dtStatFree$CX3est = dtCX3$value[ match(dtStatFree$Sample, dtCX3$Var1) ]

## Normalise by number of CNAs per sample
dtStatFree$CX2est = dtStatFree$CX2est/dtStatFree$CNAs
dtStatFree$CX5est = dtStatFree$CX5est/dtStatFree$CNAs
dtStatFree$CX3est = dtStatFree$CX3est/dtStatFree$CNAs

dtStatMeltCNAs = melt(dtStatFree, id.vars = c("Sample", "CNAHeavy"), 
                  measure.vars = c("CX2est","CX5est","CX3est"))
dtStatMeltCNAs$CNAHeavy = factor(dtStatMeltCNAs$CNAHeavy, levels = c("Linear", "MitoticCatastrophe"))
pD = ggplot(dtStatMeltCNAs, aes(x = variable, y = value, fill = CNAHeavy)) + geom_boxplot() +
  coord_capped_cart(left = "both", bottom = "both") + scale_y_continuous(labels = scales::percent_format()) +
  labs(y = "Proportion of CNAs", x = "CIN signature", fill = "CNA heavy samples") +
  scale_fill_manual(values = c("MitoticCatastrophe" = "#1f78b4", "Linear" = "#33a02c")) + 
  theme(legend.position = "none")

ggsave(file.path(OUT, "SuppFig_27_DMitotic_Catastrophe.svg"), pD, width = 90, height = 90, units = "mm")


## T-test
testPropCNAs <- dtStatMeltCNAs %>%
  group_by(variable) %>%
  t_test(value ~ CNAHeavy, var.equal = FALSE, p.adjust.method = "BH")




## Number of CNAs vs cell cycle score
CCS=file.path(BASE, "input/Lundberg_2020_Cell_cycle_score.rds")
dtCCS = data.table(readRDS(CCS))
dtCCS$Sample = substr(dtCCS$Sample, 1, 12)

## Categorial variable
dtStatMeltCNAs$TCGA = dtMeta$tcga_donor_barcode[ match(dtStatMeltCNAs$Sample, dtMeta$samplename) ]
dtStatMeltCNAs$CCS = dtCCS$CCS[ match(dtStatMeltCNAs$TCGA, dtCCS$Sample) ]

## Continuous variable
dtStatMeltCNAs$CCS_ct = dtCCS$CCS_ct[ match(dtStatMeltCNAs$TCGA, dtCCS$Sample) ]
dtStatMeltCNAs$CCS_ct = (dtStatMeltCNAs$CCS_ct-min(dtStatMeltCNAs$CCS_ct, na.rm = TRUE))/
  (max(dtStatMeltCNAs$CCS_ct, na.rm = TRUE)-min(dtStatMeltCNAs$CCS_ct, na.rm = TRUE))


dtCNAvsCCS = dtStatMeltCNAs[ ! is.na(dtStatMeltCNAs$CCS), ]

## According to Lundberg's script on github
dtCNAvsCCS$CCS = factor(as.numeric(dtCNAvsCCS$CCS), levels = c(1,2,3), labels = c("Low", "Medium", "High"))

pE = ggplot(dtCNAvsCCS, aes(x = CCS_ct, y = value, colour = CNAHeavy)) + 
  geom_point(alpha = 0.5) + facet_grid(. ~ variable) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) + scale_y_continuous(labels = scales::percent_format()) +
  coord_capped_cart(left = "both", bottom = "both") + theme(legend.position = "none") +
  labs(y = "Proportion of CNAs", x = "Cell cycle score", colour = "CNA heavy samples") +
  scale_colour_manual(values = c("MitoticCatastrophe" = "#1f78b4", "Linear" = "#33a02c"))

pF = ggplot(dtCNAvsCCS, aes(x = CCS, y = value, fill = CNAHeavy)) + 
  geom_boxplot() + facet_grid(. ~ variable) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) + scale_y_continuous(labels = scales::percent_format()) +
  coord_capped_cart(left = "both", bottom = "both") + theme(legend.position = "none") +
  labs(y = "Proportion of CNAs", x = "Cell cycle score", fill = "CNA heavy samples") +
  scale_fill_manual(values = c("MitoticCatastrophe" = "#1f78b4", "Linear" = "#33a02c"))


pEF = pE + pF

ggsave(file.path(OUT, "SuppFig_27_EF_Mitotic_Catastrophe.svg"), pEF, width = 180, height = 90, units = "mm")


## Stats for plot
dtCX2est = dtCNAvsCCS[dtCNAvsCCS$variable == "CX2est",]
pairwise.t.test(dtCX2est$value, paste0(dtCX2est$CNAHeavy, "_", dtCX2est$CCS), p.adjust.method = "BH", pool.sd = FALSE)

dtCX5est = dtCNAvsCCS[dtCNAvsCCS$variable == "CX5est",]
pairwise.t.test(dtCX5est$value, paste0(dtCX5est$CNAHeavy, "_", dtCX5est$CCS), p.adjust.method = "BH", pool.sd = FALSE)

dtCX3est = dtCNAvsCCS[dtCNAvsCCS$variable == "CX3est",]
pairwise.t.test(dtCX3est$value, paste0(dtCX3est$CNAHeavy, "_", dtCX3est$CCS), p.adjust.method = "BH", pool.sd = FALSE)




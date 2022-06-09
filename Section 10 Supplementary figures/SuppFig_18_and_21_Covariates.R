## Create summary of biological covariates

rm(list=ls(all=TRUE))

library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(lemon)
library(RColorBrewer)
library(stringr)
library(this.path)

theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))

BASE=dirname(this.path())
OUT=file.path(BASE, "output")
EXP=file.path(BASE, "input/Signature_Compendium_v5_Cosine-0.74_Activities_THRESH95_NAMESAPRIL21.rds")
METATCGA=file.path(BASE, "input/Metadata_TCGA_ASCAT_penalty70.rds")

## Load and prepare activities
exp = readRDS(EXP)
dtExp = data.table(melt(exp))

## Link TCGA and PCAWG
TCGAPCAWGLINK=file.path(BASE, "input/TCGA_PCAWG_links_plusCancer.rds")
link = readRDS(TCGAPCAWGLINK)

## Data sources
CTR=file.path(BASE, "input/PCAWG_chromothripsisOverlap.txt")
TDP=file.path(BASE, "input/Menghi2018_S3_TDP_inclTCGA.csv")
WGDGI=file.path(BASE, "input/Haase2019_TCGA.giScores.wgd.txt")
KAT=file.path(BASE, "input/Kataegis_PCAWG_20190130_calls_JD.txt")
CHRSIZES=file.path(BASE, "input/hg19.chrom.sizes.txt")
## Manually produced for this study. For whole chromosome CNAs, run script "SuppFig_Whole_chromosome_CNAs.R" first
ANEU=file.path(BASE, "input/Whole_chromosome_CNAs.rds")
LOH=file.path(BASE, "input/somatic_LOH_TCGA_ASCAT.rds")

## Single-base substitutions signatures
SBS=file.path(BASE, "input/PCAWG_sigProfiler_SBS_signatures_in_samples_waliqID.csv")
PCAWGACT=file.path(BASE, "input/PCAWG_signature_activities_THRESH095_NAMESAPRIL21.rds")
PCAWGSBS=file.path(BASE, "input/PCAWG_sigProfiler_SBS_signatures_in_samples.csv")
PCAWGSASBS=file.path(BASE, "input/SA_COMPOSITE_SNV.activity.FULL_SET.031918.txt")

# PCAWG metadata to link different IDs
METAPCAWG=file.path(BASE, "input/PCAWG_summary_table_combined_annotations_v4.txt")
METAPCAWG2=file.path(BASE, "input/pcawg_specimen_histology_August2016_v9.txt")


## Read data
ctr = fread(CTR)
tdp = fread(TDP)
wgdgi = fread(WGDGI)
kat = fread(KAT)
aneu = readRDS(ANEU)
loh = readRDS(LOH)

sbs = fread(SBS)
pcawgSBS = fread(PCAWGSBS)
pcawgAct = readRDS(PCAWGACT)
pcawgSASBS = fread(PCAWGSASBS)

meta = readRDS(METATCGA)
dtMeta = fread(METAPCAWG)
dtMeta2 = fread(METAPCAWG2)


### Functions
twoMatCorr = function(MAT1, MAT2, OUTTAB, OUTFIG, AXISLABEL) {
  
  ## Adjust order
  # All samples for which we have sigs are in the SBS data
  if(sum(! rownames(MAT1) %in% rownames(MAT2)) == 0) {
    MAT2Filt = MAT2[ match(rownames(MAT1), rownames(MAT2)), ]
    identical(rownames(MAT1), rownames(MAT2Filt))
  } else { 
    warning("Rownames of matrix 1 not all present in matrix 2.")
    
    ## Filter both matrices to contain the same samples
    MAT2Filt = MAT2[ rownames(MAT2) %in% rownames(MAT1), ]
    MAT1 = MAT1[ rownames(MAT1) %in% rownames(MAT2), ]
    
    ## Order each according to one of the matrices
    MAT2Filt = MAT2Filt[ match(rownames(MAT1), rownames(MAT2Filt)), ]
    
    identical(rownames(MAT1), rownames(MAT2Filt))
  }
  
  
  # Calculate p-values for the correlations
  allSigs = colnames(MAT1)
  allSBS = colnames(MAT2Filt)
  lSpear2 = lapply(allSigs, function(thisCS) {
    
    lSBS = lapply(allSBS, function(thisSBS) {
      corTest = cor.test(MAT2Filt[, thisSBS], MAT1[, thisCS], method = "spearman", 
                         exact = FALSE)
      out = c(thisCS, thisSBS, signif(corTest$estimate, 4), signif(corTest$p.value, 4))
      return(out)
    })
    
    dtSBS = data.table(do.call(rbind, lSBS))
    colnames(dtSBS) = c("CSsig", "SBSsig", "Spearman", "pVal")
    dtSBS$Spearman = as.numeric(dtSBS$Spearman)
    dtSBS$pVal = as.numeric(dtSBS$pVal)
    dtSBS$pAdj = p.adjust(dtSBS$pVal, method = "BH")
    
    return(dtSBS)
    
  })
  
  dtSpear2 = rbindlist(lSpear2)
  
  # Save output and plot
  write.table(dtSpear2, OUTTAB, sep = "\t", 
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  # Plot only significant relationships
  dtSpear2$Spearman[ dtSpear2$pAdj > 0.05 ] = NA
  dtSpear2$CSsig = factor(dtSpear2$CSsig, levels = allSigs)
  dtSpear2$SBSsig = factor(dtSpear2$SBSsig, levels = allSBS)
  pSBS2 = ggplot(dtSpear2, aes(x = SBSsig, y = CSsig, fill = Spearman)) + geom_tile() +
    scale_x_discrete(drop = FALSE, guide = guide_axis(angle = 90)) + scale_y_discrete(drop = FALSE) +
    theme(legend.position = "bottom", axis.text = element_text(size = 7), 
          aspect.ratio = ncol(MAT1)/ncol(MAT2Filt)) +
    scale_fill_gradient2(na.value = "white") + xlab(AXISLABEL) + ylab("CIN signatures") +
    labs(fill = "Spearman\ncorrelation") + coord_capped_cart(left = "both", bottom = "both")
  
  # cairo_pdf(paste0(OUTFIG, ".pdf"), width = 180/25.4, height = 2.5)
  # print(pSBS2); dev.off()
  ggsave(paste0(OUTFIG, ".svg"), pSBS2, width = 180/25.4, height = 2.5)
  
  return(pSBS2)
}





### Chromothripsis
ctr=data.table(table(ctr$samplename, ctr$FinalCalls))
ctr=ctr[ ctr$V2 == "Chromothripsis", ]
ctr$sample = link$TCGA[ match(ctr$V1, link$PCAWG) ]
dtExp$ctr = ctr$N[ match(substr(dtExp$Var1,1,12), ctr$sample) ]



### Tandem-duplication (score)
tdp = tdp[ tdp$Study_Name == "TCGA", ]
dtExp$tdpscore = tdp$TDP_score[ match(substr(dtExp$Var1,1,12), 
                                           substr(tdp$Sample_ID, 1, 12)) ]
### Tandem-duplications (status)
dtExp$tdpstatus = tdp$TDP_status[ match(substr(dtExp$Var1,1,12), 
                                             substr(tdp$Sample_ID, 1, 12)) ]
dtExp$tdpstatus = as.logical(dtExp$tdpstatus)



### WGD = whole genome duplication
dtExp$wgd = wgdgi$wgd[ match(substr(dtExp$Var1,1,12), wgdgi$patient) ]



### Kataegis
kat$TCGA = link$TCGA[ match(kat$sample, link$PCAWG) ]
katFilt = kat[ ! is.na(kat$TCGA) ]
# Get only significant cluster
katFilt = katFilt[ katFilt$p_streak_adj < 0.05, ]
# Transfer kataegis status. TRUE - in list, FALSE - in PCAWG and TCGA but not in list, NA - for everyone else
dtExp$Kataegis = NA
dtExp$Kataegis[ dtExp$Var1 %in% link$TCGA ] = FALSE
dtExp$Kataegis[ dtExp$Var1 %in% unique(katFilt$TCGA) ] = TRUE
# Summarise kataegis data
dfKatFilt = data.frame(table(katFilt$TCGA))
dtExp$KatNum = dfKatFilt$Freq[ match(dtExp$Var1, dfKatFilt$Var1) ]
dtExp$KatNum[ dtExp$Var1 %in% link$TCGA & ! dtExp$Var1 %in% dfKatFilt$Var1 ] = 0


## Aneuploidy - whole chromosome CNAs
dtExp$Aneu = aneu$whole[ match(dtExp$Var1, aneu$sample) ]



## LOH
# Estimate human genome length
hgLength = fread(CHRSIZES)
HGL = sum(hgLength$V2[1:24])
## Two measures: Number of LOH and proportion of LOH across the genome
dtNLoh = data.table(table(loh$sample))

loh$width = width(loh)
dtWidth = data.table(aggregate(width ~ sample, data = loh, sum))
dtWidth$Num = dtNLoh$N[ match(dtNLoh$V1, dtWidth$sample) ]

dtWidth$Prop = dtWidth$width / HGL
dtWidth$Cancer = meta$cancer_type[ match(dtWidth$sample, meta$name) ]

dtExp$LOH = dtWidth$Prop[ match(dtExp$Var1, dtWidth$sample) ]


## Loop over samples and test
allSigs = levels(dtExp$Var2)
wgdCohortProp = signif(mean(dtExp$wgd[ dtExp$Var2 == "CX1" ]), 4)
lTests = lapply(allSigs, function(thisSig) {
  
  thisExp = dtExp[ dtExp$Var2 == thisSig, ]
  
  ## Continuous vars
  corCTR = cor.test(thisExp$value, thisExp$ctr, method = 'spearman', use = 'pairwise.complete.obs')
  corTDPscore = cor.test(thisExp$value, thisExp$tdpscore, method = 'spearman', use = 'pairwise.complete.obs')
  corKatNum = cor.test(thisExp$value, thisExp$KatNum, method = 'spearman', use = 'pairwise.complete.obs')
  corAneu = cor.test(thisExp$value, thisExp$Aneu, method = 'spearman', use = 'pairwise.complete.obs')
  corLOH = cor.test(thisExp$value, thisExp$LOH, method = 'spearman', use = 'pairwise.complete.obs')
  
  ## Categorial vars
  tTDPstatus = t.test(thisExp$value[thisExp$tdpstatus], thisExp$value[! thisExp$tdpstatus], var.equal = FALSE)
  tKat = t.test(thisExp$value[thisExp$Kataegis], thisExp$value[! thisExp$Kataegis], var.equal = FALSE)
  
  ## Whole genome duplication
  # Proportion test at 0
  zeros = thisExp[ thisExp$value == 0, ]
  pWGD = prop.test(x = sum(zeros$wgd), n = nrow(zeros), p = wgdCohortProp)
  wgd = signif(mean(zeros$wgd), 4)
  # t-Test above zero exposures
  tWGD = t.test(thisExp$value[thisExp$wgd & thisExp$value != 0], 
                thisExp$value[(! thisExp$wgd) & thisExp$value != 0 ], var.equal = FALSE)
  
  ## Prepare output
  out = c(corCTR$estimate, corCTR$p.value, 
        corTDPscore$estimate, corTDPscore$p.value,
        corKatNum$estimate, corKatNum$p.value,
        corAneu$estimate, corAneu$p.value,
        corLOH$estimate, corLOH$p.value,
        tTDPstatus$estimate[1]-tTDPstatus$estimate[2], tTDPstatus$p.value,
        tKat$estimate[1]-tKat$estimate[2], tKat$p.value,
        pWGD$estimate - wgdCohortProp, pWGD$p.value,
        tWGD$estimate[1]-tWGD$estimate[2], tWGD$p.value )
  out = signif(out, 3)
  return(out)
  
})

dtTests = data.table("Sig" = allSigs, do.call(rbind, lTests))
colnames(dtTests) = c("Signature", "CTRcor", "CTRpval", "TDPscorecor", "TDPscorepval", 
                      "NumKatcor", "NumKatpval", "Aneucor", "Aneupval", "LOHcor", "LOHpval",
                      "TDPstatusdiff", "TDPstatuspval", "Katdiff", "Katpval", 
                      "WGDpropdiff", "WGDproppval", "WGDvaldiff", "WGDvalpval")

## Correct all p-values together
allPvals = c(dtTests$CTRpval, dtTests$TDPscorepval, dtTests$NumKatpval, dtTests$Aneupval, 
             dtTests$LOHpval, dtTests$TDPstatuspval, dtTests$Katpval,
             dtTests$WGDproppval, dtTests$WGDvalpval)
allPadj = p.adjust(allPvals, method = "BH")

dtTests$CTRpadj = allPadj[1:17]
dtTests$TDPscorepadj = allPadj[18:34]
dtTests$NumKatpadj = allPadj[35:51]
dtTests$Aneupadj = allPadj[52:68]
dtTests$LOHpadj = allPadj[69:85]
dtTests$TDPstatuspadj = allPadj[86:102]
dtTests$Katpadj = allPadj[103:119]
dtTests$WGDproppadj = allPadj[120:136]
dtTests$WGDvalpadj = allPadj[137:153]


## Write output for Matt and shiny app
lOut = list("Chromothripsis" = data.table("Signature" = dtTests$Signature, 
                                   "Covariate" = "Chromothripsis", 
                                   "Correlation" = dtTests$CTRcor, 
                                   "pAdj" = dtTests$CTRpadj),
     "Kataegis" = data.table("Signature" = dtTests$Signature, 
                             "Covariate" = "Kataegis", 
                             "Correlation" = dtTests$NumKatcor, 
                             "pAdj" = dtTests$NumKatpadj),
     "LOH" = data.table("Signature" = dtTests$Signature, 
                             "Covariate" = "LOH", 
                             "Correlation" = dtTests$LOHcor, 
                             "pAdj" = dtTests$LOHpadj),
     "WholeCNA" = data.table("Signature" = dtTests$Signature, 
                        "Covariate" = "Whole-chromosome CNAs", 
                        "Correlation" = dtTests$Aneucor, 
                        "pAdj" = dtTests$Aneupadj),
     "WGDProp" = data.table("Signature" = dtTests$Signature, 
                             "Covariate" = "Whole-genome duplication (proportion test)", 
                             "Difference" = dtTests$WGDpropdiff, 
                             "pAdj" = dtTests$WGDproppadj),
     "WGDtTest" = data.table("Signature" = dtTests$Signature, 
                            "Covariate" = "Whole-genome duplication (t-test)", 
                            "Difference" = dtTests$WGDvaldiff, 
                            "pAdj" = dtTests$WGDvalpadj))
saveRDS(lOut, file.path(OUT, "Covariates_CtrKatLOHAneuWGD.rds"))

# Save output
write.table(dtExp, file.path(OUT, "SM_8_Signature_activities_and_biological_covariates.txt"),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(dtTests, file.path(OUT, "SM_8_Tests_biological_covariates.txt"),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

## Plot output
## Convert stupid table format
dtTestsMelt = melt(dtTests, id.vars = "Signature")

dtTestsMeltVals = dtTestsMelt[ ! (grepl("padj",dtTestsMelt$variable) | 
                                    grepl("pval",dtTestsMelt$variable)), ]
colnames(dtTestsMeltVals) = c("Signature", "Covariate", "Correlation")

dtTestsMeltAdj = dtTestsMelt[ grep("padj",dtTestsMelt$variable), ]
colnames(dtTestsMeltAdj) = c("Sig", "CovPAdj", "PAdj")

dtTestsPlot = cbind(dtTestsMeltVals, dtTestsMeltAdj)
dtTestsPlot$Sig = NULL

dtTestsPlot$PlotPadj = -log(dtTestsPlot$PAdj)
dtTestsPlot$Significant = dtTestsPlot$PAdj < 0.05


## Not needed as plotted individually
# pOut = ggplot(dtTestsPlot, aes(x = Correlation, y = PlotPadj, colour = Significant)) +
#   facet_wrap(Covariate ~ ., scales = "free", ncol = 2) + 
#   geom_vline(xintercept = 0, colour = "grey30", linetype = "dashed") + 
#   geom_point() + geom_text_repel(aes(label = Signature))


# cairo_pdf(file.path(OUTFIGURES, "EDF_xx_Covariates_table.pdf"), width = 290/25.4, height = 15)
# print(pOut); dev.off()
# ggsave(file.path(OUTPLOTS, "Figure_3_Aetiologies.svg"), pOut2, width = 67.7/25.4, height = 117.5/25.4)


## Chromothripsis
pCTR = ggplot(dtTestsPlot[ dtTestsPlot$Covariate == "CTRcor", ], 
       aes(x = Correlation, y = PlotPadj, colour = Significant)) + 
  geom_point() + geom_text_repel(aes(label = Signature), force = 100) + 
  coord_capped_cart(left = "both", bottom = "both") +
  labs(x = "Spearman's Rho", y = "-log(q-value)", title = "Chromothripsis") + 
  geom_vline(xintercept = 0, colour = "grey30", linetype = "dashed") +
  theme(legend.position = c(0.3, 0.8)) +
  scale_colour_manual(values = c("grey80", "black"))

# cairo_pdf(file.path(OUT, "SuppFig_Covariates_Chromothripsis.pdf"), width = 90/25.4, height = 90/25.4)
# print(pCTR); dev.off()
ggsave(file.path(OUT, "SuppFig_18_Covariates_Chromothripsis.svg"), pCTR, width = 90/25.4, height = 90/25.4)


## Kataegis
pKAT = ggplot(dtTestsPlot[ dtTestsPlot$Covariate == "NumKatcor", ], 
              aes(x = Correlation, y = PlotPadj, colour = Significant)) + 
  geom_point() + geom_text_repel(aes(label = Signature), force = 100) + 
  coord_capped_cart(left = "both", bottom = "both") +
  labs(x = "Spearman's Rho", y = "-log(q-value)", title = "Kataegis") + 
  geom_vline(xintercept = 0, colour = "grey30", linetype = "dashed") +
  theme(legend.position = c(0.3, 0.8)) +
  scale_colour_manual(values = c("grey80", "black"))

# cairo_pdf(file.path(OUT, "SuppFig_Covariates_Kataegies.pdf"), width = 90/25.4, height = 90/25.4)
# print(pKAT); dev.off()
ggsave(file.path(OUT, "SuppFig_18_Covariates_Kataegies.svg"), pKAT, width = 90/25.4, height = 90/25.4)


## LOH
pLOH = ggplot(dtTestsPlot[ dtTestsPlot$Covariate == "LOHcor", ], 
              aes(x = Correlation, y = PlotPadj, colour = Significant)) + 
  geom_point() + geom_text_repel(aes(label = Signature), force = 100) + 
  coord_capped_cart(left = "both", bottom = "both") +
  labs(x = "Spearman's Rho", y = "-log(q-value)", title = "Loss of heterozygosity") + 
  geom_vline(xintercept = 0, colour = "grey30", linetype = "dashed") +
  theme(legend.position = c(0.3, 0.8)) +
  scale_colour_manual(values = c("grey80", "black"))

# cairo_pdf(file.path(OUT, "SuppFig_Covariates_LOH.pdf"), width = 90/25.4, height = 90/25.4)
# print(pLOH); dev.off()
ggsave(file.path(OUT, "SuppFig_18_Covariates_LOH.svg"), pLOH, width = 90/25.4, height = 90/25.4)


## Whole-chromosome CNAs
pAN = ggplot(dtTestsPlot[ dtTestsPlot$Covariate == "Aneucor", ], 
              aes(x = Correlation, y = PlotPadj, colour = Significant)) + 
  geom_point() + geom_text_repel(aes(label = Signature), force = 100) + 
  coord_capped_cart(left = "both", bottom = "both") +
  labs(x = "Spearman's Rho", y = "-log(q-value)", title = "Whole-chromosome CNAs") + 
  geom_vline(xintercept = 0, colour = "grey30", linetype = "dashed") +
  theme(legend.position = c(0.3, 0.8)) +
  scale_colour_manual(values = c("grey80", "black"))

# cairo_pdf(file.path(OUT, "SuppFig_Covariates_Whole-chr_CNAs.pdf"), width = 90/25.4, height = 90/25.4)
# print(pAN); dev.off()
ggsave(file.path(OUT, "SuppFig_18_Covariates_Whole-chr_CNAs.svg"), pAN, width = 90/25.4, height = 90/25.4)



## WGD - val diff / prop diff
# Proportion at zero
pWGD1 = ggplot(dtTestsPlot[ dtTestsPlot$Covariate == "WGDpropdiff", ], 
              aes(x = Correlation, y = PlotPadj, colour = Significant)) + 
  geom_point() + geom_text_repel(aes(label = Signature), force = 100) + 
  coord_capped_cart(left = "both", bottom = "both") +
  labs(x = "Proportion of WGD samples at zero compared to background", y = "-log(q-value)", title = "Whole-genome duplication: Proportion test at zero") + 
  geom_vline(xintercept = 0, colour = "grey30", linetype = "dashed") +
  theme(legend.position = c(0.3, 0.8)) +
  scale_colour_manual(values = c("grey80", "black"))

# cairo_pdf(file.path(OUT, "SuppFig_Covariates_WGD1.pdf"), width = 90/25.4, height = 90/25.4)
# print(pWGD1); dev.off()
ggsave(file.path(OUT, "SuppFig_18_Covariates_WGD1.svg"), pWGD1, width = 90/25.4, height = 90/25.4)


# Difference for values above zero
pWGD2 = ggplot(dtTestsPlot[ dtTestsPlot$Covariate == "WGDvaldiff", ], 
               aes(x = Correlation, y = PlotPadj, colour = Significant)) + 
  geom_point() + geom_text_repel(aes(label = Signature), force = 100) + 
  coord_capped_cart(left = "both", bottom = "both") +
  labs(x = "Difference in activities between WGD and non-WGD samples", y = "-log(q-value)", 
       title = "Whole-genome duplication: T test for samples with activities above zero") + 
  geom_vline(xintercept = 0, colour = "grey30", linetype = "dashed") +
  theme(legend.position = c(0.3, 0.8)) +
  scale_colour_manual(values = c("grey80", "black"))

# cairo_pdf(file.path(OUT, "SuppFig_Covariates_WGD2.pdf"), width = 90/25.4, height = 90/25.4)
# print(pWGD2); dev.off()
ggsave(file.path(OUT, "SuppFig_18_Covariates_WGD2.svg"), pWGD2, width = 90/25.4, height = 90/25.4)




# ### SBS signatures
# 
# # Transfer TCGA names to PCAWG SBS exposures
# sbs$TCGA = link$TCGA[ match(sbs$`Sample Name`, link$PCAWG) ]
# sbsFilt = sbs[ sbs$TCGA %in% substr(rownames(exp),1,12), ]
# sbsNames = sbsFilt$TCGA
# sbsFilt$TCGA = NULL
# sbsFilt = sbsFilt[,-(1:3)]
# mSBS = as.matrix(sbsFilt)
# rownames(mSBS) = sbsNames
# 
# mSBS = t(apply(mSBS,1,function(x) x/sum(x)))
# 
# mExp = exp[ match(rownames(mSBS), substr(rownames(exp),1,12)), ]
# 
# identical(rownames(mExp), rownames(mSBS))
# 
# # The Magic
# corSpearman = cor(mSBS, mExp, method = 'spearman')
# dfSpear = melt(corSpearman)
# 
# # Calculate p-values for the correlations
# dfSpear$pVals = apply(dfSpear, 1, function(thisRow) {
#   
#   if(is.na(thisRow[3])) { return(NA) }
#   thisSBS = as.character(thisRow[1])
#   thisCS = as.character(thisRow[2])
#   suppressWarnings(cor.test(mSBS[, thisSBS], mExp[, thisCS], method = "spearman")$p.value)
#   
# })
# 
# dfSpear$pAdj = p.adjust(dfSpear$pVals, method = "BH")
# colnames(dfSpear) = c("SBSsig", "CSsig", "Spearman", "pVal", "pAdj")
# 
# # Save output and plot
# write.table(dfSpear, file.path(OUTTABLES, "EDF_xx_SBS_Correlation.txt"), sep = "\t", 
#             row.names = FALSE, col.names = TRUE, quote = FALSE)
# 
# # Plot only significant relationships
# dfSpear$Spearman[ dfSpear$pAdj > 0.05 ] = NA
# pSBS = ggplot(dfSpear, aes(x = SBSsig, y = CSsig, fill = Spearman)) + geom_tile() +
#   scale_x_discrete(drop = FALSE, guide = guide_axis(angle = 90)) + scale_y_discrete(drop = FALSE) +
#   theme(legend.position = "bottom", axis.text = element_text(size = 7), aspect.ratio = 17/65) +
#   scale_fill_gradient2(na.value = "white") + xlab("SBS signature") + ylab("CIN signatures") +
#   labs(fill = "Spearman\ncorrelation") + coord_capped_cart(left = "both", bottom = "both")
# 
# cairo_pdf(file.path(OUT, "EDF_xx_SBS_Signatures.pdf"), width = 180/25.4, height = 2.5)
# print(pSBS); dev.off()
# ggsave(filename = file.path(OUT, "EDF_xx_SBS_Signatures.svg"), pSBS,
#        width = 180/25.4, height = 2.5)




###### Redo SBS with all of PCAWG

#### SignatureProfiler first
## First: Map names from sample names (08fgfx-...) to specimen IDs (SPxxx)
dtIDs = data.table("barcodes" = rownames(pcawgAct))
dtIDs$SA = dtMeta$icgc_sample_id[ match(dtIDs$barcodes, dtMeta$samplename) ]
dtIDs$SP = dtMeta2$`# icgc_specimen_id`[ match(dtIDs$SA, dtMeta2$icgc_sample_id)]

## Replace names
rownames(pcawgAct) = dtIDs$SP[ match(rownames(pcawgAct), dtIDs$barcodes) ]

## Prepare SBS sigs
mSBSPCAWG = as.matrix(pcawgSBS[,4:ncol(pcawgSBS)])
rownames(mSBSPCAWG) = pcawgSBS$`Sample Names`
mSBSPCAWG = t(apply(mSBSPCAWG,1,function(x) x/sum(x)))

pSP_SBS = twoMatCorr(MAT1 = pcawgAct, MAT2 = mSBSPCAWG, 
                     OUTTAB = file.path(OUT, "SuppFig_21_SBS_Correlation_allPCAWG_SP.txt"), 
                     OUTFIG = file.path(OUT, "SuppFig_21_SBS_Signatures_allPCAWG_SP"), 
                     AXISLABEL = "SP SBS signatures")




#### SP: DBS Sigs
SPDBS = file.path(BASE, "input/PCAWG_sigProfiler_DBS_signatures_in_samples.csv")
dtSPDBS = fread(SPDBS)

mSPDBS = as.matrix(dtSPDBS[,4:ncol(dtSPDBS)])
rownames(mSPDBS) = dtSPDBS$`Sample Names`
mSPDBS = t(apply(mSPDBS,1,function(x) x/sum(x)))

pSP_DBS = twoMatCorr(MAT1 = pcawgAct, MAT2 = mSPDBS, 
                     OUTTAB = file.path(OUT, "SuppFig_21_DBS_Correlation_allPCAWG_SP.txt"), 
                     OUTFIG = file.path(OUT, "SuppFig_21_DBS_Signatures_allPCAWG_SP"), 
                     AXISLABEL = "SP DBS signatures")


#### SP: IND Sigs
SPID = file.path(BASE, "input/PCAWG_sigProfiler_ID_signatures_in_samples.csv")
dtSPID = fread(SPID)

mSPID = as.matrix(dtSPID[,4:ncol(dtSPID)])
rownames(mSPID) = dtSPID$`Sample Names`
mSPID = t(apply(mSPID,1,function(x) x/sum(x)))

pSP_ID = twoMatCorr(MAT1 = pcawgAct, MAT2 = mSPID, 
                     OUTTAB = file.path(OUT, "SuppFig_21_ID_Correlation_allPCAWG_SP.txt"), 
                     OUTFIG = file.path(OUT, "SuppFig_21_ID_Signatures_allPCAWG_SP"), 
                     AXISLABEL = "SP ID signatures")



#### SA SBS
## Fix sample names
colnames(pcawgSASBS) = sapply(str_split(colnames(pcawgSASBS), "__"), function(x) x[2])
theseRowNames = pcawgSASBS[,"V1"]
pcawgSASBS[,"V1"] = NULL

## Convert to matrix and transpose
mSA = t(as.matrix(pcawgSASBS))

## Add simple Signature names
colnames(mSA) = sapply(str_split(theseRowNames$V1, "\\_"), function(x) x[4])
          
## Normalise
mSA = t(apply(mSA,1,function(x) x/sum(x)))

## Correlate and plot
pSA_SBS = twoMatCorr(MAT1 = pcawgAct, MAT2 = mSA, 
                     OUTTAB = file.path(OUT, "SuppFig_21_SBS_Correlation_allPCAWG_SA.txt"),
                     OUTFIG = file.path(OUT, "SuppFig_21_SBS_Signatures_allPCAWG_SA"),
                     AXISLABEL = "SA SBS signatures")



#### SA DBS
SADBS = file.path(BASE, "input/SignatureAnalyzer_DBS.activity.FULL_SET.012718.txt")
dtSADBS = fread(SADBS)

## Fix sample names
colnames(dtSADBS) = sapply(str_split(colnames(dtSADBS), "__"), function(x) x[2])
theseRowNames = dtSADBS[,"V1"]
dtSADBS[,"V1"] = NULL

## Convert to matrix and transpose
mSADBS = t(as.matrix(dtSADBS))

## Add simple Signature names
colnames(mSADBS) = sapply(str_split(theseRowNames$V1, "\\_"), function(x) x[2])

## Normalise
mSADBS = t(apply(mSADBS,1,function(x) x/sum(x)))

## Correlate and plot
pSA_DBS = twoMatCorr(MAT1 = pcawgAct, MAT2 = mSADBS, 
                     OUTTAB = file.path(OUT, "SuppFig_21_DBS_Correlation_allPCAWG_SA.txt"),
                     OUTFIG = file.path(OUT, "SuppFig_21_DBS_Signatures_allPCAWG_SA"),
                     AXISLABEL = "SA DBS signatures")



#### SA DBS
SAID = file.path(BASE, "input/SignatureAnalyzer_ID.activity.FULL_SET.012718.txt")
dtSAID = fread(SAID)

## Fix sample names
colnames(dtSAID) = sapply(str_split(colnames(dtSAID), "__"), function(x) x[2])
theseRowNames = dtSAID[,"V1"]
dtSAID[,"V1"] = NULL

## Convert to matrix and transpose
mSAID = t(as.matrix(dtSAID))

## Add simple Signature names
colnames(mSAID) = sapply(str_split(theseRowNames$V1, "\\_"), function(x) x[2])

## Normalise
mSADBS = t(apply(mSAID,1,function(x) x/sum(x)))

## Correlate and plot
pSA_ID = twoMatCorr(MAT1 = pcawgAct, MAT2 = mSAID, 
                     OUTTAB = file.path(OUT, "SuppFig_21_ID_Correlation_allPCAWG_SA.txt"),
                     OUTFIG = file.path(OUT, "SuppFig_21_ID_Signatures_allPCAWG_SA"),
                     AXISLABEL = "SA ID signatures")


          
          
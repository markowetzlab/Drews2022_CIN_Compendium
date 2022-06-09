# Compare activities for the same samples from different high-throughput technologies

rm(list = ls(all=TRUE))

library(this.path)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(lsa)
library(patchwork)
library(lemon)
library(RColorBrewer)

theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))

## Function
plotDefinitionComparisons = function(mWithNorm, vFiles, vLevels, vLabels, OUT1) {
  
  lNoNorm = list()
  for(i in 1:length(vFiles)) {
    
    theseSigs = fread(vFiles[i])
    vFeats = theseSigs$V1
    mSigs = t(as.matrix(theseSigs[,-1]))
    colnames(mSigs) = vFeats
    rownames(mSigs) = paste0(rownames(mSigs), "_", names(vFiles[i]))
    lNoNorm[[length(lNoNorm)+1]] = mSigs
  }
  names(lNoNorm) = names(vFiles)
  
  ### Compare cosine similarities
  lCos = list()
  for(i in 1:length(lNoNorm)) {
    
    thisNoNorm = lNoNorm[[i]]
    mWithNorm = mWithNorm[ ,match(colnames(thisNoNorm), colnames(mWithNorm)) ]
    
    vCos = sapply(1:nrow(mWithNorm), function(i) { 
      max(cosine(mWithNorm[i,], t(thisNoNorm)))
    })
    dtCos = data.table("sig" = rownames(mWithNorm), "maxcos" = vCos, "penalty" = names(lNoNorm[i]))
    lCos[[length(lCos) + 1]] = dtCos
    
  }
  
  dtCos = rbindlist(lCos)
  dtCos$penalty = factor(dtCos$penalty, 
                         levels = vLevels, 
                         labels = vLabels)
  dtCos$sig = factor(dtCos$sig, levels = rev(unique(dtCos$sig)))
  
  # Categorise max cosine similarity
  dtCos$category = "<0.5"
  dtCos$category[ dtCos$maxcos >= 0.5 & dtCos$maxcos < 0.8 ] = "<0.8"
  dtCos$category[ dtCos$maxcos >= 0.8 & dtCos$maxcos < 0.95 ] = "<0.95"
  # dtCos$category[ dtCos$maxcos >= 0.9 & dtCos$maxcos < 0.95 ] = "<0.95"
  dtCos$category[ dtCos$maxcos >= 0.95 ] = ">0.95"
  dtCos$category = factor(dtCos$category, levels = c("<0.5", "<0.8", "<0.95", ">0.95"))
  
  p1 = ggplot(dtCos, aes(x = penalty, y = sig, fill = category)) + geom_tile() +
    scale_x_discrete(drop=FALSE) + xlab("ASCAT Penalty") + ylab("Gold standard pan-cancer sigs") + 
    labs(fill = "Max cosine\nsimilarity") +
    # scale_fill_gradient2(low = "white", mid = "white", midpoint = 0.5, high = "#cb181d", limits = c(0,1))
    scale_fill_manual(values = c("#ffffcc", "#fed976", "#fd8d3c", "#e31a1c"), drop = FALSE)
  
  cairo_pdf(OUT1, height = 4, width = 5); print(p1); dev.off()

  return(p1)
}


## Paths
BASE=dirname(this.path())

# Load gold standard (comparison signature definitions)
mWithNorm = readRDS(file.path(BASE, "input/Definitions/Pancancer_signatures_Oct2020.rds"))


## EDF 6: WES on-target
OUT1=file.path(BASE, "output/EDF_6_bottomrow_Definitions_cel478_vs_WESontarget_Koptimal.pdf")
vFiles = c("penalty35" = file.path(BASE, "input/Definitions/6_Signatures_UI511A.txt"),
           "penalty50" = file.path(BASE, "input/Definitions/6_Signatures_5DC58S.txt"),
           "penalty70" = file.path(BASE, "input/Definitions/6_Signatures_I2YCXE.txt"),
           "penalty100" = file.path(BASE, "input/Definitions/6_Signatures_H56S47.txt"),
           "penalty140" = file.path(BASE, "input/Definitions/6_Signatures_SWWY4F.txt"))
vLevels = names(vFiles)
vLabels = c(35, 50, 70, 100, 140)

catch1 = plotDefinitionComparisons(mWithNorm, vFiles, vLevels, vLabels, OUT1)



## EDF 6: WES off-target
OUT1=file.path(BASE, "output/EDF_6_bottomrow_Definitions_cel478_vs_WESofftarget_Koptimal.pdf")
vFiles = c("alpha0.001" = file.path(BASE, "input/Definitions/6_Signatures_303RA8.txt"),
           "alpha0.01" = file.path(BASE, "input/Definitions/6_Signatures_3OUSCT.txt"),
           "alpha0.05" = file.path(BASE, "input/Definitions/6_Signatures_0X801D.txt"))
vLevels = names(vFiles)
vLabels = c(0.001, 0.01, 0.05)

catch2 = plotDefinitionComparisons(mWithNorm, vFiles, vLevels, vLabels, OUT1)



## EDF 6: shallow WGS
OUT1=file.path(BASE, "output/EDF_6_bottomrow_Definitions_cel478_vs_sWGS_K10.pdf")
vFiles = c("alpha0.001" = file.path(BASE, "input/Definitions/6_Signatures_3SEUB9.txt"),
           "alpha0.01" = file.path(BASE, "input/Definitions/6_Signatures_6ZWIJT.txt"),
           "alpha0.05" = file.path(BASE, "input/Definitions/6_Signatures_GI0J7B.txt"))
vLevels = names(vFiles)
vLabels = c(0.001, 0.01, 0.05)

catch3 = plotDefinitionComparisons(mWithNorm, vFiles, vLevels, vLabels, OUT1)



## EDF 6: WGS downsampled to SNP6 positions
OUT1=file.path(BASE, "output/EDF_6_bottomrow_Definitions_cel478_vs_wgssnp6_K10.pdf")
vFiles = c("penalty50" = file.path(BASE, "input/Definitions/6_Signatures_E8610I.txt"),
           "penalty70" = file.path(BASE, "input/Definitions/6_Signatures_G25CO5.txt"),
           "penalty100" = file.path(BASE, "input/Definitions/6_Signatures_GTGB8P.txt"),
           "penalty140" = file.path(BASE, "input/Definitions/6_Signatures_FAQTWZ.txt"))
# 35 had no solution with K=10
vLevels = c("penalty35", "penalty50", "penalty70", "penalty100", "penalty140")
vLabels = c(35, 50, 70, 100, 140)

catch4 = plotDefinitionComparisons(mWithNorm, vFiles, vLevels, vLabels, OUT1)



## EDF 6: SNP6 without matched normal
OUT1=file.path(BASE, "output/EDF_6_bottomrow_Definitions_cel478_vs_cel478NoNorm_K10_ASCATcl.pdf")
vFiles = c("penalty35" = file.path(BASE, "input/Definitions/6_Signatures_K0KUA5.txt"),
           "penalty50" = file.path(BASE, "input/Definitions/6_Signatures_A2H2K6.txt"),
           "penalty70" = file.path(BASE, "input/Definitions/6_Signatures_5N498Y.txt"),
           "penalty100" = file.path(BASE, "input/Definitions/6_Signatures_K6DGUK.txt"),
           "penalty140" = file.path(BASE, "input/Definitions/6_Signatures_HJUVT9.txt"))
vLevels = c("penalty35", "penalty50", "penalty70", "penalty100", "penalty140")
vLabels = c(35, 50, 70, 100, 140)

catch5 = plotDefinitionComparisons(mWithNorm, vFiles, vLevels, vLabels, OUT1)



## SuppFig 33: WGS all data
OUT1=file.path(BASE, "output/SuppFig_33_c_Definitions_cel478_vs_WGS_Full_data_Koptimal.pdf")
vFiles = c("ABSOLUTE" = file.path(BASE, "input/Definitions/6_Signatures_GWAUZB.txt"),
           "Battenberg" = file.path(BASE, "input/Definitions/6_Signatures_A6DFCQ.txt"),
           "Sclust" = file.path(BASE, "input/Definitions/6_Signatures_ICEU3T.txt"),
           "Consensus" = file.path(BASE, "input/Definitions/6_Signatures_OFK7V5.txt"),
           "ASCAT SNP6 pos." = file.path(BASE, "input/Definitions/6_Signatures_G25CO5.txt"))
vLevels = names(vFiles)
vLabels = names(vFiles)

catch6 = plotDefinitionComparisons(mWithNorm, vFiles, vLevels, vLabels, OUT1)



## SuppFig 33: WES all data (two callers)
OUT1=file.path(BASE, "output/SuppFig_33_c_Definitions_cel478_vs_TCGA_WES_Two_callers.pdf")
vFiles = c("ASCAT70" = file.path(BASE, "input/Definitions/6_Signatures_I2YCXE.txt"),
           "Sequenca" = file.path(BASE, "input/Definitions/6_Signatures_A8O3B7.txt"),
           "CNVkit_Standard" = file.path(BASE, "input/Definitions/6_Signatures_W61UFI.txt"),
           "CNVkit_Refitted" = file.path(BASE, "input/Definitions/6_Signatures_T2IIIX.txt"))
vLevels = c("ASCAT70", "Sequenca", "CNVkit_Standard", "CNVkit_Refitted")
vLabels = c("ASCAT70", "Sequenca", "CNVkit_Standard", "CNVkit_Refitted")

catch7 = plotDefinitionComparisons(mWithNorm, vFiles, vLevels, vLabels, OUT1)


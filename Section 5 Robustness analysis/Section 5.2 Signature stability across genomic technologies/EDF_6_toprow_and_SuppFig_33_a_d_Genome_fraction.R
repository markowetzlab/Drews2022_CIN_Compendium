## Genome fraction 

rm(list=ls(all=TRUE))

library(this.path)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(Cairo)
library(GenomicRanges)
library(lemon)

theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))

BASE=dirname(this.path())
GOLD=file.path(BASE, "input/TCGA_478_Samples_SNP6_GOLD.rds")
OUT=file.path(BASE, "output")
## Technically not needed as I round the data now.
WIGGLE=0.1

dtGold = readRDS(GOLD)
dtGold$segVal = round(dtGold$segVal)

## Function
calcGenomicFraction = function(dtGold, dtSamp, NAME, WIGGLE = 0.1) {
  
  allSamps = levels(dtGold$sample)
  lGenomFrac = lapply(allSamps, function(thisSample) {
    
    ## Catch if sample is not present in this condition
    if(nrow(dtSamp[ dtSamp$sample == thisSample, ]) == 0) {
      out = c(thisSample, NA)
      return(out)
    }
    
    ## Extract for each sample
    thisGold = makeGRangesFromDataFrame(dtGold[ dtGold$sample == thisSample, ], keep.extra.columns = TRUE)
    thisSamp = makeGRangesFromDataFrame(dtSamp[ dtSamp$sample == thisSample, ], keep.extra.columns = TRUE)
    
    ## Identify how much this sample has been covered by SNP6 (this is our 100%)
    genomeCov = sum(width(thisGold))
    
    ## Identify where the two overlap
    hits = findOverlaps(query = thisGold, subject = thisSamp)
    
    ## First, extend each granges object to the number of overlaps => gr[queryHits(hits)]
    ## Second, find the parallel (pairwise?) intersection to adjust the segment boundaries to the actual
    ## overlapping region
    overlaps = pintersect(thisGold[queryHits(hits)], thisSamp[subjectHits(hits)])
    
    ## Third, since we start with the Gold as query, we need to add the segval from the subject
    overlaps$newSegVal = thisSamp$segVal[subjectHits(hits)]
    
    ## Test whether the new segments is close to the gold standard ones
    overlaps$overlap = overlaps$segVal - WIGGLE < overlaps$newSegVal & 
      overlaps$segVal + WIGGLE > overlaps$newSegVal
    gFrac = signif(sum( width(overlaps[ overlaps$overlap ]) ) / genomeCov, 4)
    
    out = c(thisSample, gFrac)
    return(out)
    
  })
  
  dtGFrac = data.table(do.call(rbind, lGenomFrac))
  colnames(dtGFrac) = c("Sample", "GenomeFraction")
  dtGFrac$GenomeFraction = as.numeric(dtGFrac$GenomeFraction)
  dtGFrac$Name = NAME
  
  return(dtGFrac)
   
}

cohohrtGFrac = function(LISTSAMPLE, NAMES, PLOTNAMES) {
  
  # Read data
  lWGSNP6 = lapply(LISTSAMPLE, function(thisSamp) {
    dtSamp = readRDS(thisSamp) 
    dtSamp$segVal = round(dtSamp$segVal)
    return(dtSamp)
  })
  names(lWGSNP6) = NAMES
  
  # Calc genome fraction
  lGFrac1 = lapply(1:length(lWGSNP6), function(i) {
    
    print(i)
    dtSamp = lWGSNP6[[i]]
    thisName = names(lWGSNP6)[i]
    calcGenomicFraction(dtGold, dtSamp, NAME = thisName)
    
  })
  
  # Plot boxplot
  dtGFrac1 = rbindlist(lGFrac1)
  dtGFrac1$Name = factor(dtGFrac1$Name, levels = names(lWGSNP6),
                         labels = PLOTNAMES)
  
  pBox = ggplot(dtGFrac1, aes(x = Name, y = GenomeFraction)) + geom_boxplot()
  pBoxZoom = ggplot(dtGFrac1, aes(x = Name, y = GenomeFraction)) + geom_boxplot() +
    coord_cartesian(ylim = c(0.75, 1))
  pDens = ggplot(dtGFrac1, aes(x = GenomeFraction, group = Name, colour = Name)) + geom_density()
  
  ## Put in a list and return
  lOut = list(data = dtGFrac1, plotBox = pBox, plotBoxZoom = pBoxZoom, plotDensity = pDens)
  return(lOut)
  
}


## Testing 
# 
# dtSamp = readRDS(SAMPLE)
# dtSamp$segVal = round(dtSamp$segVal)
# 
# dtGFrac = calcGenomicFraction(dtGold, dtSamp, NAME = "Penalty70")
# 
# ggplot(dtGFrac, aes(x = GenomeFraction)) + geom_histogram()
# ggplot(dtGFrac, aes(x = Name, y = GenomeFraction)) + geom_boxplot()



#### SNP6 - no matched normal
LISTSAMPLE = c(file.path(BASE, "input/NewASCAT_CELnoNorm/NewASCAT_478TCGA_CELnoNorm_penalty35.rds"),
               file.path(BASE, "input/NewASCAT_CELnoNorm/NewASCAT_478TCGA_CELnoNorm_penalty50.rds"),
               file.path(BASE, "input/NewASCAT_CELnoNorm/NewASCAT_478TCGA_CELnoNorm_penalty70.rds"),
               file.path(BASE, "input/NewASCAT_CELnoNorm/NewASCAT_478TCGA_CELnoNorm_penalty100.rds"),
               file.path(BASE, "input/NewASCAT_CELnoNorm/NewASCAT_478TCGA_CELnoNorm_penalty140.rds"))
NAMES = c("Penalty35", "Penalty50", "Penalty70", "Penalty100", "Penalty140")
PLOTNAMES = c(35, 50, 70, 100, 140)

lSNP6NoNorm = cohohrtGFrac(LISTSAMPLE, NAMES, PLOTNAMES)
saveRDS(lSNP6NoNorm, file.path(OUT, "1_SNP6_NoNorm.rds"))



#### WGS - downsampled to SNP6 positions
LISTSAMPLE = c(file.path(BASE, "input/PCAWG_WGS_downSNP6/478_PCAWG_WGS_downSNP6_withNormals_penalty35.rds"),
               file.path(BASE, "input/PCAWG_WGS_downSNP6/478_PCAWG_WGS_downSNP6_withNormals_penalty50.rds"),
               file.path(BASE, "input/PCAWG_WGS_downSNP6/478_PCAWG_WGS_downSNP6_withNormals_penalty70.rds"),
               file.path(BASE, "input/PCAWG_WGS_downSNP6/478_PCAWG_WGS_downSNP6_withNormals_penalty100.rds"),
               file.path(BASE, "input/PCAWG_WGS_downSNP6/478_PCAWG_WGS_downSNP6_withNormals_penalty140.rds"))
NAMES = c("Penalty35", "Penalty50", "Penalty70", "Penalty100", "Penalty140")
PLOTNAMES = c(35, 50, 70, 100, 140)

lWGSNP6 = cohohrtGFrac(LISTSAMPLE, NAMES, PLOTNAMES)
saveRDS(lWGSNP6, file.path(OUT, "2_WGS_downSNP6.rds"))



#### WGS - downsamples to shallow WGS
LISTSAMPLE = c(file.path(BASE, "input/ASCATsc_WGS_downsWGS/ASCAT.sc_WGS_downsWGS_alpha0.001.rds"),
               file.path(BASE, "input/ASCATsc_WGS_downsWGS/ASCAT.sc_WGS_downsWGS_alpha0.01.rds"),
               file.path(BASE, "input/ASCATsc_WGS_downsWGS/ASCAT.sc_WGS_downsWGS_alpha0.05.rds"))
NAMES = c("Alpha0.001", "Alpha0.01", "Alpha0.05")
PLOTNAMES = c(0.001, 0.01, 0.05)

lWGSsWGS = cohohrtGFrac(LISTSAMPLE, NAMES, PLOTNAMES)
saveRDS(lWGSsWGS, file.path(OUT, "3_WGS_downsWGS.rds"))



#### WES - on-target reads
LISTSAMPLE = c(file.path(BASE, "input/ASCAT_WES_pilot_ontarget/ASCAT_TCGA_WES_478samples_penalty35.rds"),
               file.path(BASE, "input/ASCAT_WES_pilot_ontarget/ASCAT_TCGA_WES_478samples_penalty50.rds"),
               file.path(BASE, "input/ASCAT_WES_pilot_ontarget/ASCAT_TCGA_WES_478samples_penalty70.rds"),
               file.path(BASE, "input/ASCAT_WES_pilot_ontarget/ASCAT_TCGA_WES_478samples_penalty100.rds"),
               file.path(BASE, "input/ASCAT_WES_pilot_ontarget/ASCAT_TCGA_WES_478samples_penalty140.rds"))
NAMES = c("Penalty35", "Penalty50", "Penalty70", "Penalty100", "Penalty140")
PLOTNAMES = c(35, 50, 70, 100, 140)

lWESontarget = cohohrtGFrac(LISTSAMPLE, NAMES, PLOTNAMES)
saveRDS(lWESontarget, file.path(OUT, "4_WES_ontarget.rds"))



#### WES - off-target
LISTSAMPLE = c(file.path(BASE, "input/ASCATsc_WES_off_target/ASCAT.sc_WES_offtarget_30k_alpha0.001.rds"),
               file.path(BASE, "input/ASCATsc_WES_off_target/ASCAT.sc_WES_offtarget_30k_alpha0.01.rds"),
               file.path(BASE, "input/ASCATsc_WES_off_target/ASCAT.sc_WES_offtarget_30k_alpha0.05.rds"))
NAMES = c("Alpha0.001", "Alpha0.01", "Alpha0.05")
PLOTNAMES = c(0.001, 0.01, 0.05)

lWESofftarget = cohohrtGFrac(LISTSAMPLE, NAMES, PLOTNAMES)
saveRDS(lWESofftarget, file.path(OUT, "5_WES_offtarget.rds"))


#### WGS - Full data 
LISTSAMPLE = c(file.path(BASE, "input/PCAWG_WGS_full_data/PCAWG_ABSOLUTE_CN_profiles.rds"),
               file.path(BASE, "input/PCAWG_WGS_full_data/PCAWG_Battenberg_CN_profiles.rds"),
               file.path(BASE, "input/PCAWG_WGS_full_data/PCAWG_Sclust_CN_profiles.rds"),
               file.path(BASE, "input/PCAWG_WGS_full_data/PCAWG_Consensus_CN_profiles.rds"),
               file.path(BASE, "input/PCAWG_WGS_full_data/478_PCAWG_WGS_downSNP6_withNormals_penalty70.rds"))
NAMES = c("ABSOLUTE", "Battenberg", "Sclust", "Consensus", "ASCAT SNP6 pos.")
PLOTNAMES = c("ABSOLUTE", "Battenberg", "Sclust", "Consensus", "ASCAT SNP6 pos.")

lWGSABSOLUTE = cohohrtGFrac(LISTSAMPLE, NAMES, PLOTNAMES)
saveRDS(lWGSABSOLUTE, file.path(OUT, "5_WGS_Full_data.rds"))


#### WES - two callers
LISTSAMPLE = c(file.path(BASE, "input/ShallowWGS_Two_callers/ASCAT_TCGA_WES_478samples_penalty70.rds"),
               file.path(BASE, "input/ShallowWGS_Two_callers/TCGA_WES_Sequenca.rds"),
               file.path(BASE, "input/ShallowWGS_Two_callers/TCGA_WES_CNVkit_Standard.rds"),
               file.path(BASE, "input/ShallowWGS_Two_callers/TCGA_WES_CNVkit_Refitted.rds"))
NAMES = c("ASCAT70", "Sequenca", "CNVkit_standard", "CNVkit_refitted")
PLOTNAMES = c("ASCAT70", "Sequenca", "CNVkit_standard", "CNVkit_refitted")

lWESTwoCallers = cohohrtGFrac(LISTSAMPLE, NAMES, PLOTNAMES)
saveRDS(lWESTwoCallers, file.path(OUT, "5_WES_Two_callers.rds"))


## Heatmap bins
lAll = list(lSNP6NoNorm, lWGSNP6, lWGSsWGS, lWESontarget, lWESofftarget)
# lAll = list(lWGSABSOLUTE)
# lAll = list(lWESTwoCallers)
lHeatmapPlots = lapply(lAll, function(thisList) {
  
  dtDat = thisList$data
  dtDat$Cat = cut(dtDat$GenomeFraction, breaks = seq(0, 1, by = 0.1),
                  labels = paste0("<", seq(0.1, 1, by = 0.1)*100, "%"))
  
  dtSummary = data.table(melt(table(dtDat$Name, dtDat$Cat)))
  dtSummary$Var1 = factor(dtSummary$Var1, levels = levels(dtDat$Name))
  
  dtSummary$Count = cut(dtSummary$value,breaks=c(-1,10,50,100,200,300,max(dtSummary$value,na.rm=T)),
                  labels=c("0-10", "10-50","50-100","100-200","200-300",">300"))
  dtSummary$Count = factor(as.character(dtSummary$Count), levels = rev(levels(dtSummary$Count)))
  
  ggplot(dtSummary, aes(x = Var1, y = Var2, fill = Count)) + geom_tile(colour = "white") + 
    scale_fill_manual(values=c("#d53e4f","#f46d43","#fdae61","#fee08b","#e6f598","#abdda4"),
                      na.value = "grey90", drop=FALSE) +
    labs(x = "ASCAT Penalty", y = "Genome Agreement", fill = "Sample\nNumber")
  
})

heatmapRow = lHeatmapPlots[[1]] + lHeatmapPlots[[2]] + lHeatmapPlots[[3]] +
  lHeatmapPlots[[4]] + lHeatmapPlots[[5]] + plot_layout(nrow = 1, guides = "collect")

cairo_pdf(file.path(OUT, "6_EDF_6_Plot_Heatmap_bins.pdf"), width = 180/25.4, height = 50/25.4)
print(heatmapRow); dev.off()

ggsave(file.path(OUT, "6_EDF_6_Plot_Heatmap_bins.png"), heatmapRow, width = 180, height = 50, units = "mm")
ggsave(file.path(OUT, "6_EDF_6_Plot_Heatmap_bins.svg"), heatmapRow, width = 180, height = 50, units = "mm")


## Additional figures for reviews
lAll2 = list(lWGSABSOLUTE, lWESTwoCallers)
lHeatmapPlots2 = lapply(lAll, function(thisList) {
  
  dtDat = thisList$data
  dtDat$Cat = cut(dtDat$GenomeFraction, breaks = seq(0, 1, by = 0.1),
                  labels = paste0("<", seq(0.1, 1, by = 0.1)*100, "%"))
  
  dtSummary = data.table(melt(table(dtDat$Name, dtDat$Cat)))
  dtSummary$Var1 = factor(dtSummary$Var1, levels = levels(dtDat$Name))
  
  dtSummary$Count = cut(dtSummary$value,breaks=c(-1,10,50,100,200,300,max(dtSummary$value,na.rm=T)),
                        labels=c("0-10", "10-50","50-100","100-200","200-300",">300"))
  dtSummary$Count = factor(as.character(dtSummary$Count), levels = rev(levels(dtSummary$Count)))
  
  ggplot(dtSummary, aes(x = Var1, y = Var2, fill = Count)) + geom_tile(colour = "white") + 
    scale_fill_manual(values=c("#d53e4f","#f46d43","#fdae61","#fee08b","#e6f598","#abdda4"),
                      na.value = "grey90", drop=FALSE) +
    labs(x = "ASCAT Penalty", y = "Genome Agreement", fill = "Sample\nNumber")
  
})


ggsave(file.path(OUT, "6_SuppFig_33_a_Plot_Heatmap_bins_PCAWG_Full_data.svg"), lHeatmapPlots2[[1]], width = 80, height = 50, units = "mm")
ggsave(file.path(OUT, "6_SuppFig_33_a_Plot_Heatmap_bins_PCAWG_Full_data.png"), lHeatmapPlots2[[1]], width = 80, height = 50, units = "mm")

ggsave(file.path(OUT, "6_SuppFig_33_d_Plot_Heatmap_bins_TCGA_WES_Two_callers.svg"), lHeatmapPlots2[[2]], width = 80, height = 50, units = "mm")
ggsave(file.path(OUT, "6_SuppFig_33_d_Plot_Heatmap_bins_TCGA_WES_Two_callers.png"), lHeatmapPlots2[[2]], width = 80, height = 50, units = "mm")


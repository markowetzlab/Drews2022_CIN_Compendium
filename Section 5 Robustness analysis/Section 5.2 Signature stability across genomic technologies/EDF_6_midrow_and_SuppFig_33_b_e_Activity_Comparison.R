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
# From https://stackoverflow.com/questions/52297978/decrease-overal-legend-size-elements-and-text
addSmallLegend <- function(myPlot, pointSize = 0.5, textSize = 3, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}


plotActivityComparisons = function(dtWithNorm, vFiles, vNames, vLabels, OUTNAMEBASE) {
  
  # Load files
  lNoNorm = lapply(vFiles, function(thisFile) readRDS(thisFile))
  names(lNoNorm) = vNames
  
  ### Compare cosine similarities
  lCos = sapply(lNoNorm, function(thisNoNorm) {
    
    mWithMatched = mWithNorm[ match(rownames(thisNoNorm), rownames(mWithNorm)), ]
    vCos = sapply(1:nrow(mWithMatched), function(i) { cosine(mWithMatched[i,], thisNoNorm[i,]) })
    
    
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  dtCos = data.table(melt(lCos))
  # CEL files noNorm
  dtCos$L1 = factor(dtCos$L1,
                    levels = vNames,
                    labels = vLabels)
  
  # Not needed
  # p2a = ggplot(dtCos, aes(x = value, colour = L1)) + geom_density(alpha = 0.5, adjust = 0.5) + 
  #   xlab("Cosine similarity to matched normal") + ylab("Density") + labs(colour = "Penalty") + 
  #   ggtitle("Density adjust parameter = 0.5") + theme(legend.position = "none") + ylim(0, 15) +
  #   scale_colour_brewer(palette = "Dark2") + coord_capped_cart(left = "both", bottom = "both")
  # p2b = ggplot(dtCos, aes(x = value, colour = L1)) + geom_density(alpha = 0.5, adjust = 1) + 
  #   xlab("Cosine similarity to matched normal") + ylab("Density") + labs(colour = "Penalty") + 
  #   ggtitle("Density adjust parameter = 1 (default)") + theme(legend.position = "none") + ylim(0, 15) +
  #   scale_colour_brewer(palette = "Dark2") + coord_capped_cart(left = "both", bottom = "both")
  p2c = ggplot(dtCos, aes(x = value, colour = L1)) + geom_density(alpha = 0.5, adjust = 2) +
    xlab("Cosine similarity to matched normal") + ylab("Density") + labs(colour = "Penalty") +
    scale_y_continuous(breaks = seq(0,15, 5), limits = c(0, 15)) +
    scale_x_continuous(breaks = c(0.5, 0.75, 1), limits = c(0.5, 1)) +
    coord_capped_cart(left = "both", bottom = "both")
    
  # Make legend smaller and move
  p2c = addSmallLegend(p2c, pointSize = 6, textSize = 6, spaceLegend = 0.4) + theme(legend.position = c(0.35, 0.5))
  
  # Add correct colour scale
  if(length(levels(dtCos$L1)) == 5 & "35" %in% levels(dtCos$L1)) {
    vCols = c("35" = "#5377AE", "50" = "#6F9440", "70" = "#B3522B", "100" = "#D68E31", "140" = "#5C4E94")
  } else if(length(levels(dtCos$L1)) == 3 & "alpha 0.01" %in% levels(dtCos$L1)) {
    vCols = c("alpha 0.001" = "#A6CEE3", "alpha 0.01" = "#B2DF8A", "alpha 0.05" = "#FB9A99")
  } else if("ABSOLUTE" %in% levels(dtCos$L1)) {
    vCols = c("ABSOLUTE" = "#F8766D", "Battenberg" = "#7CAE00", "Sclust" = "#00BFC4", "Consensus" = "#C77CFF", "SNP6 Gold" = "#000000")
  } else if("Sequenca" %in% levels(dtCos$L1)) {
    vCols = c("Sequenca" = "#7CAE00", "CNVkit_standard" = "#00BFC4", "CNVkit_refitted" = "#C77CFF", "ASCAT70" = "#000000")
  } else {
    warning("No pre-set colour scale matched the data. Ignore, if you use new data.")
  }
  p2c = p2c + scale_colour_manual(values = vCols)
  
  # Output plots into multiple devices
  # p2 = p2a + p2b + p2c + plot_layout(widths = c(0.3, 0.3, 0.4))
  cairo_pdf(paste0(OUTNAMEBASE, ".pdf"), height = 5.6/2.54, width = 4.6/2.54); print(p2c); dev.off()
  ggsave(paste0(OUTNAMEBASE, ".svg"), p2c, height = 56, width = 46, units = "mm")
  
  return(p2c)
}


## Paths
BASE=dirname(this.path())

# Load gold standard (comparison signature activities)
mWithNorm = readRDS(file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures.rds"))
dtWithNorm = data.table(melt(mWithNorm))


# EDF 6 midrow: WGS down to shallow WGS
vFiles = c(file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_WGS_downsWGS_alpha0.01.rds"),
           file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_WGS_downsWGS_alpha0.05.rds"),
           file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_WGS_downsWGS_alpha0.001.rds"))
vNames = c("alpha0.01", "alpha0.05", "alpha0.001")
vLabels = c("alpha 0.01", "alpha 0.05", "alpha 0.001")
OUTNAMEBASE=file.path(BASE, "output/EDF_6_midrow_Activities_cel478_vs_WGS_downsWGS_threealphas_CosineSim")

catch1 = plotActivityComparisons(dtWithNorm, vFiles, vNames, vLabels, OUTNAMEBASE)



# EDF 6 midrow: WES on-target
vFiles = c(file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_WES_ontarget_penalty35.rds"),
           file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_WES_ontarget_penalty50.rds"),
           file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_WES_ontarget_penalty70.rds"),
           file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_WES_ontarget_penalty100.rds"),
           file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_WES_ontarget_penalty140.rds"))
vNames=c("penalty35", "penalty50", "penalty70", 
  "penalty100", "penalty140")
vLabels=c("35", "50", "70", "100", "140")
OUTNAMEBASE=file.path(BASE, "output/EDF_6_midrow_Activities_cel478_vs_WESontarget_CosineSim")

catch2 = plotActivityComparisons(dtWithNorm, vFiles, vNames, vLabels, OUTNAMEBASE)



# EDF 6 midrow: WES off-target
vFiles = c(file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_WES_offtarget_30K_alpha0.001.rds"), 
           file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_WES_offtarget_30K_alpha0.01.rds"),
           file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_WES_offtarget_30K_alpha0.05.rds"))
vNames=c("alpha0.001", "alpha0.01", "alpha0.05")
vLabels=c("alpha 0.001", "alpha 0.01", "alpha 0.05")
OUTNAMEBASE=file.path(BASE, "output/EDF_6_midrow_Activities_cel478_vs_WES_offtarget_30K_threealphas_CosineSim")

catch3 = plotActivityComparisons(dtWithNorm, vFiles, vNames, vLabels, OUTNAMEBASE)




# EDF 6 midrow: SNP6 arrays without matched normal
vFiles = c(file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_NewASCAT_CELnoNorm_penalty35.rds"), 
           file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_NewASCAT_CELnoNorm_penalty50.rds"),
           file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_NewASCAT_CELnoNorm_penalty70.rds"),
           file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_NewASCAT_CELnoNorm_penalty100.rds"),
           file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_NewASCAT_CELnoNorm_penalty140.rds"))
vNames=c("penalty35", "penalty50", "penalty70", "penalty100", "penalty140")
vLabels=c("35", "50", "70", "100", "140")
OUTNAMEBASE=file.path(BASE, "output/EDF_6_midrow_Activities_cel478_vs_cel478NoNorm_CosineSim")

catch4 = plotActivityComparisons(dtWithNorm, vFiles, vNames, vLabels, OUTNAMEBASE)



# EDF 6 midrow: WGS downsampled to SNP6
vFiles = c(file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_WGSdownSNP6_withNormals_penalty35.rds"), 
           file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_WGSdownSNP6_withNormals_penalty50.rds"),
           file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_WGSdownSNP6_withNormals_penalty70.rds"),
           file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_WGSdownSNP6_withNormals_penalty100.rds"),
           file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_WGSdownSNP6_withNormals_penalty140.rds"))
vNames=c("penalty35", "penalty50", "penalty70", "penalty100", "penalty140")
vLabels=c("35", "50", "70", "100", "140")
OUTNAMEBASE=file.path(BASE, "output/EDF_6_midrow_Activities_cel478_vs_WGS478downSNP6_withNorm_CosineSim")

catch5 = plotActivityComparisons(dtWithNorm, vFiles, vNames, vLabels, OUTNAMEBASE)



## SuppFig 33 b) WGS all data included
vFiles = c(file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_PCAWG_ABSOLUTE_CN_profiles.rds"), 
           file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_PCAWG_Battenberg_CN_profiles.rds"),
           file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_PCAWG_Sclust_CN_profiles.rds"),
           file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_PCAWG_Consensus_CN_profiles.rds"),
           file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_WGSdownSNP6_withNormals_penalty70.rds"))
vNames=c("ABSOLUTE", "Battenberg", "Sclust", "Consensus", "SNP6")
vLabels=c("ABSOLUTE", "Battenberg", "Sclust", "Consensus", "SNP6 Gold")
OUTNAMEBASE=file.path(BASE, "output/SuppFig_33_b_midrow_Activities_cel478_vs_WGS478downSNP6_withNorm_CosineSim")

catch6 = plotActivityComparisons(dtWithNorm, vFiles, vNames, vLabels, OUTNAMEBASE)



## SuppFig 33 e) WES (two callers) and all data included
vFiles = c(file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_WES_ontarget_penalty70.rds"), 
           file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_TCGA_WES_Sequenca.rds"),
           file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_TCGA_WES_CNVkit_Standard.rds"),
           file.path(BASE, "input/Activities/4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_TCGA_WES_CNVkit_Refitted.rds"))
vNames=c("ASCAT70", "Sequenca", "CNVkit_standard", "CNVkit_refitted")
vLabels=c("ASCAT70", "Sequenca", "CNVkit_standard", "CNVkit_refitted")
OUTNAMEBASE=file.path(BASE, "output/SuppFig_33_e_Activities_cel478_vs_WES_two_callers_CosineSim")

catch7 = plotActivityComparisons(dtWithNorm, vFiles, vNames, vLabels, OUTNAMEBASE)



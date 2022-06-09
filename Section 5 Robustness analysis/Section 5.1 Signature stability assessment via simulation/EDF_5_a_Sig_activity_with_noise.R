#### Plot signature-specific thresholds

rm(list=ls(all=TRUE))

library(this.path)
library(data.table)
library(reshape2)
library(mclust)
library(ggplot2)
library(ggthemes)
library(lemon)
library(RColorBrewer)
library(ggrastr)

theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6))

BASE=dirname(this.path())
ACTIVITIES=file.path(BASE, "input/Signature_Compendium_v5_Cosine-0.74_Activities_NAMESAPRIL21.rds")
## Currently too large for github - please download from Markowetz lab and CRUK servers:
## https://content.cruk.cam.ac.uk/fmlab/drews2022/2_Activities_fullTCGA_1000sims_10pGaussian_10pSamplePoisson_NAMESAPRIL21_RUNNOV20.rds
DAT = file.path(BASE, "input/2_Activities_fullTCGA_1000sims_10pGaussian_10pSamplePoisson_NAMESAPRIL21_RUNNOV20.rds")
OUTPUT=file.path(BASE, "output")

dtOri = data.table(melt(readRDS(ACTIVITIES)))
lSignatures = readRDS(DAT)


# Bring list of signature activities into one data table
lMelt = list()
allSigs = length(lSignatures)
for(i in 1:allSigs) {
  dtSignatures = data.table(melt(lSignatures[[i]] ))
  dtSignatures$iter = i
  lMelt[[length(lMelt)+1]] = dtSignatures
}

dtSigs = rbindlist(lMelt)

# Plot boxplot for each sample and signature
allSigs = levels(dtSigs$Var2)
lPlots = lapply(allSigs, function(thisSig) {
  
  dtCS2 = dtSigs[ dtSigs$Var2 == thisSig, ]
  
  # prepare original values
  dtOriCS2 = dtOri[ dtOri$Var2 == thisSig, ]
  
  # Sort samples by original signature value or by median simulation value
  newOrder = as.character(dtOriCS2$Var1)[ order(dtOriCS2$value, decreasing = TRUE) ]
  
  # Reorder factors
  dtOriCS2$Var1 = factor(as.character(dtOriCS2$Var1), levels = newOrder)
  dtCS2$Var1 = factor(as.character(dtCS2$Var1), levels = newOrder)
  
  # Plot IQR
  p1 = ggplot(dtCS2, aes(x = Var1, y = value)) + rasterise(geom_boxplot(outlier.shape = NA, coef = 0, lwd = 0.1), dpi = 900) +
    rasterise(geom_point(data = dtOriCS2, aes(x = Var1, y = value), colour = "red", size = 0.1), dpi = 900) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank()) + labs(y = paste("Activity of", thisSig)) + 
    scale_x_discrete(breaks = NULL)
  
  cairo_pdf(file.path(OUTPUT, paste0("Simulation_", thisSig, "_test.pdf")), width = 50/25.4, height = 45/25.4)
  print(p1); dev.off()
  ggsave(file.path(OUTPUT, paste0("Simulation_", thisSig, ".svg")), p1, width = 50/25.4, height = 45/25.4)
  
  return(p1)
})





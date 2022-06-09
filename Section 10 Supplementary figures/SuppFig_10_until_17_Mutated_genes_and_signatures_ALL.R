#### Plot mutated genes vs activity

rm(list=ls(all=TRUE))

library(this.path)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(lemon)
library(ggrastr)

theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))

BASE=dirname(this.path())
OUTFIGURES=file.path(BASE, "output")
OUTTABLES=file.path(BASE, "output")
MUTDATA=file.path(BASE, "input/Mutations_and_GoI_for_Supplementary_Figure.rds")
RESGENES=file.path(BASE, "input/Table_tTest_Genes_of_Interest_THRESH95_NAMESAPRIL2021.rds")

## Load and prepare activities
dtMut = readRDS(MUTDATA)
dtRes = readRDS(RESGENES)

## Annotate genes
dtHigh = dtRes[ dtRes$MeanDiff > 0.4 & dtRes$pAdjTTest < 0.005, ]
allSigs = unique(dtHigh$Signature)

dtMut$Significant = FALSE
for(thisSig in allSigs) {
  mutGenes = dtHigh$Driver[ dtHigh$Signature == thisSig ]
  dtMut$Significant[ dtMut$Signature == thisSig & dtMut$Gene %in% mutGenes ] = TRUE
}
dtMut$Significant = factor(dtMut$Significant)

## Example plot
ggplot(dtMut[dtMut$Gene=="JAK2", ], aes(x = Driver, y = Exposure, fill = Significant)) + 
  geom_boxplot(outlier.colour = NA, outlier.shape = NA) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.2) +
  facet_wrap(Signature ~ Gene, scales = "free_y", nrow = 1) + scale_fill_discrete(drop = FALSE)

GENESPERPAGE=6
allGenes = unique(dtMut$Gene)
lPlots = lapply(0:floor(length(allGenes)/GENESPERPAGE), function(thisPage) {
  
  theseGenes = allGenes[((thisPage*GENESPERPAGE)+1):((thisPage*GENESPERPAGE)+GENESPERPAGE)]
  theseMut = dtMut[ dtMut$Gene %in% theseGenes,]
  
  ggplot(theseMut, aes(x = Driver, y = Exposure, colour = Significant, fill = Significant)) + 
    rasterise(geom_jitter(width = 0.1, height = 0, alpha = 0.2, size = 0.1), dpi = 300) +
    geom_boxplot(outlier.colour = NA, outlier.shape = NA, alpha = 0.6) +
    labs(y = "Activity", x = "Mutation") + scale_x_discrete(guide = guide_axis(angle = 90)) +
    facet_wrap(Gene ~ Signature, scales = "free_y", nrow = GENESPERPAGE, ncol = 17) + scale_fill_discrete(drop = FALSE)
  
  
})

catchAll = lapply(1:length(lPlots), function(thisPage) {
  
  # 85 for last page (40mm per row)
  png(file.path(OUTFIGURES, paste0("SuppFig_10_to_17_Genes_and_Muts_Page_", thisPage, ".png")), res = 400, width = 200, height = 245, units = "mm")
  print(lPlots[[thisPage]]); dev.off()  
  
})


## Heatmap sig activities vs cancer types

rm(list=ls(all=TRUE))

library(this.path)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(lemon)
library(RColorBrewer)
library(gridExtra)

theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))

BASE=dirname(this.path())
OUTPLOTS=file.path(BASE, "output")
EXP=file.path(BASE, "input/Signature_Compendium_v5_Cosine-0.74_Activities_THRESH95_NAMESAPRIL21.rds")
META=file.path(BASE, "input/Metadata_TCGA_ASCAT_penalty70.rds")

ACTIVITYTHRESHOLD = 0.5

## Load data
exp=readRDS(EXP)
meta=readRDS(META)

## Prepare metadata file
metaStudy = meta[ ! is.na(meta$dCIN), ]
dtMeta = data.table(table(metaStudy$cancer_type))

## Prepare activities
dtExp = data.table(melt(exp))
colnames(dtExp) = c("Sample", "Signature", "Exposure")
dtExp$Cancer = meta$cancer_type[ match(dtExp$Sample, meta$name) ]

## Summarise active samples per cancer type
dtActive = dtExp[ dtExp$Exposure > 0,] 
dtSummary = data.table(table(dtActive$Signature, dtActive$Cancer))


### Plot Number 1: Binary Heatmap
dtSummary$Prop = dtSummary$N / dtMeta$N[ match(dtSummary$V2, dtMeta$V1) ]
dtSummary$Active = dtSummary$Prop > ACTIVITYTHRESHOLD

## Prepare plot
dtSummary$V1 = factor(dtSummary$V1, levels = rev(colnames(exp)))

pOut = ggplot(dtSummary, aes(x = V2, y = V1, fill = Active)) + geom_tile(width = 0.9, height = 0.9) +
  scale_x_discrete(position = "top", guide = guide_axis(angle = -90)) +
  scale_fill_manual(values = c("white", "grey30")) + 
  theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.line.x = element_line(size = 0.5), axis.ticks = element_line(size = 0.5), 
        axis.ticks.length.x = unit(.1, "cm"), plot.margin = unit(c(0, 0, 0, 0), "null")) +
  labs(x = "TCGA cancer types") + coord_capped_cart(top = "both") 

## 116 x 67.7 mm => current space in figure 2
## Not needed anymore
# cairo_pdf(file.path(OUTPLOTS, "Figure_2_Aetiologies_heatmap_binary.pdf"), width = 67.7/25.4, height = 116/25.4)
# print(pOut); dev.off()
# ggsave(file.path(OUTPLOTS, "Figure_2_Aetiologies_heatmap_binary.svg"), pOut, width = 67.7/25.4, height = 116/25.4)



### Plot Number 2: Four activity categories
BREAKS = c(0, 0.05, 0.25, 0.5, 0.75, 1)
LABELS = c("0-5%", "5-25%", "25-50%", "50-75%", ">75%" )
dtSummary$Cut = cut(dtSummary$Prop, breaks = BREAKS, labels = LABELS, include.lowest = TRUE)

# # Cluster cancers based on discretisation
# dtSummary$CutNum = as.numeric(dtSummary$Cut)
# mTest = dcast(V2 ~ V1, value.var = "CutNum", data = dtSummary)
# rownames(mTest) = mTest$V2
# mTest$V2 = NULL
# mTest = as.matrix(mTest)
# hTest = hclust(dist(mTest, method = "manhattan"), method = "ward.D2")
# 
# dtSummary$V2 = factor(as.character(dtSummary$V2), 
#                          levels = unique(dtSummary$V2)[hTest$order])


pOut2 = ggplot(dtSummary, aes(x = V2, y = V1, fill = Cut)) + geom_tile(width = 0.9, height = 0.9) +
  scale_x_discrete(position = "top", guide = guide_axis(angle = -90)) +
  # scale_fill_manual(values = c("white", "lightskyblue1", "deepskyblue1", "deepskyblue3", "deepskyblue4")) +
  scale_fill_manual(values = c("white", brewer.pal(4, "Blues"))) +
  # scale_fill_brewer(palette = "Blues") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.line.x = element_line(size = 0.5), axis.ticks = element_line(size = 0.5), 
        axis.ticks.length.x = unit(.1, "cm"), plot.margin = unit(c(0, 0, 0, 0), "null"),) +
  labs(x = "TCGA cancer types") + coord_capped_cart(top = "both")

## Put legend in its own figure to not mess with figure dimensions
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
getLegend = function(a.gplot) {
  tmp = ggplot_gtable(ggplot_build(a.gplot))
  leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend = tmp$grobs[[leg]]
  
  return(legend)
}

myLegend = getLegend(pOut2)
pLegend = grid.arrange(arrangeGrob(myLegend))

# Now remove legend
pOut2 = pOut2 + theme(legend.position = "none")

cairo_pdf(file.path(OUTPLOTS, "Figure_2_Aetiologies_heatmap.pdf"), width = 67.7/25.4, height = 117.5/25.4)
print(pOut2); dev.off()
ggsave(file.path(OUTPLOTS, "Figure_2_Aetiologies_heatmap.svg"), pOut2, width = 67.7/25.4, height = 117.5/25.4)

ggsave(file.path(OUTPLOTS, "Figure_2_Aetiologies_heatmap_legend.svg"), pLegend, width = 2, height = 1)


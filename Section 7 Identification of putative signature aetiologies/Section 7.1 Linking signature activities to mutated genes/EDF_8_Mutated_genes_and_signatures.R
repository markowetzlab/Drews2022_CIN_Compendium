# Standard packages used throughout repo

rm(list=ls(all=TRUE))

library(this.path)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
# For dealing with strings
library(stringr)
# For building panel C
library(patchwork)
library(scales)

theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))

## Alpha of rectangle pattern
ALPHA=0.1

## Paths
BASE=dirname(this.path())
OUT=file.path(BASE, "output")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## Panel C: Genes of Interest (GoI)
TGOI=file.path(BASE, "input/Table_tTest_Genes_of_Interest_THRESH95_NAMESAPRIL2021.rds")
GENEDETAILS=file.path(BASE, "input/Reactome_and_Driver_genes_and_their_IDs.txt")
REACTOMEDETAILS=file.path(BASE, "input/Reactome_genes_pathway_dictionary.txt")

# Even though they might not be mutated, we show them to the reader.
MOSTIMPORTANT=c("KRAS", "BRAF", "NOTCH1", "NF1", "SMAD4", "ARID1A", "JAK2", "POLE", 
                "MYCN", "EGFR", "ATM", "BRCA1_rescued") # "BRCA1_noRescue"
REMOVEBRCA1NOTRESCUED = TRUE


## Load GoIs and determine genes to be plotted
dtGenes = fread(GENEDETAILS)
## Convert labels in lookup data frame
dtGenes$Status[ grepl("oncogene", dtGenes$Status) ] = "Oncogene"
dtGenes$Status[ grepl("tsg", dtGenes$Status) ] = "TSG"
dtGenes$Status[ dtGenes$Status == "" ] = "Either"

dtReact = fread(REACTOMEDETAILS)
dtGoI = readRDS(TGOI)

## Prepare dtGoI
## Simple filtering
# dtSigGoI = dtGoI[ dtGoI$MeanDiff > 0 & dtGoI$pAdjTTest < PVALTHRESH, ]

## High-quality filtering - only positive, only highly significant results
dtSigGoI = dtGoI[ dtGoI$MeanDiff > 0.4 & dtGoI$pAdjTTest < 0.005, ]
allGenes = unique(c(dtSigGoI$Driver, MOSTIMPORTANT))

## Remove BRCA1 if wanted
if(REMOVEBRCA1NOTRESCUED) {
  allGenes = allGenes[ allGenes != "BRCA1_noRescue" ]
  dtSigGoI = dtSigGoI[ dtSigGoI$Driver != "BRCA1_noRescue", ]
}

## Save sorted output
dtSorted = dtSigGoI[ order(dtSigGoI$Signature, -dtSigGoI$MeanDiff), ]
write.table(dtSorted, file.path(OUT, "Table_Positive_significant_genes_per_signature.txt"), 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Add Reactome source
dtSubset = dtReact[ dtReact$New %in% allGenes, c("New", "ReactomePathway")]
dtAllGenes = rbind(dtSubset, data.table("New" = allGenes[ ! allGenes %in% dtSubset$New ], "ReactomePathway" = NA))
colnames(dtAllGenes) = c("GeneName", "Pathway")

# Add driver status - Source: Bailey et al. 2018
dtAllGenes$Status = dtGenes$Status[ match(dtAllGenes$GeneName, dtGenes$name) ]
dtAllGenes$Status[ dtAllGenes$GeneName == 'MYCN' ] = "Oncogene"
if(! REMOVEBRCA1NOTRESCUED) {
  dtAllGenes$Status[ dtAllGenes$GeneName == 'BRCA1_noRescue' ] = "TSG"  
}
dtAllGenes$Status[ dtAllGenes$GeneName == 'BRCA1_rescued' ] = "TSG"

# Rename BRCA1
if(! REMOVEBRCA1NOTRESCUED) {
  dtAllGenes$GeneName[ dtAllGenes$GeneName == "BRCA1_noRescue" ] = "BRCA1 not rescued"
}
dtAllGenes$GeneName[ dtAllGenes$GeneName == "BRCA1_rescued" ] = "BRCA1 rescued"

# Rename Reactome pathways
dtAllGenes$Pathway[ dtAllGenes$Pathway == "DNA-Replication" & 
                      ! is.na(dtAllGenes$Pathway) ] = "DNA Replication"
dtAllGenes$Pathway[ dtAllGenes$Pathway == "DNA-Repair" & 
                      ! is.na(dtAllGenes$Pathway) ] = "DNA Repair"
dtAllGenes$Pathway[ dtAllGenes$Pathway == "Chromatin-Organisation" & 
                      ! is.na(dtAllGenes$Pathway) ] = "Chromatin Organisation"
dtAllGenes$Pathway[ dtAllGenes$Pathway == "Cell-Cycle" & 
                      ! is.na(dtAllGenes$Pathway) ] = "Cell Cycle"

# Rename TSG to tumour suppre...
dtAllGenes$Status[ dtAllGenes$Status == "TSG" ] = "Tumour suppressor"

## Plot
# Set order for categories (gene status and Reactome pathways) and then genes
mAllGenes = melt(dtAllGenes, id.vars = "GeneName")
mAllGenes$value = factor(mAllGenes$value, 
                         levels = c("Tumour suppressor", "Oncogene", "Either", 
                                    "Unknown", "DNA Replication", 
                                    "DNA Repair", "Chromatin Organisation", "Cell Cycle"))

datGeneOrder = rev(unique(as.character(mAllGenes$GeneName[ order(as.character(mAllGenes$value), 
                                                             as.character(mAllGenes$GeneName), 
                                                             decreasing = TRUE) ])))
mAllGenes$GeneName = factor(mAllGenes$GeneName, levels = datGeneOrder)

# Set colours
colsStatus = c("Oncogene" = "#e41a1c", "Tumour suppressor" = "#4daf4a", "Unknown" = "#377eb8",
               "Either" = "#ff7f00")

# Only driver status in pA1
numGenes = length(unique(mAllGenes$GeneName[ mAllGenes$variable == "Status" ]))
numCats = length(unique(mAllGenes$value[ mAllGenes$variable == "Status" ]))


pC1 = ggplot(mAllGenes[ mAllGenes$variable == "Status", ], 
             aes(y = value, x = GeneName, fill = value)) + 
  geom_blank() +
  annotate("rect", 
           xmin = (seq(1, numGenes, by = 2) - 0.5), 
           xmax = (seq(1, numGenes, by = 2) + 0.5), 
           ymin = 0.5, ymax = numCats + 0.5, alpha = ALPHA) + 
  annotate("rect", xmin = 0,xmax = numGenes + 0.5, 
           ymin = (seq(1, numCats, by = 2) -0.5), 
           ymax = (seq(1, numCats, by = 2) +0.5), alpha = ALPHA) + 
  geom_tile() +
  theme(legend.position = "none", axis.title.y = element_blank(), plot.margin = margin( t = 5, r = 5, b = 0, l = 5, "pt")) +
  scale_x_discrete(position = "top", drop = FALSE, guide = guide_axis(angle = 90)) + scale_fill_manual(values = colsStatus) +
  coord_equal() + xlab("Significant genes")


# # Rotate plot
# dtOrderRot = mAllGenes[mAllGenes$variable == "Status",]
# as.character(dtOrderRot$GeneName)[ order(dtO) ]
# datGeneOrderRot = unique(as.character(mAllGenes$GeneName)[ order(ordered(mAllGenes$value), 
#                                                              as.character(mAllGenes$GeneName), 
#                                                              decreasing = TRUE) ])
# mAllGenes$GeneName2 = factor(as.character(mAllGenes$GeneName), levels = datGeneOrderRot)
# pC1Rot = ggplot(mAllGenes[ mAllGenes$variable == "Status", ], 
#              aes(x = value, y = GeneName2, fill = value)) + 
#   geom_blank() + 
#   annotate("rect", 
#            ymin = (seq(1, numGenes, by = 2) - 0.5), 
#            ymax = (seq(1, numGenes, by = 2) + 0.5), 
#            xmin = 0.5, xmax = numCats + 0.5, alpha = ALPHA) + 
#   annotate("rect", ymin = 0, ymax = numGenes + 0.5, 
#            xmin = (seq(1, numCats, by = 2) -0.5), 
#            xmax = (seq(1, numCats, by = 2) +0.5), alpha = ALPHA) + 
#   geom_tile() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none",
#         axis.title.y = element_blank(), plot.margin = margin( t = 5, r = 5, b = 0, l = 5, "pt")) +
#   scale_y_discrete(drop = FALSE) + scale_fill_manual(values = colsStatus) +
#   coord_equal() + xlab("Significant genes")
# 


## Only Reactome origin in pA2
# Give BRCA1 special cases also the same Reactome pathways
addon1 = mAllGenes[ mAllGenes$GeneName == "BRCA1", ]
addon1$GeneName = "BRCA1 rescued"
if(! REMOVEBRCA1NOTRESCUED) {
  addon2 = mAllGenes[ mAllGenes$GeneName == "BRCA1", ]
  addon2$GeneName = "BRCA1 not rescued"
}

mAllGenes = mAllGenes[ ! mAllGenes$GeneName %in% c("BRCA1 rescued", "BRCA1 not rescued"), ]
if(! REMOVEBRCA1NOTRESCUED) {
  mAllGenes = rbind(mAllGenes, addon1, addon2)
} else {
  mAllGenes = rbind(mAllGenes, addon1)
}

# Remove NA and only keept Reactome pathway entries
mPlotC2 = mAllGenes[ mAllGenes$variable == "Pathway" & ! is.na( mAllGenes$value ), ]
colsReactome = c("Cell Cycle" = "#8dd3c7", "Chromatin Organisation" = "#ffffb3",
                 "DNA Repair" = "#bebada", "DNA Replication" = "#fb8072")
numCats = length(unique(mPlotC2$value))
pC2 = ggplot(mPlotC2, aes(y = value, x = GeneName, fill = value)) +
  geom_blank() +
  annotate("rect",
           xmin = (seq(1, numGenes, by = 2) - 0.5),
           xmax = (seq(1, numGenes, by = 2) + 0.5),
           ymin = 0.5, ymax = numCats + 0.5, alpha = ALPHA) +
  annotate("rect", xmin = 0,xmax = numGenes + 0.5,
           ymin = (seq(1, numCats, by = 2) - 0.5),
           ymax = (seq(1, numCats, by = 2) + 0.5), alpha = ALPHA) +
  geom_tile() + coord_equal() +
    theme(axis.title.y = element_blank(), legend.position = "none", axis.ticks.x = element_blank(),
          axis.title.x = element_blank(), axis.text.x = element_blank(),
          plot.margin = margin( t = 5, r = 5, b = 0, l = 5, "pt")) +
  scale_fill_manual(values = colsReactome) + scale_x_discrete(drop = FALSE)


## Signature vs genes of interest
# Change BRCA1 names in dtSigGoI (comes directly from hardrive) as well
dtSigGoI$Driver[ dtSigGoI$Driver == "BRCA1_noRescue" ] = "BRCA1 not rescued"
dtSigGoI$Driver[ dtSigGoI$Driver == "BRCA1_rescued" ] = "BRCA1 rescued"
dtSigGoI$Status = dtAllGenes$Status[ match(dtSigGoI$Driver, dtAllGenes$GeneName) ]
dtSigGoI$Driver = factor(dtSigGoI$Driver, levels = datGeneOrder)
dtSigGoI$Signature = factor(dtSigGoI$Signature, levels = rev(paste0("CX", 1:17)))

## Version with mean difference in colour and point shape indicating number of mutated samples
dtSigGoI$SampleCat = cut(dtSigGoI$SampsMut, breaks = c(0, 10, 100, 500, 6000), labels = c("<10", "<100", "<500", ">500"), include.lowest = TRUE)


# Some mutated genes are associated with mean shifts of larger than 1. Put in upper limit.
dtPlotC3 = dtSigGoI
dtPlotC3$MeanDiff[ dtPlotC3$MeanDiff > 1 ] = 1.1
numCats = length(levels(dtPlotC3$Signature))
pC3 = ggplot(dtPlotC3, aes(y = Signature, x = Driver, colour = MeanDiff, shape = SampleCat)) + 
  geom_blank() +
  annotate("rect", 
           xmin = (seq(1, numGenes, by = 2) - 0.5), 
           xmax = (seq(1, numGenes, by = 2) + 0.5), 
           ymin = 0.5, ymax = numCats + 0.5, alpha = ALPHA) + 
  annotate("rect", xmin = 0,xmax = numGenes + 0.5, 
           ymin = (seq(1, numCats, by = 2) - 0.5), 
           ymax = (seq(1, numCats, by = 2) + 0.5), alpha = ALPHA) +
  geom_point(size = 2.4) + scale_shape_manual( values = c(16,17,15,7) ) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
        plot.margin = margin( t = 5, r = 5, b = 5, l = 5, "pt")) +
  scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE) +
  scale_colour_viridis_c(rescaler = function(x, to = c(0, 1), from = NULL) {
    ifelse(x<1, scales::rescale(x, to = to,from = c(min(x, na.rm = TRUE), 1)), 1)},
    breaks = seq(0,1, by = 0.25), labels = c(seq(0,0.75,by = 0.25), ">1")) +
   guides(fill = guide_colourbar(barwidth = 0.5, barheight = 15)) + coord_equal() + 
  labs(colour = "Mean\ndifference", shape = "Number of\nsamples", y = "CIN signatures")

saveRDS(dtPlotC3, file.path(OUT, "Supp_Table_Mutated_genes.rds"))

pC = pC1 / pC2 / pC3

## Save plots
cairo_pdf(file.path(OUT, "EDF_8_a_Significant_genes.pdf"), width = 190/25.4, height = 120/25.4)
print(pC); dev.off()

ggsave(file.path(OUT, "EDF_8_a_Significant_genes.svg"), pC, width = 8.25, height = 5)


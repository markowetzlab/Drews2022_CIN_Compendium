## Compare activities of simulated samples with their produced signature

rm(list=ls(all=TRUE))

## Libraries
library(this.path)
library(ggplot2)
library(reshape2)
library(stringr)
library(Polychrome)
library(patchwork)
library(ggthemes)
library(lemon)

theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))

#### Paths and functions
thisPath=dirname(this.path())
IN=file.path(thisPath, "input")
OUT=file.path(thisPath, "output")

# This file is created by running "Simulation_Calc_activities.R" script from the CINGenomeSimulation repo.
mAct=readRDS(file.path(IN, "Out_6_Activities_20each_N240.rds"))
dfOverview=readRDS(file.path(IN, "Out_2_Step_3_Overview_simulation_20each_N240.rds"))
SCALE=TRUE

if(SCALE) {
  mAct = scale(mAct, center = TRUE, scale = TRUE)
}

## Names that start with a number got an additional X as prefix (part of R naming matrices)
## Repair names
dfAct = reshape2::melt(mAct)
colnames(dfAct) = c("Samples", "Signature", "Activity")
dfAct$Samples = as.character(dfAct$Samples)

vNamesRepair = sapply(dfAct$Samples, function(thisName) {
  if(nchar(thisName) == 7) {
    return(substr(thisName, 2, 7))
  } else {
    return(thisName)
  }
})

dfAct$Samples = vNamesRepair

# Match order for plotting
thisOrder = dfOverview$Samples[ order(dfOverview$Signature) ]

dfAct$Samples = factor(vNamesRepair, levels = rev(thisOrder))
dfAct$Simulation = dfOverview$Signature[ match(as.character(dfAct$Samples), dfOverview$Samples) ]

## Plot heatmap
p1 = ggplot(dfAct, aes(x = 1, y = Samples, fill = Simulation)) + geom_tile() + labs(x = "Simulation") +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) + scale_fill_manual(values = unname(glasbey.colors(18)))
p2 = ggplot(dfAct, aes(x = Signature, y = Samples, fill = Activity)) + geom_tile() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), axis.line.y = element_blank()) + scale_fill_gradient2(low = "blue", mid = "white", midpoint = 0, high = "red")
pOut1 = p1 + p2 + plot_layout(nrow = 1, widths = c(0.05, 0.95), guides = "collect")

ggsave(file.path(OUT, "SuppFig_5_Activities_20each_N240_with_Annotation_Heatmap.png"), 
       pOut1, width = 16, height = 16, units = "cm", bg = "white")
ggsave(file.path(OUT, "SuppFig_5_Activities_20each_N240_with_Annotation_Heatmap.svg"), 
       pOut1, width = 16, height = 16, units = "cm", bg = "white")



# Simplify names: remove for WGD the last entry - not dominant signature anyways
lGroup = lapply(unique(dfAct$Simulation), function(thisGroup) {
  
  # print(thisGroup)
  dfGroup = dfAct[ dfAct$Simulation == thisGroup, ]
  
  # Trim down three element names because only two are dominant
  splitGroup = str_split(thisGroup, "_")[[1]]
  newGroupName = ifelse(length(splitGroup) == 3, paste0(splitGroup[1:2], collapse = "_"), thisGroup)
  
  # Replace Chr with CHR
  newGroupName = gsub("Chr", "CHR", newGroupName)
  
  dfGroup$ShortSim = newGroupName
  return(dfGroup)
})

dfShortGroup = do.call(rbind, lGroup)
dfShortGroup$Simulation = dfShortGroup$ShortSim
dfShortGroup$ShortSim = NULL

##### Labeling top XX dominant samples for each each signature
topXX = c("HRDLST" = 100, "ecDNA" = 100, "CHR" = 100, "WGDearly" = 96, "WGDlate" = 24)

# Label top XX for each signature
allSigs = levels(dfShortGroup$Signature)
lSig = lapply(allSigs, function(thisSig) {

  dfSig = dfShortGroup[ dfShortGroup$Signature == thisSig, ]
  dfSig$TopXX = "No"

  # Actually label samples
  numTop = topXX[names(topXX) == thisSig]
  dfSig$TopXX[ order(dfSig$Activity, decreasing = TRUE)[1:numTop] ] = "Yes"

  return(dfSig)
})

dfSig = do.call(rbind, lSig)


## Show of the top 100 samples how many have been predicted right.

## Correct HRDLST to LST
dfSig$Signature = factor(as.character(dfSig$Signature), levels = allSigs,
                         labels = c("LST", "ecDNA", "CHR", "WGDearly", "WGDlate"))

## Loop over each signature and get all samples with "Yes" or "No" (top XX samples)
allSigs = levels(dfSig$Signature)
lTab = lapply(allSigs, function(thisSig) {
  
  dfThisSig = dfSig[ dfSig$Signature == thisSig, ]
  dfThisSig$RightGroup = ifelse(grepl(thisSig, dfThisSig$Simulation), "RightGroup", "WrongGroup")
  
  dfThisSig$RightGroup = factor(dfThisSig$RightGroup, levels = c("RightGroup", "WrongGroup"))
  dfThisSig$TopXX = factor(dfThisSig$TopXX, levels = c("No", "Yes"))
  
  # Call melt explicitly as reshape2 is now deprecated.
  dfSummary = reshape2::melt(table(dfThisSig$TopXX, dfThisSig$RightGroup))
  dfSummary$Signature = thisSig 
  
  return(dfSummary)
  
})

dfTab = do.call(rbind, lTab)
colnames(dfTab) = c("DominantSample", "Simulation", "N", "Signature")

dfRight = dfTab[dfTab$Simulation=="RightGroup",]
p3 = ggplot(dfRight, aes(x = Signature, y = N, fill = DominantSample)) + 
  geom_bar(stat = "identity", position = "fill") + 
  labs(y = "Proportion of correctly identified samples", fill = "Identified\nsamples with\ndominant signature?") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("Yes" = "#1b9e77", "No" = "#d95f02" ))
p4 = ggplot(dfRight, aes(x = Signature, y = N, fill = DominantSample)) + 
  geom_bar(stat = "identity") + 
  labs(y = "Sample number", fill = "Identified\nsamples with\ndominant signature?") +
  scale_fill_manual(values = c("Yes" = "#1b9e77", "No" = "#d95f02" ))

## Calculate average accuracy
# Number of times a dominant signature was modelled
allDoms = sum(dfRight$N)

# Number of times we got it right
allRights = sum(dfRight$N[ dfRight$DominantSample == "Yes" ])
allRights/allDoms

pOut3 = p3 + p4 + plot_layout(guides = "collect")

ggsave(file.path(OUT, "SuppFig_6_Activities_20each_N240_with_Annotation_Barplot.png"), 
       pOut3, width = 16, height = 6, units = "cm", bg = "white")
ggsave(file.path(OUT, "SuppFig_6_Activities_20each_N240_with_Annotation_Barplot.svg"), 
       pOut3, width = 16, height = 6, units = "cm", bg = "white")


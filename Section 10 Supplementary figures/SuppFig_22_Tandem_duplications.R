## Plot tandem dup subgroups

rm(list=ls(all=TRUE))

library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(lemon)
library(patchwork)
library(scales)
library(viridisLite)
library(ggrepel)
library(this.path)

theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))

## Functions
# Source: https://stackoverflow.com/questions/52297978/decrease-overal-legend-size-elements-and-text
addSmallLegend = function(myPlot, pointSize = 0.5, textSize = 6, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}



## Boilerplate
BASE=dirname(this.path())
OUT=file.path(BASE, "output")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)


## Signature files
ACT=file.path(BASE, "input/Signature_Compendium_v5_Cosine-0.74_Activities_THRESH95_NAMESAPRIL21.rds")
TANDEMDUPS=file.path(BASE, "input/Menghi2018_S3_TDP_inclTCGA.csv")
ZSCORES=TRUE


## Load data
tandups = fread(TANDEMDUPS)
tandups = tandups[grepl("TCGA", tandups$Sample_ID),]
tandups$Name = substr(tandups$Sample_ID, 1, 12)

act=readRDS(ACT)
if(ZSCORES) {
  dtAct = melt(scale(act, center = TRUE, scale = TRUE))  
} else {
  dtAct = melt(act)  
}

colnames(dtAct) = c("Sample", "Signature", "Exposure")


## Transfer number of TDPs, TDP score, TDP status and group
dtAct$NumTDP = tandups$No_TD[ match(substr(dtAct$Sample,1,12), tandups$Name) ]
dtAct$ScoreTDP = tandups$TDP_score[ match(substr(dtAct$Sample,1,12), tandups$Name) ]
dtAct$StatusTDP = tandups$TDP_status[ match(substr(dtAct$Sample,1,12), tandups$Name) ]
dtAct$GroupTDP = tandups$TDP_group[ match(substr(dtAct$Sample,1,12), tandups$Name) ]

dtTDP = dtAct[ ! is.na(dtAct$NumTDP), ]
dtTDP$GroupTDP = factor(dtTDP$GroupTDP, levels = unique(dtTDP$GroupTDP))



## Correlate with sigs
lTDP = split(dtTDP, dtTDP$Signature)
lResults = lapply(lTDP, function(thisSig) {
  
  # Test for correlation with exposure and number of tdps
  corNums = cor(thisSig$Exposure, thisSig$NumTDP, method = "spearman")
  corNumsP = cor.test(thisSig$Exposure, thisSig$NumTDP, method = "spearman")$p.value
  
  # ... TDP score
  corScore = cor(thisSig$Exposure, thisSig$ScoreTDP, method = "spearman")
  corScoreP = cor.test(thisSig$Exposure, thisSig$ScoreTDP, method = "spearman")$p.value
  
  # ... status of TDP
  dfStatus = aggregate(Exposure ~ StatusTDP, data = thisSig, mean)
  meanDiff = dfStatus$Exposure[ dfStatus$StatusTDP == 1 ] - 
    dfStatus$Exposure[ dfStatus$StatusTDP == 0 ]
  meanDiffP = pairwise.t.test(thisSig$Exposure, thisSig$StatusTDP, pool.sd = FALSE,
                              p.adjust.method = "none")$p.value
  
  # ... group of TDP
  dfGroup = aggregate(Exposure ~ GroupTDP, data = thisSig, mean)
  vGroups = dfGroup$Exposure
  names(vGroups) = dfGroup$GroupTDP
  
  # When working with z-scores, sometimes groups have no variance but values unequal to 0.
  # This breaks the pairwise.t.test function. Therefore adding tiny amount of jitter to
  # the groups that have zero variance so the function doesn't fail. And then later assign
  # a p-value of 1, so result is basically discarded.
  if(ZSCORES) {
    dfVar = aggregate(Exposure ~ GroupTDP, data = thisSig, var)
    toBeFixed = as.character(dfVar$GroupTDP[ dfVar$Exposure == 0 ])
    valsToBeFixed = thisSig$Exposure[ thisSig$GroupTDP %in% toBeFixed ]
    valsToBeFixed = valsToBeFixed + rnorm(length(valsToBeFixed), mean = 0, sd = 1e-05)
    thisSig$Exposure[ thisSig$GroupTDP %in% toBeFixed ] = valsToBeFixed
  }
  pairGroup = pairwise.t.test(thisSig$Exposure, thisSig$GroupTDP, pool.sd = FALSE,
                              p.adjust.method = "none")
  ## We only compare non-TDP to other groups and not all groups with each other
  pairPVals = pairGroup$p.value[,"Non TDP"]
  pairPVals = p.adjust(pairPVals, method = "BH")
  
  ## Set p-value to 1 for the groups without variance (= few samples with identical values).
  if(ZSCORES) {
    pairPVals[ names(pairPVals) %in% toBeFixed ] = 1
  }
  
  # Collect all numbers, round them and output them
  vNums = c(corNums, corNumsP, corScore, corScoreP, meanDiff, meanDiffP, vGroups, pairPVals)
  vNums = signif(vNums, 4)
  
  return(vNums)
  
})

mResults = do.call(rbind, lResults)
colnames(mResults) = c("CorNumTDP", "PvalCorNumTDP", "CorTDPScore", "PvalCorTDPScore",
                        "MeanDiff", "PvalMeandiff", "MeanNonTDP", "MeanTDP1", "MeanTDP12", 
                        "MeanTDP2", "MeanTDP13", "MeanTDP23", "MeanTDP3", "PvalMeanTDP1",
                        "PvalMeanTDP12", "PvalMeanTDP2", "PvalMeanTDP13", "PvalMeanTDP23",
                        "PvalMeanTDP3")
dfResults = data.frame(mResults)
dfResults$PadjCorNumTDP = p.adjust(as.numeric(dfResults$PvalCorNumTDP), method = "BH")
dfResults$PadjCorTDPScore = p.adjust(as.numeric(dfResults$PvalCorTDPScore), method = "BH")

# dfResults$Signature = rownames(dfResults)

# write.table(x = dfResults, file = file.path(OUT, "Table_TandemDup_Results.txt"), 
            # row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)


## Plot the means and the p-values separately
# Prepare number of category
dtCat = dtTDP[ dtTDP$Signature == "CX1", ]
table(dtCat$GroupTDP)
## Now a bit un-computational but using the number of samples in a manual way down below.
CLASSLABELS = c("Non-TDP (525 samples)", "Class 1 (36 samples)", 
                "Class 1/2mix (12 samples)", "Class 2 (42 samples)",
                "Class 1/3mix (5 samples)", "Class 2/3mix (6 samples)",
                "Class 3 (5 samples)")

#### Plot 1: Correlation with number of tandem duplications
dfPlot1 = dfResults
dfPlot1$Plot = -log(dfPlot1$PadjCorNumTDP)
dfPlot1$Sig = factor(dfPlot1$Plot > 3, levels = c("TRUE", "FALSE"))
p1 = ggplot(dfPlot1, aes(x = CorNumTDP, y = Plot, colour = Sig)) + geom_point() +
  geom_text_repel(label = rownames(dfPlot1), force = 25) + 
  scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "grey80")) +
  labs(colour = "Significant\nSpearman corr.", x = "Spearman correlation", y = "-log(q-value)",
       title = "Signature activity vs. number of tandem duplications") + 
  coord_capped_cart(bottom = "both", left = "both") + theme(legend.position = c(0.55, 0.85))
  

#### Plot 2: Correlation with tandem duplication score
dfPlot2 = dfResults
dfPlot2$Plot = -log(dfPlot2$PadjCorTDPScore)
dfPlot2$Sig = factor(dfPlot2$Plot > 3, levels = c("TRUE", "FALSE"))
p2 = ggplot(dfPlot2, aes(x = CorTDPScore, y = Plot, colour = Sig)) + geom_point() +
  geom_text_repel(label = rownames(dfPlot1), force = 25) + 
  scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "grey80")) +
  labs(colour = "Significant\nSpearman corr.", x = "Spearman correlation", y = "-log(q-value)",
       title = "Signature activity vs. tandem duplication score") + 
  coord_capped_cart(bottom = "both", left = "both") + theme(legend.position = c(0.55, 0.85))


#### Plot 3: Per tandem group t-test
## Split data into two data.frames. One for heatmap and other for geom_point (p-value)
# Prepare TDP classes
colsToPlot = c("MeanNonTDP", "MeanTDP1", "MeanTDP12", "MeanTDP2", "MeanTDP13", "MeanTDP23", "MeanTDP3")
dfTDP = dfResults[, colsToPlot]
dfTDP$Signature = rownames(dfTDP)
mTDP = melt(dfTDP, id.vars = "Signature")

mTDP$Signature = factor(mTDP$Signature, levels = paste0("CX", 1:17))
mTDP$variable = factor(mTDP$variable, levels = colsToPlot, labels = CLASSLABELS)

# Prepare p-values
colsToPlot = c("PvalMeanTDP1", "PvalMeanTDP12", "PvalMeanTDP2", "PvalMeanTDP13", 
               "PvalMeanTDP23", "PvalMeanTDP3")
dfPVals = dfResults[, colsToPlot]
dfPVals$Signature = rownames(dfPVals)
mPVals = melt(dfPVals, id.vars = "Signature")
mPVals$Signature = factor(mPVals$Signature, levels = paste0("CX", 1:17))
mPVals$variable = factor(mPVals$variable, levels = colsToPlot, labels = CLASSLABELS[-1])

# Values are already q-values (see lapply above)
mPVals$QValue = cut(mPVals$value, breaks = c(0, 0.01, 0.05, 0.15, 1),
                    labels = c("<0.01", "<0.05", "<0.15", ">0.15"))

# Remove non-significant results
mPVals = mPVals[ mPVals$QValue != ">0.15", ]

## Plot
p3 = ggplot(mTDP, aes(x = Signature, y = variable, fill = value)) + geom_tile() + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  xlab("CIN signature") + ylab("Tandem Duplicator Class") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + coord_equal() +
  geom_point(data = mPVals, aes(x = Signature, y = variable, shape = QValue), size = 2) +
  labs(fill = "Difference in\nmean signature\nactivity", shape = "q-value")

p3 = addSmallLegend(p3, pointSize = 2, textSize = 7, spaceLegend = 0.5)

ggsave(file.path(OUT, "SuppFig_22_Tandem_Dups_Panel_A.svg"), p1, width = 90/25.4, height =  90/25.4)
ggsave(file.path(OUT, "SuppFig_22_Tandem_Dups_Panel_B.svg"), p2, width = 90/25.4, height =  90/25.4)
ggsave(file.path(OUT, "SuppFig_22_Tandem_Dups_Panel_C.svg"), p3, width = 180/25.4, height =  90/25.4)


#### Save output for web portal
lOut = list("NumTD" = data.table("Signature" = rownames(dfResults), 
                                          "Covariate" = "Number of tandem duplications", 
                                          "Correlation" = dfResults$CorNumTDP, 
                                          "pAdj" = dfResults$PadjCorNumTDP),
            "ScoreTD" = data.table("Signature" = rownames(dfResults), 
                                 "Covariate" = "Tandem duplication score", 
                                 "Correlation" = dfResults$CorTDPScore, 
                                 "pAdj" = dfResults$PadjCorTDPScore),
            "TDClass" = data.table("Signature" = mTDP$Signature, 
                                 "Covariate" = "Tandem duplication class", 
                                 "Class" = mTDP$variable,
                                 "Difference" = mTDP$value),
            "TDClassPadj" = data.table("Signature" = mPVals$Signature, 
                                       "Covariate" = "Tandem duplication class", 
                                       "Class" = mPVals$variable,
                                       "pAdj" = mPVals$value,
                                       "pAdjCat" = mPVals$QValue))
# saveRDS(lOut, file.path(OUTTABLES, "Covariates_TandemDups.rds"))





## Test BRCA1 and BRCA2 mutants for signatures CX2 and CX5 with correction for CX3

rm(list=ls(all=TRUE))

library(this.path)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(lemon)
library(RColorBrewer)
library(car)

theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))

BASE=dirname(this.path())
OUTFIGURES=file.path(BASE, "output")

EXP=file.path(BASE, "input/Signature_Compendium_v5_Cosine-0.74_Activities_THRESH95_NAMESAPRIL21.rds")
SURV=file.path(BASE, "input/Survival_and_BRCA_Status_fullTCGA.rds")

## Load and prepare activities
exp = readRDS(EXP)
surv = readRDS(SURV)

## CX3/CX2 are already in surv. Add CX5
surv$CX5 = exp[match(surv$Name, rownames(exp)),"CX5"]


## Scale activities
surv$sCX3 = scale(surv$CX3)
surv$sCX2 = scale(surv$CX2)
surv$sCX5 = scale(surv$CX5)


## Remove double WT category
surv$Status[ surv$Status == "WT BRCA1/2 TCGA" ] = "WT BRCA1/2"


#### All cancers
## Control group is: WT BRCA1/2
## Test groups: germline BRCA1+LOH, germline BRCA2+LOH, BRCA1 Hypermethyl., 
## RAD51C Hypermethyl., somatic BRCA1+LOH, somatic BRCA2+LOH

## Prepare data (filter samples for correct BRCA1/2 status, simplify for lm)
testAll = surv[ surv$Status %in% c("WT BRCA1/2", "germline BRCA1+LOH", "germline BRCA2+LOH",
                                   "BRCA1 Hypermethyl.", "RAD51C Hypermethyl.", "somatic BRCA1+LOH", "somatic BRCA2+LOH"), ]

#### Testing for difference between groups with correcting for the other two IHR sigs
## CX2
lmCX2 = lm(sCX2 ~ 1 + Status + Cancer +sCX3 + sCX5 + Status*sCX3 + Status*sCX5, testAll)
summary(lmCX2)

## CX5
lmCX5 = lm(sCX5 ~ 1 + Status + Cancer + sCX3 + sCX2 + Status*sCX3 + Status*sCX2, testAll)
summary(lmCX5)

## CX3
lmCX3 = lm(sCX3 ~ 1 + Status + Cancer + sCX2 + sCX5 + Status*sCX2 + Status*sCX5, testAll)
summary(lmCX3)

## ==> Single categories (Status) always at least two significant for the three sigs

## Plots showing multivariate regression
avPlots(lmCX2)
avPlots(lmCX5)
avPlots(lmCX3)

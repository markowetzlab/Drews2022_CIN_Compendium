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

theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))

## Windows
BASE=dirname(this.path())
OUTFIGURES=file.path(BASE, "output")

EXP=file.path(BASE, "input/Signature_Compendium_v5_Cosine-0.74_Activities_THRESH95_NAMESAPRIL21.rds")
SURV=file.path(BASE, "input/Survival_and_BRCA_Status_fullTCGA.rds")
CLASSFILE=file.path(BASE, "input/CX3CX2_Clinical_classifier.rds")



## Little function for applying the classifier to a data frame with columns "CX3", "CX2" and "Name"
applyClinClass = function(dtSurvOV, lModel) {
  
  ## Scale the activities
  mTCGAOV = as.matrix(dtSurvOV[ , c("CX3", "CX2")])
  rownames(mTCGAOV) = dtSurvOV$Name
  smTCGAOV = sweep(sweep(mTCGAOV, 2, lModel$mean, FUN = '-'), 2, lModel$scale, FUN = "/")
  
  dtSurvOV$sCX3 = smTCGAOV[,"CX3"]
  dtSurvOV$sCX2 = smTCGAOV[,"CX2"]
  
  ## Apply classifier
  if(identical(dtSurvOV$Name, rownames(smTCGAOV))) {
    dtSurvOV$Classifier = ifelse(dtSurvOV$sCX3 >= dtSurvOV$sCX2, "Predicted sensitive", "Predicted resistant")  
  } else {
    stop("Something went wrong.")
  }
  
  return(dtSurvOV)
  
}




## Load and prepare activities
exp = readRDS(EXP)
surv = readRDS(SURV)
lModel = readRDS(CLASSFILE)


## Step 1: Scale according to classifier
## Step 2: Classify
surv = applyClinClass(surv, lModel)


## General stats
table(surv$Classifier)

prop.table(table(surv$Classifier))

prop.table(table(surv$Classifier, surv$Cancer), 2)

table(surv$Classifier, surv$Cancer)

## This script correlates age at diagnosis and number of CNAs attributable to a signature

rm(list=ls(all=TRUE))

library(this.path)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(lemon)
# Robust linear regression
library(MASS)
# F (Wald) test
library(sfsmisc)

## Paths
BASE=dirname(this.path())
OUT=file.path(BASE, "output")

# Cancer type
META=file.path(BASE, "input/Metadata_TCGA_ASCAT_penalty70.rds")
# Age
CLIN=file.path(BASE, "input/Survival_and_BRCA_Status_fullTCGA.rds")
MUTS=file.path(BASE, "input/Mutations_SNVs_AMPs_DELs_in_TCGA.rds")
# Estimated number of CNAs per signature
CNAS=file.path(BASE, "input/Estimated_CNAs_per_signature_and_sample.rds")
# Prior knowledge about cancer development times
GERSTUNG=file.path(BASE, "input/Gerstung2020_Figure5B_Age_PCAWG_Cancers.txt")

## Load data
# exp = readRDS(EXP)
# raw = readRDS(RAW)
meta = readRDS(META)
clin = readRDS(CLIN)
muts=readRDS(MUTS)
cnas = readRDS(CNAS)
gerstung=fread(GERSTUNG)


#### Do it like Alexandrov et al. did in the recent PCAWG Mutational Signature paper 
### Step 1: Data cleaning
## "Before evaluating the association between age and the activity of a mutational signature, 
## all outliers for both age and numbers of mutations attributed to a signature in a cancer type 
## were removed from the data. An outlier was defined as any value outside three standard 
## deviations from the mean value."

### Step 2: Robust linear regression for each cancer type and signature
## "A robust linear regression model that estimated the slope of the line and whether this 
## slope was significantly different from zero (F test; P value < 0.05) was performed using 
## the MATLAB function robustfit (https://www.mathworks.com/help/stats/robustfit.html) with 
## default parameters. The P values from the F tests were corrected using the Benjamini–Hochberg 
## procedure for false discovery rates. Results are available at syn12030687 and syn20317940."

## Prepare
# meta => Cancer type
# clin => Age 
# cnas => Number of attributed CNAs
meta$Age = as.numeric(clin$Age[ match(meta$name, clin$Name) ])
metaFilt = meta[ meta$dCIN,]
# Remove samples with no age (NA) => 116 samples
metaFilt = metaFilt[ ! is.na(metaFilt$Age), ]
# Transfer estimated number of CNAs for CX1
metaFilt$AttrCNAs = cnas[ match(metaFilt$name, rownames(cnas)),"CX1"]


### Filter cancer types
## Because CNAs need time and lots of samples, we take only those from Gerstung et al. (Nature, 2020),
## that took longer than 10 years to develop for 1x (see Figure 5B). 10 years for accumulation of 
## CNAs and 1x because that has the most direct link to age (see Figure 5A).
## Also we take cancers with more than 200 samples (random large number) to ensure we are powered
## enough. => 9 cancer types survive => 2 of those are significant

## Additional information on Alexandrov et al. method:
## 1) Age and SBS1 is surprisingly complicated. Alexandrov et al. (Nature, 2020) used two approaches to
## calculate signatures and exposures for PCAWG: Signature Analyzer (SA) and SigProfiler (SP). Both
## approaches yield different significant results with cancer types:
## SP: Breast, Medullo, PiloAstro, Lymph-BNHL, Ovary, Thy-Adenoca => 6
## SA: Breast, Medullo, PiloAsto, Lymph-BNHL, Ovary, Thy-Adenoca (THCA), GBM, Oligo, Colorect, 
## Cervix (CESC), Kidney, Panc-Adenoca, Panc-Endocrine, Prost, Skin => 15
## Names PCAWG: https://www.nature.com/articles/s41586-020-1969-6/tables/1
## 
## 2) SBS1 and CX1 is also quite diverse:
## Only cancer-specific significant correlation is BRCA. Other, non-significant but mostly
## around the BRCA correlation strength are STAD (stomach), READ (rectal), HNSC (head), SARC, BLCA.
## Some like SKCM (skin) or kidneys have strong negative correlations but not significant.
## OV has medium strong positive and LUSC almost zero correlation with SBS1.
## Only overlap to SP is BRCA.
## 
## 3) Age is extremely tricky. Only OV and LUSC are significant. UCEC is towards significance.
## BRCA is strongly negative but not significant. THCA, CESC, Kidneys, COAD, ESCA, PRAD have all
## negative slopes.
## UVM, THCA, CESC and LIHC (all the ones with strong negative slopes) are known to have distinct 
## and recurrent chromosomal imbalances.

cancers10years = unlist(strsplit(gerstung$TCGA[ ! is.na(gerstung$TCGA) ], ","))
tCancers = table(metaFilt$cancer_type)

cancersOfInterest = tCancers[ names(tCancers) %in% cancers10years & tCancers > 200 ]
metaFiltCoI = metaFilt[ metaFilt$cancer_type %in% names(cancersOfInterest), ]

## Just check whether age correlates to number of CNAs with correcting for cancer type
summary(lm(Age ~ AttrCNAs + cancer_type, data = metaFiltCoI))

rlmAge = rlm(Age ~ AttrCNAs + cancer_type, data = metaFiltCoI)
# Robust F-test (Wald test) to test for significance
ftestAge = f.robftest(rlmAge, var = "AttrCNAs")
# => positive but not significant



## Loop over cancer type and then over each signature
## On a cancer level: Remove age outliers
## On a signature level: 
## 1) Estimate number of CNAs attributed by the signature, remove outliers
## 2) Build robust linear regression

allCancers = unique(metaFiltCoI$cancer_type)
allSigs = colnames(cnas)
SDTHRESH = 3
lResult = lapply(allCancers, function(thisCancer) {
  
  # print(thisCancer)
  ## Age outliere removal
  metaCancer = metaFiltCoI[ metaFiltCoI$cancer_type == thisCancer, ]
  upperBound = mean(metaCancer$Age) + SDTHRESH*sd(metaCancer$Age)
  lowerBound = mean(metaCancer$Age) - SDTHRESH*sd(metaCancer$Age)
  metaCancer = metaCancer[ metaCancer$Age > lowerBound & metaCancer$Age < upperBound, ]
  
  ## Loop over signatures - not needed. Currently only interested in CX1
  lRLM = lapply("CX1", function(thisSig) {
  # lRLM = lapply(allSigs, function(thisSig) {
    
    ## CNA outlier removal  
    upperBound = mean(metaCancer$AttrCNAs) + SDTHRESH*sd(metaCancer$AttrCNAs)
    lowerBound = mean(metaCancer$AttrCNAs) - SDTHRESH*sd(metaCancer$AttrCNAs)
    metaCancer = metaCancer[ metaCancer$AttrCNAs > lowerBound & metaCancer$AttrCNAs < upperBound, ]
    
    ## Escape mechanism in some circumstances where there is no signature exposure
    # if(sum(metaCancer$RelExp) == 0) {
    #   out = c(thisCancer, thisSig, rep(0, 6), 1)
    #   return(out)
    # }
    
    ## Build robust linear regression
    # Using Huber weights as standard
    rlmAge = rlm(Age ~ AttrCNAs, data = metaCancer)
    # Robust F-test (Wald test) to test for significance
    ftestAge = f.robftest(rlmAge, var = "AttrCNAs")
    out = c(thisCancer, thisSig, coefficients(summary(rlmAge))[,1],
            coefficients(summary(rlmAge))[,3], ftestAge$statistic, ftestAge$df[2],
            ftestAge$p.value)
    
    return(out)
  })
  
  dtCancer = data.table(do.call(rbind, lRLM))
  colnames(dtCancer) = c("Cancer", "Signature", "CoeffIntercept", "CoeffAttrCNAs",
                         "tIntercept", "tExposure", "F", "df", "pVal")
  dtCancer$CoeffIntercept = as.numeric(dtCancer$CoeffIntercept)
  dtCancer$CoeffAttrCNAs = as.numeric(dtCancer$CoeffAttrCNAs)
  dtCancer$tIntercept = as.numeric(dtCancer$tIntercept)
  dtCancer$tExposure = as.numeric(dtCancer$tExposure)
  dtCancer$F = as.numeric(dtCancer$F)
  dtCancer$df = as.numeric(dtCancer$df)
  dtCancer$pVal = as.numeric(dtCancer$pVal)
  dtCancer$pAdj = p.adjust(dtCancer$pVal, method = "BH")

  return(dtCancer)
  
})

dtAge = rbindlist(lResult)
dtAge$pAdj = p.adjust(dtAge$pVal, method = "BH")
dtRes = dtAge[ dtAge$pAdj < 0.15, ]
dtCX1 = dtAge[ dtAge$Signature == "CX1", ]

## Save table
write.table(dtAge, file.path(OUT, "SuppInfo_Age.txt"),
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
## Save table
write.table(dtAge[ dtAge$pAdj < 0.05, ], file.path(OUT, "SuppInfo_Age_Significant.txt"),
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

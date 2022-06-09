#### Test genes for cancer-specific mutation behaviour

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
OUTTABLE=file.path(BASE, "output")
dir.create(OUTTABLE, showWarnings = FALSE, recursive = TRUE)

RAW=file.path(BASE, "input/Signature_Compendium_v5_Cosine-0.74_Activities_raw_THRESH95_NAMESAPRIL21.rds")
META=file.path(BASE, "input/Metadata_TCGA_ASCAT_penalty70.rds")
GENENAMES=file.path(BASE, "input/Table_tTest_Test_GoI_Enrichment_SigFilt_THRESH95_NAMESAPRIL2021.rds")
MUTS=file.path(BASE, "input/Mutations_SNVs_AMPs_DELs_in_TCGA.rds")

# Run special gene expression analysis for VHL?
SPECIAL=FALSE


## Load
mRaw = readRDS(RAW)
meta = readRDS(META)
dtSig = readRDS(GENENAMES)
dtMuts = readRDS(MUTS)


## Combine data
dtScale = data.table(melt(scale(mRaw)))
colnames(dtScale) = c("Sample", "Sig", "Activity")
dtScale$Cancer = factor(meta$cancer_type[ match(dtScale$Sample, meta$name)  ])

## Build multivariate model
allSigsWithMuts = unique(dtSig$Signature)
lAllSigs = lapply(allSigsWithMuts, function(thisSig) {
  
  print(thisSig)
  
  dtThis = dtScale[ dtScale$Sig == thisSig, ]
  
  allDriver = dtSig$Driver[ dtSig$Signature == thisSig ]
  # Remove BRCA1_rescue and BRCA1_noRescue for now
  allDriver = allDriver[ ! grepl("_", allDriver) ]
  
  
  lAllDriver = lapply(allDriver, function(thisDriver) {
    
    ## Prepare data
    # print(thisDriver)
    theseMutants = unique(as.character(dtMuts$sample[ dtMuts$gene == thisDriver ]))
    dtThis$Status = ifelse(dtThis$Sample %in% theseMutants, "Mutant", "WT")
    dtThis$Status = factor(dtThis$Status, levels = c("WT", "Mutant"))
    
    ## Fit the model and output relevant results
    lmModel = lm(Activity ~ 1 + Status + Cancer, dtThis)
    out = c(thisSig, thisDriver, summary(lmModel)$coefficients["StatusMutant",])
    
    return(out)
    
  } )
  
  ## Convert results
  dtAllDriver = data.table(do.call(rbind, lAllDriver))
  colnames(dtAllDriver) = c("Sig", "Driver", "Estimate", "StdError", "t", "pVal")
  dtAllDriver$Estimate = as.numeric(dtAllDriver$Estimate)
  dtAllDriver$StdError = as.numeric(dtAllDriver$Estimate)
  dtAllDriver$t = as.numeric(dtAllDriver$t)
  dtAllDriver$pVal = as.numeric(dtAllDriver$pVal)
  
  return(dtAllDriver)
  
})

## Merge results and correct p-values
dtResults = rbindlist(lAllSigs)
dtResults$pAdj = p.adjust(dtResults$pVal, method = "BH")
dtResults$Significant = dtResults$pAdj < 0.05

write.table(dtResults, file.path(OUTTABLE, "SuppMat_Significant_Mutant_Genes.txt"), quote = FALSE, sep = '\t', 
            row.names = FALSE, col.names = TRUE)


#### Special analysis for VHL - needs dtThis with "CX1" and "VHL"
if(SPECIAL) {
  
  # Please download file from CRUK CI server to the input folder. Address here:
  # --To be added--
  dtGeneExpr = readRDS("input/CosmicCompleteGeneExpression_Signature_Samples_only.rds")
  
  # Filter for VHL
  dtVHL = dtGeneExpr[ dtGeneExpr$gene == "VHL", ]
  
  # Match expression to dtThis 
  dtThis$VHL = dtGeneExpr$`z-score`[ match(substr(dtThis$Sample,1,12), dtGeneExpr$sample) ]
  dtFilt = dtThis[ ! is.na(dtThis$VHL), ]
  dtFilt$Cancer = factor(as.character(dtFilt$Cancer))
  
  # Just check that enough samples in other cancer types are left
  table(dtFilt$Cancer)
  
  # Correlate
  
  # Linear model correcting for tumour type
  lmVHL = lm(formula = Activity ~ 1 + VHL + Cancer, data = dtFilt)
  # p=0.531270
  summary(lmVHL)
  
  rlmVHL = rlm(formula = Activity ~ 1 + VHL + Cancer, data = dtFilt)
  # p-value = 0.3678
  fTestVHL = f.robftest(rlmVHL, var = "VHL")
  
}







rm(list=ls(all=TRUE))

library(this.path)
library(data.table)
library(reshape2)

## Paths
BASE=dirname(this.path())
OUTTABLE=file.path(BASE, "output")
dir.create(OUTTABLE, showWarnings = FALSE, recursive = TRUE)

RAW=file.path(BASE, "input/Signature_Compendium_v5_Cosine-0.74_Activities_raw_THRESH95_NAMESAPRIL21.rds")

COLS=file.path(BASE, "input/TCGA_colour_scheme.txt")
META=file.path(BASE, "input/Metadata_TCGA_ASCAT_penalty70.rds")


## Driver
CNSTATE=file.path(BASE, "input/CN_State_GoIs_TCGA.rds")
SNVMUTS=file.path(BASE, "input/SNV_Status_GoIs_TCGA.rds")
GENELOOKUP=file.path(BASE, "input/Reactome_and_Driver_genes_and_their_IDs.txt")

# Use to introduce additional knowledge for germline BRCA1 and BRCA2
BRCASAMPLES=file.path(BASE, "input/TCGA_Exposures_and_BRCA_Status_plusGene.rds")

# For 2nd round of reviews do additional analysis on CCNE1
REVIEW=FALSE




## Load files
raw = readRDS(RAW)
meta = readRDS(META)

# Driver
cns = readRDS(CNSTATE)
dels = cns[ cns$consequence == "DEL", ]
amps = cns[ cns$consequence == "AMP", ]
snv = readRDS(SNVMUTS)


# For providing additional data on BRCA1/2, RAD51C mutated samples
brcaMuts = readRDS(BRCASAMPLES)
brca1 = unique(as.character(brcaMuts$Sample[ grepl("[germline|somatic] BRCA1", 
                                                   brcaMuts$Status) ]))
brca2 = unique(as.character(brcaMuts$Sample[ grepl("[germline|somatic] BRCA2", 
                                                   brcaMuts$Status) ]))
rad51c = unique(as.character(brcaMuts$Sample[ grepl("RAD51C", brcaMuts$Status) ]))


## For checking oncogenes and tsgs with amps and dels
dtGeneLookup = fread(GENELOOKUP)





## Convert labels in lookup data frame
dtGeneLookup$Status[ grepl("oncogene", dtGeneLookup$Status) ] = "Oncogene"
dtGeneLookup$Status[ grepl("tsg", dtGeneLookup$Status) ] = "TSG"
dtGeneLookup$Status[ dtGeneLookup$Status == "" ] = "Ambivalent"

## Merge mutations
# Depending on type of gene, I will either accept deletions or amplifications
# From Bailey et al -> Oncogenes -> AMPs and not frame shift SNVs (others can convey a LOF but these most likely not)
# From Bailey et al -> TSG -> DELs and all SNVs
# Otherwise -> potential TSG as I added them because they have a cellular function -> DELs and all SNVs
allGenes = unique(snv$gene)


## Oncogenes (Potentially refined with transcriptomic data)
oncogenes = dtGeneLookup$name[ dtGeneLookup$Status == "Oncogene" ]

# For reviews add CCNE1 and rerun this script
if(REVIEW) { oncogenes = c(oncogenes, "CCNE1")}

mutsOnc = snv[ snv$gene %in% oncogenes & ( ! grepl("Frame_Shift", snv$consequence) ), ]
ampsOnc = amps[ amps$gene %in% oncogenes, ]

snvONCS = rbind(mutsOnc, ampsOnc)


## TSGs - treat every other gene as tsg
tsgs = dtGeneLookup$name[ dtGeneLookup$Status != "Oncogene" ]

# For reviews add CCNE1 and rerun this script
if(REVIEW) { tsgs = tsgs[ tsgs != "CCNE1" ] }

mutsTSG = snv[ snv$gene %in% tsgs, ]
delsTSG = dels[ dels$gene %in% tsgs, ]

snvDELS = rbind(mutsTSG, delsTSG)


## Merge oncogenes and tsgs
muts = rbind(snvONCS, snvDELS)

## Save output muts
write.table(muts, file.path(OUTTABLE, "Mutations_SNVs_AMPs_DELs_in_TCGA.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(muts, file.path(OUTTABLE, "Mutations_SNVs_AMPs_DELs_in_TCGA.rds"))

## Save summary table
mutSumm = unclass(table(muts$gene, muts$consequence))
write.table(mutSumm, file.path(OUTTABLE, "Mutations_SNVs_AMPs_DELs_in_TCGA_Summary.txt"), 
            sep = "\t", quote = FALSE, row.names = TRUE)
saveRDS(mutSumm, file.path(OUTTABLE, "Mutations_SNVs_AMPs_DELs_in_TCGA_Summary.rds"))




#### PART A: fisher test drivers

scExp = scale(raw, center = TRUE, scale = TRUE)
dtAct = data.table(melt(scExp))
colnames(dtAct) = c("Samples", "Signature", "Exposure")

## Using muts and dtSigQuant from above.
allDrivers = unique(muts$gene)
allSigs = levels(dtAct$Signature)

# Only doing assessments for samples for which we actually have coding mutation information.
# Caveat: Not every sample will have a SNV with high impact.
# Assumption: Every sample has at least one non-synonymous mutation
dtSigMuts = dtAct[ dtAct$Samples %in% unique(muts$sample), ]

## Add cancer type (only needed for additional analyses)
# dtSigMuts$Cancer = meta$cancer_type[ match(dtSigMuts$Samples, meta$name) ]

lSignatureTests = lapply(allSigs, function(sig) {
  
  print(sig)
  
  thisSig = dtSigMuts[ dtSigMuts$Signature == sig, ]
  
  ## For all normal drivers
  lAllDrivers = lapply(allDrivers, function(driver) { 
    
    # Transfer mutation status
    thisSig$driver = FALSE
    mutatedSamples = unique(muts$sample[ muts$gene == driver ])
    # Add prior knowledge
    if(driver == "BRCA1") {  mutatedSamples = unique(c(mutatedSamples, brca1)) }
    if(driver == "BRCA2") {  mutatedSamples = unique(c(mutatedSamples, brca2)) }
    if(driver == "RAD51C") {  mutatedSamples = unique(c(mutatedSamples, rad51c)) }
    
    thisSig$driver[ thisSig$Samples %in% mutatedSamples ] = TRUE
    thisSig$driver = as.logical(thisSig$driver)
    
    # For documentation
    # cohortProp = signif(mean(thisSig$driver), 4)
    # zeroProp = signif(mean(thisSig$driver[ thisSig$Exposure == 0 ]), 4)
    
    # t-Test above zero exposures
    positive = thisSig$Exposure[ thisSig$driver ]
    negative = thisSig$Exposure[ ! thisSig$driver ]
    meanPos = ifelse(length(positive) > 0, signif(mean(positive), 4), 0)
    meanNeg = signif(mean(negative), 4)
    meanDiff = signif(meanPos-meanNeg, 4)
    if(length(positive) > 1) { 
      ## T-test between groups
      tTest = t.test(positive, negative, var.equal = FALSE)
      tTestPval = signif(tTest$p.value, 4)
      
      ## Multivariate regression
      # thisSig$driver = factor(thisSig$driver, levels = c("FALSE", "TRUE"))
      # lmModel = lm(Exposure ~ 1 + driver + Cancer, thisSig)
      # meanDiff = signif(summary(lmModel)$coefficients["driverTRUE",1], 4)
      # tTestPval = signif(summary(lmModel)$coefficients["driverTRUE",4],)
      
      
    } else {
      tTestPval = 1
    }
    
    out = c(sig, driver, length(positive), length(negative), meanPos, 
            meanNeg, meanDiff, tTestPval)
    
    # out = c(sig, driver, zeroProp, cohortProp, signif(zeroProp-cohortProp, 4), 
    # signif(fTest$p.value, 4), meanPos, meanNeg, signif(meanPos-meanNeg, 4), tTestPval)
    return(out)
    
    
  } )
  
  
  ### Two additional tests regarding BRCA1 whether it was rescued or not
  ## For BRCA1 without TB53BP1, RIF1 and MAD2L2 mutation
  # Transfer mutation status
  thisSig$driver = FALSE
  brcaHit = unique(c(muts$sample[muts$gene == "BRCA1"], brca1))
  rescueHit = unique(muts$sample[muts$gene %in% c("TB53BP1", "RIF1", "MAD2L2")])
  thisSig$driver[ thisSig$Samples %in% brcaHit & 
                    ! ( thisSig$Samples %in% rescueHit ) ] = TRUE
  thisSig$driver = as.logical(thisSig$driver)
  
  # t-Test above zero exposures
  positive = thisSig$Exposure[ thisSig$driver ]
  negative = thisSig$Exposure[ ! thisSig$driver ]
  meanPos = ifelse(length(positive) > 0, signif(mean(positive),4), 0)
  meanNeg = signif(mean(negative),4)
  if(length(positive) > 1) { 
    tTest = t.test(positive, negative, var.equal = FALSE) 
    tTestPval = signif(tTest$p.value,4)
  } else {
    tTestPval = 1
  }
  
  out = c(sig, "BRCA1_noRescue", length(positive), length(negative), meanPos, meanNeg, 
          signif(meanPos-meanNeg,4), tTestPval)
  lAllDrivers[[length(lAllDrivers)+1]] = out
  
  ## For BRCA1 with TP53BP1, RIF1 or MAD2L2
  # Transfer mutation status
  thisSig$driver = NULL
  thisSig$driver = FALSE
  thisSig$driver[ thisSig$Samples %in% brcaHit & thisSig$Samples %in% rescueHit ] = TRUE
  thisSig$driver = as.logical(thisSig$driver)
  
  # t-Test above zero exposures
  positive = thisSig$Exposure[ thisSig$driver ]
  negative = thisSig$Exposure[ ! thisSig$driver ]
  meanPos = ifelse(length(positive) > 0, signif(mean(positive),4), 0)
  meanNeg = signif(mean(negative),4)
  if(length(positive) > 1) { 
    tTest = t.test(positive, negative, var.equal = FALSE) 
    tTestPval = signif(tTest$p.value,4)
  } else {
    tTestPval = 1
  }
  
  ## Note: The t.test is not significant but if I remove the one outlier 
  ## from the "positives" the test becomes significant. A wilcox.test is 
  ## slightly significant but plus p-value correction it is not anymore.
  ## Too little data to call it.
  out = c(sig, "BRCA1_rescued", length(positive), length(negative), meanPos, meanNeg, 
          signif(meanPos-meanNeg,4), tTestPval)
  lAllDrivers[[length(lAllDrivers)+1]] = out
  
  
  ## Wrap up
  dtDriverTests = data.table(do.call(rbind, lAllDrivers))
  colnames(dtDriverTests) = c("Signature", "Driver", "SampsMut", "SampsNonMut", 
                              "MeanDriverPos", "MeanDriverNeg", "MeanDiff", "pValTTest")
  dtDriverTests$pAdjTTest = p.adjust(as.numeric(dtDriverTests$pValTTest), method = "BH")
  
  return(dtDriverTests)
  
} )

dtDriverExps = rbindlist(lSignatureTests)
colnames(dtDriverExps) = c("Signature", "Driver", "SampsMut", "SampsNonMut", "MeanDriverPos",
                           "MeanDriverNeg", "MeanDiff", "pValTTest", "pAdjTTest")
dtDriverExps$SampsMut = as.numeric(dtDriverExps$SampsMut)
dtDriverExps$SampsNonMut = as.numeric(dtDriverExps$SampsNonMut)
dtDriverExps$MeanDriverPos = as.numeric(dtDriverExps$MeanDriverPos)
dtDriverExps$MeanDriverNeg = as.numeric(dtDriverExps$MeanDriverNeg)
dtDriverExps$MeanDiff = as.numeric(dtDriverExps$MeanDiff)
dtDriverExps$pValTTest = as.numeric(dtDriverExps$pValTTest)
dtDriverExps$pAdjTTest = as.numeric(dtDriverExps$pAdjTTest)

## Save raw table
write.table(x = dtDriverExps, 
            file = file.path(OUTTABLE,
                             "Table_tTest_Genes_of_Interest_THRESH95_NAMESAPRIL2021.txt"), 
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
saveRDS(dtDriverExps, file.path(OUTTABLE,
                                "Table_tTest_Genes_of_Interest_THRESH95_NAMESAPRIL2021.rds"))


if(REVIEW) {
  dtCCNE = dtDriverExps[ dtDriverExps$Driver == "CCNE1", ]
  write.table(x = dtCCNE, 
              file = file.path(OUTTABLE,
                               "Table_tTest_CCNE1_only.txt"), 
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  
  dtCCNE[ dtCCNE$MeanDiff > 0.4 & dtCCNE$pAdjTTest < 0.005, ]
}


## Save two summary tables - all positives significant and after special filtering
dtOut = dtDriverExps[dtDriverExps$pAdjTTest < 0.05,]
dtOutPrint = dtOut[ dtOut$MeanDiff > 0, ]

write.table(x = dtOutPrint, 
            file = file.path(OUTTABLE,
                             "Table_tTest_Test_GoI_Enrichment_THRESH95_NAMESAPRIL2021.txt"), 
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
saveRDS(dtOutPrint, file.path(OUTTABLE, 
                              "Table_tTest_Test_GoI_Enrichment_THRESH95_NAMESAPRIL2021.rds"))


dtFilt = dtOutPrint[ dtOutPrint$MeanDiff > 0.4 & dtOutPrint$pAdjTTest < 0.005, ]
write.table(dtFilt, file.path(OUTTABLE,
                              "Table_tTest_Test_GoI_Enrichment_SigFilt_THRESH95_NAMESAPRIL2021.txt"), 
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
saveRDS(dtFilt, file.path(OUTTABLE,
                          "Table_tTest_Test_GoI_Enrichment_SigFilt_THRESH95_NAMESAPRIL2021.rds"))





#### PART B: Multivariate regression

## Same process just with cancer type correction => For review process only

## Add cancer type (only needed for additional analyses)
dtSigMuts$Cancer = meta$cancer_type[ match(dtSigMuts$Samples, meta$name) ]
lSignatureTestsCorr = lapply(allSigs, function(sig) {
  
  print(sig)
  
  thisSig = dtSigMuts[ dtSigMuts$Signature == sig, ]
  
  ## For all normal drivers
  lAllDrivers = lapply(allDrivers, function(driver) { 
    
    # Transfer mutation status
    thisSig$driver = FALSE
    mutatedSamples = unique(muts$sample[ muts$gene == driver ])
    # Add prior knowledge
    if(driver == "BRCA1") {  mutatedSamples = unique(c(mutatedSamples, brca1)) }
    if(driver == "BRCA2") {  mutatedSamples = unique(c(mutatedSamples, brca2)) }
    if(driver == "RAD51C") {  mutatedSamples = unique(c(mutatedSamples, rad51c)) }
    
    thisSig$driver[ thisSig$Samples %in% mutatedSamples ] = TRUE
    thisSig$driver = as.logical(thisSig$driver)
    
    # Gather stats
    positive = thisSig$Exposure[ thisSig$driver ]
    negative = thisSig$Exposure[ ! thisSig$driver ]
    if(length(positive) > 1) { 
      
      ## Multivariate regression
      thisSig$driver = factor(thisSig$driver, levels = c("FALSE", "TRUE"))
      lmModel = lm(Exposure ~ 1 + driver + Cancer, thisSig)
      meanDiff = signif(summary(lmModel)$coefficients["driverTRUE",1], 4)
      tTestPval = signif(summary(lmModel)$coefficients["driverTRUE",4],)
    } else {
      meanDiff = 0
      tTestPval = 1
    }
    
    out = c(sig, driver, length(positive), length(negative), meanDiff, tTestPval)
    return(out)
  } )
  
  ## Wrap up
  dtDriverTests = data.table(do.call(rbind, lAllDrivers))
  colnames(dtDriverTests) = c("Signature", "Driver", "SampsMut", "SampsNonMut", 
                              "MeanDiff", "pValTTest")
  return(dtDriverTests)
  
} )

dtDriverExpsCorr = rbindlist(lSignatureTestsCorr)
colnames(dtDriverExpsCorr) = c("Signature", "Driver", "SampsMut", "SampsNonMut", 
                               "MeanDiff", "pValLM")
dtDriverExpsCorr$SampsMut = as.numeric(dtDriverExpsCorr$SampsMut)
dtDriverExpsCorr$SampsNonMut = as.numeric(dtDriverExpsCorr$SampsNonMut)
dtDriverExpsCorr$MeanDiff = as.numeric(dtDriverExpsCorr$MeanDiff)
dtDriverExpsCorr$pValLM = as.numeric(dtDriverExpsCorr$pValLM)


## Filter and correct p-values
dtDriverFiltCorr = dtDriverExpsCorr[ dtDriverExpsCorr$MeanDiff > 0.2, ]
dtDriverFiltCorr$pAdjLM = p.adjust(as.numeric(dtDriverFiltCorr$pValLM), method = "BH")
dtDriverSigCorr = dtDriverFiltCorr[ dtDriverFiltCorr$pAdjLM < 0.05, ]

write.table(dtDriverSigCorr, file.path(OUTTABLE,
                                       "Table_tTest_Test_GoI_Enrichment_SigFilt_THRESH95_NAMESAPRIL2021_Multivar.txt"), 
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


## Which ones are not in the first list
lMissed = lapply(unique(dtFilt$Signature), function(sig) {
  
  dtOld = dtFilt[ dtFilt$Signature == sig, ]
  dtNew = dtDriverSigCorr[ dtDriverSigCorr$Signature == sig, ]
  
  dtMiss = dtOld[ ! dtOld$Driver %in% dtNew$Driver, ]
  return(dtMiss)
} )

dtMissed = rbindlist(lMissed)

write.table(dtMissed, file.path(OUTTABLE,
                                "Table_tTest_Test_GoI_Enrichment_SigFilt_THRESH95_NAMESAPRIL2021_Multivar_Missed.txt"), 
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)





#### PART C: gather data for gene-by-signature plot

## Rerun a few tests for supplementary figures showing the signature activity per significant gene
## Date: 18/05/21

## Exposures => scaled raw values => dtSigMuts from above
IMPORTANTGENES = c("JAK2", "POLE", "AKT1", "BRAF", "CCND1", "CDK4", "CUL1", "EGFR", 
                   "ERBB2", "ERBB3", "ERCC2", "H3F3A", "IDH1", "IDH2", "KRAS", "MAPK1", 
                   "MYC", "MYCN", "NFE2L2", "PCBP1", "PIK3R2", "PMS1", "PPP2R1A", 
                   "RAC1", "SOS1", "SPOP", "U2AF1", "ARID1A", "ATM", "BRCA1", "BRCA1 rescued", 
                   "BRCA2", "CDKN2C", "CIC", "FBXW7", "NF1", "NOTCH1", "PBRM1", 
                   "RPL22", "SMAD4", "SOX9", "TP53", "VHL", "PSMA4")


lSignatureTests = lapply(allSigs, function(sig) {
  
  print(sig)
  
  thisSig = dtSigMuts[ dtSigMuts$Signature == sig, ]
  
  ## For all normal drivers
  lAllDrivers = lapply(IMPORTANTGENES, function(driver) { 
    
    # Transfer mutation status
    thisSig$Driver = FALSE
    mutatedSamples = unique(muts$sample[ muts$gene == driver ])
    # Add prior knowledge
    if(driver == "BRCA1") {  mutatedSamples = unique(c(mutatedSamples, brca1)) }
    if(driver == "BRCA2") {  mutatedSamples = unique(c(mutatedSamples, brca2)) }
    if(driver == "RAD51C") {  mutatedSamples = unique(c(mutatedSamples, rad51c)) }
    
    thisSig$Driver[ thisSig$Samples %in% mutatedSamples ] = TRUE
    thisSig$Driver = as.logical(thisSig$Driver)
    
    thisSig$Gene = driver
    return(thisSig)
    
  } )
  
  ## For BRCA1 with TP53BP1, RIF1 or MAD2L2
  # Transfer mutation status
  thisSig$Driver = FALSE
  
  brcaHit = unique(c(muts$sample[muts$gene == "BRCA1"], brca1))
  rescueHit = unique(muts$sample[muts$gene %in% c("TB53BP1", "RIF1", "MAD2L2")])
  thisSig$Driver[ thisSig$Samples %in% brcaHit & thisSig$Samples %in% rescueHit ] = TRUE
  thisSig$Driver = as.logical(thisSig$Driver)
  
  thisSig$Gene = "BRCA1 rescued"
  lAllDrivers[[length(lAllDrivers)+1]] = thisSig
  
  
  dtDriverTests = data.table(do.call(rbind, lAllDrivers))
  return(dtDriverTests)
  
} )

dtSigDrivers = rbindlist(lSignatureTests)

saveRDS(dtSigDrivers, file.path(OUTTABLE, "Mutations_and_GoI_for_Supplementary_Figure.rds"))


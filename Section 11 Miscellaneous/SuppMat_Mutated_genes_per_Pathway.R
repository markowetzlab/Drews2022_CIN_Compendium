## Output fraction of samples with mutation

rm(list=ls(all=TRUE))

library(this.path)
library(data.table)
library(reshape2)
# For outptutting excel files -> Easier as this goes into supplement
library(writexl)


## Paths
BASE=dirname(this.path())
OUTTABLE=file.path(BASE, "output")


## Files
ACT=file.path(BASE, "input/Signature_Compendium_v5_Cosine-0.74_Activities_THRESH95_NAMESAPRIL21.rds")
MUTS=file.path(BASE, "input/Mutations_SNVs_AMPs_DELs_in_TCGA.rds")
EXTRABRCA=file.path(BASE, "input/Survival_and_BRCA_Status_fullTCGA.rds")

## Which genes to check
lCheck = list("CX2" = c("BRCA1", "BRCA2", "RAD51C"), "CX3" = c("BRCA1", "BRCA2", "RAD51C"), 
     "CX4" = c("PIK3R2", "AKT1", "PIK3CA"), "CX5" = c("BRCA1", "BRCA2", "RAD51C"),
     "CX6" = c("PSMA4", "CUL1", "RAC1"), "CX10" = c("FBXW7"))
INTERVALS=seq(0.1, 1, by = 0.1)


## Load data
dtAct = data.table(melt(readRDS(ACT)))
colnames(dtAct) = c("Sample", "Sig", "Activity")
dtMuts = readRDS(MUTS)
dtExtra = readRDS(EXTRABRCA)


# Check that all gBRCA1/2 and RAD51C hypermeth is listed correctly
for(thisGene in c("BRCA1", "BRCA2", "RAD51C")) {
  
  # Identify samples with extra mutational information
  allBRCA1 = dtExtra$Name[ (! grepl("WT", dtExtra$Status)) &  
                               (grepl(thisGene, dtExtra$Status)) ]
  # Convert to same format as muts data table.
  dtBRCA1 = data.table(data.frame(thisGene, allBRCA1, "Extra_mutation"))
  colnames(dtBRCA1) = colnames(dtMuts)
  # Add them to the muts data table. 
  # Duplicates are no problem as names are filtered later.
  dtMuts = rbind(dtMuts, dtBRCA1)
}


## Assign intervals to samples
lAct = split(dtAct, dtAct$Sig)
lAssign = lapply(lAct, function(thisSig) {
  
  # Identify breakpoints according to user-supplied intervals
  vBreaks = c(0, 1e-16, quantile(thisSig$Activity[ thisSig$Activity > 0 ], probs = INTERVALS))
  vLabels = c("Inactive", paste0("<", INTERVALS*100, "%"))
  
  thisSig$Quant = cut(thisSig$Activity, breaks = vBreaks, labels = vLabels, include.lowest = TRUE)
  return(thisSig)
  
} )


## Loop over signatures and pathways and count occurrences
lSigsAndGenes = lapply(names(lCheck), function(thisSig) {
  
  print(thisSig)
  
  ## Prepare data
  theseAct = lAssign[[ thisSig ]]
  theseGenes = lCheck[[ thisSig ]]
  
  ## Label samples with mutation in gene
  theseAct$GoI = "WT"
  for(thisGene in theseGenes) {
    
    mutatedSamples = unique(dtMuts$sample[ dtMuts$gene %in% thisGene ])
    theseAct$GoI[ theseAct$Sample %in% mutatedSamples ] = thisGene
    
  }
  
  ## Factorise to preserve order
  theseAct$GoI = factor(theseAct$GoI, levels = c(theseGenes, "WT"))
  
  ## Count instances and format for easy xlsx export
  dfTab = as.data.frame.matrix(table(theseAct$GoI, theseAct$Quant))
  dfTab = cbind("Gene" = rownames(dfTab), dfTab)
  
  ## Remove WT because will be added again in a second
  dfTab = dfTab[ dfTab$Gene != "WT", ]
  
  

  ## Count pathway and WT
  theseAct$Pathway = "WT"
  theseAct$Pathway[ theseAct$GoI != "WT" ] = "Pathway"
  
  # Count instances and format for easy xlsx export
  dfTab2 = as.data.frame.matrix(table(theseAct$Pathway, theseAct$Quant))
  dfTab2 = cbind("Gene" = rownames(dfTab2), dfTab2)
  
  dfOut = rbind(dfTab, dfTab2)
  
  
  ## Add summary column of active samples
  countsPerGene = table(theseAct$GoI[ theseAct$Activity > 0 ])
  countsPerGene = countsPerGene[ names(countsPerGene) != "WT" ]
  Active = c(countsPerGene, 
           sum(dfOut[ dfOut$Gene == "Pathway", c(-1,-2) ]), 
           sum(dfOut[ dfOut$Gene == "WT", c(-1,-2) ]))
  dfOut = cbind(dfOut, Active)
  
  ## Reorder
  vLabels = c("Gene", "Inactive", "Active", paste0("<", INTERVALS*100, "%"))
  dfOut = dfOut[, vLabels]
  
  ## Add proportions in percent
  dfProps = signif(dfOut["Pathway",-1]/(dfOut["Pathway",-1]+dfOut["WT",-1]), 2)
  dfProps = cbind("Gene" = "Mutant proportion", dfProps)
  dfReturn = rbind(dfOut, dfProps)
  
  return(dfReturn)
})

names(lSigsAndGenes) = names(lCheck)
write_xlsx(lSigsAndGenes,file.path(OUTTABLE, "SuppMat_Mutated_genes_per_Pathway.xlsx"), format_headers = TRUE)


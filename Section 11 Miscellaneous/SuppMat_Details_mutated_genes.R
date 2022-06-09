## Merge mutated genes per signature with number of mutations per mutation type

rm(list=ls(all=TRUE))

library(this.path)
library(data.table)
library(writexl)

# Paths and loading files
BASE=dirname(this.path())
OUTTABLE=file.path(BASE, "output")
dtSig=readRDS(file.path(BASE, "input/Table_tTest_Test_GoI_Enrichment_SigFilt_THRESH95_NAMESAPRIL2021.rds"))
dtSig2=fread(file.path(BASE, "input/Table_tTest_Test_GoI_Enrichment_SigFilt_THRESH95_NAMESAPRIL2021_Multivar.txt"))
dtMuts=readRDS(file.path(BASE, "input/Mutations_SNVs_AMPs_DELs_in_TCGA_Summary.rds"))


# Merge
dtMerge = cbind(dtSig, dtMuts[ match(dtSig$Driver, rownames(dtMuts)), ])

# Reorder
dtMerge$Signature = factor(dtMerge$Signature, levels = paste0("CX", 1:17))
dtMerge = dtMerge[ order(dtMerge$Signature, -dtMerge$MeanDiff), ]

write_xlsx(dtMerge,file.path(OUTTABLE, "SuppMat_Mutated_Details_mutated_genes_tTest.xlsx"), format_headers = TRUE)



# Merge #2
dtMerge2 = cbind(dtSig2, dtMuts[ match(dtSig2$Driver, rownames(dtMuts)), ])

# Reorder #2
dtMerge2$Signature = factor(dtMerge2$Signature, levels = paste0("CX", 1:17))
dtMerge2 = dtMerge2[ order(dtMerge2$Signature, -dtMerge2$MeanDiff), ]

write_xlsx(dtMerge2,file.path(OUTTABLE, "SuppMat_Mutated_Details_mutated_genes_lm.xlsx"), format_headers = TRUE)



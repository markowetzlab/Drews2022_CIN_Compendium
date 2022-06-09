# Convert PCAWG SV data into one file with added microhomology identification

library(vcfR)
library(data.table)

## Please download structural variants files for the PCAWG project from the 
## respective data portal.
SVFOLDER="/Path/To/Your/PCAWGSV/Folder"

## Read in VCFs and identify distribution of microhomologies for each sample
allVCFs = list.files(SVFOLDER, pattern = "\\.vcf\\.gz$", full.names = TRUE)
lHomLen = lapply(allVCFs, function(thisVCFFile) {
  
  ## Sample name
  sampleName = strsplit(basename(thisVCFFile), "\\.")[[1]][1]
  print(sampleName)
  
  thisVCF = read.vcfR(thisVCFFile, verbose = FALSE )
  if(nrow(thisVCF) == 0) { return(NULL) }
  ## Only info we need is in HOMLEN but for future proof we use the whole INFO part of the vcf
  ## Too slow
  # dfInfo = INFO2df(thisVCF)
  tbInfo = vcfR2tidy(thisVCF, info_only = TRUE)$fix
  
  
  ## No info available from snowman about how HOMSEQ or HOMLEN were obtained and what the values mean
  ## Especially the value NA.
  
  ## Assumption: Length 0 when HOMSEQ==NA and HOMLEN==NA
  ## What about HOMLEN==NA but HOMSEQ!=NA?
  ## Works only with tibble as tibble converst string NAs into ""
  tbInfo$HOMLEN[ is.na(tbInfo$HOMLEN) & tbInfo$HOMSEQ == "" ] = 0
  
  ## Also carry over NAs (probably data not good enough to determine MHs) for getting number of SVs
  dtOut = data.table("Sample" = sampleName, "HOMLEN" = tbInfo$HOMLEN)
  return(dtOut)
  
})

dtHom = rbindlist(lHomLen)
saveRDS(dtHom, "~/CINsignatures/2_intermediate_files/PCAWG_Microhomologies_at_SVs.rds")

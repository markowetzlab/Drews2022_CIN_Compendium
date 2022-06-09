## This script can be run after signatures have been derived from the sum-of-posterior matrix.

rm(list=ls(all=TRUE))

library(this.path)
library(ggplot2)
library(reshape2)
library(patchwork)
library(ggthemes)
library(lemon)

theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))


NORMPERFEATURE=TRUE

#### Load solutions
thisPath=dirname(this.path())
INPUT=file.path(thisPath, "input")
OUT=file.path(thisPath, "output")
# Manually adjust wich signature file should be loaded.
mDefs = readRDS(file.path(INPUT, "Out_6_Signatures_XK5TRD_normalised.rds"))
NAME=file.path(OUT, "SuppFig_4_Signature_defs")

# Per signature and per feature
FEATS=c("segsize", "changepoint", "bp10MB", "bpchrarm", "osCN")
NUMCOMP = ncol(mDefs)

if(NORMPERFEATURE) {
  
  ## First normalise per signature
  theseNewSigs = apply(mDefs, 1, function(thisSig) {
    
    lSig = sapply(FEATS, function(thisFeat) {
      
      theseVals = thisSig[ grepl(thisFeat, names(thisSig)) ]
      theseNew = theseVals / sum(theseVals)
      # Catch edge case with only zeros and no weights (produce NaN)
      if(is.nan(sum(theseNew))) {
        altNew = rep(0, length(theseNew))
        names(altNew) = names(theseNew)
        theseNew = altNew
      }
      
      # Scale by numbers of components => Final sum of vector should be five (for five features)
      theseNew = theseNew * (length(theseNew)/NUMCOMP)
    } )
    
    vSig = unlist(lSig)
    return(vSig)
    
  } )
  
  mDefs = t(theseNewSigs)
  colnames(mDefs) = sapply(strsplit(colnames(mDefs), "\\."), function(x) x[[2]])
  
}


# Sort manually to match order of simulated sigs
# Ruben Param 2 plus late WGD (2x)
# N=240 K=5 (optimal K) Repr. Solution
mDefs = mDefs[c(3,5,1,4,2),]
rownames(mDefs) = c("HRD/LST", "ECDNA", "CHR", "WGD early", "WGD late")

dfMelt = melt(t(mDefs))
colnames(dfMelt)<-c("component", "sig", "weight")

dfMelt$feature = gsub('[0-9]+', '', dfMelt$component)
dfMelt$feature[ dfMelt$feature == "bpMB" ] = "bp10MB"
dfMelt$feature = factor(dfMelt$feature, levels = c("segsize", "changepoint",
                                                   "bp10MB", "bpchrarm",
                                                   "osCN"))

p0 = ggplot(dfMelt, aes(x = component, y = weight, fill = feature))+
  geom_bar(stat = "identity") + facet_wrap(. ~ sig, ncol = 1, scales = "free_y") + 
  coord_capped_cart(left = "both", bottom = "both") +
  labs(x = "Feature components", y = "Identified weight") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none") 

ggsave(paste0(NAME, ".png"), p0, width = 12, height = 10, units = "cm", bg = "white")

ggsave(paste0(NAME, ".svg"), p0, width = 12, height = 10, units = "cm", bg = "white")

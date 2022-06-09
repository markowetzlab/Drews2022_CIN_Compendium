## Plot signature definitions in barplot style

rm(list=ls(all=TRUE))

library(this.path)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(lemon)

theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))


## Paths
BASE=dirname(this.path())
OUT=file.path(BASE, "output")
## At least one should be false:
# FALSE, TRUE => Normalised per sig and then per feat
# TRUE, FALSE => Normalised per component
# FALSE, FALSE => Normalised per signature
# TRUE, TRUE => Normalised per component and then per feat (funny?)
NORMALISEPERCOMP = FALSE
NORMALISEPERFEAT = TRUE 
DEFS=file.path(BASE, "input/Signature_Compendium_v5_Cosine-0.74_Signatures_NAMESAPRIL21.rds")

# Load and prepare data
mDefs = readRDS(DEFS)
if(NORMALISEPERCOMP) { mDefs = apply(mDefs, 2, function(x) x/sum(x)) }

# Per signature and per feature
FEATS=c("segsize", "changepoint", "bp10MB", "bpchrarm", "osCN")
NUMCOMP = ncol(mDefs)
if(NORMALISEPERFEAT) {
  
  ## First normalise per signature (to remove the comparable signature)
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


# Melt to data.frame
dtDefs = data.table(melt(mDefs))


# Prepare plot and produce the figure
dtDefs$Col = gsub('[[:digit:]]+', '', dtDefs$Var2 )
dtDefs$Col = factor(dtDefs$Col, levels = c("segsize", "changepoint", "bpMB", "bpchrarm", "osCN"),
                    labels = c("Segment size", "Change point", "Breakpoints per 10MB", 
                               "Breakpoints per arm", "Lengths of chains of osc. CN"))

p1 = ggplot(dtDefs, aes(x = Var2, y = value, fill = Col)) + geom_col() + facet_wrap(. ~ Var1, ncol = 2) +
  scale_fill_manual(values = c("#ffbe0b", "#fb5607", "#ff006e", "#8338ec", "#3a86ff")) +
  scale_x_discrete(guide = guide_axis(angle = 90)) + 
  labs(x = "Feature components", y = "Weights", fill = "Feature") +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom") + 
  coord_capped_cart(left = "both", bottom = "both") +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ggsave(file.path(OUT, "SuppFig_2_Definitions_Barplot.svg"), p1, width = 180, height = 220, units = "mm", bg = "white")
ggsave(file.path(OUT, "SuppFig_2_Definitions_Barplot.png"), p1, width = 180, height = 220, units = "mm", bg = "white")


#### Test whole chromosomes and sig 1 exposure
rm(list=ls(all=TRUE))

library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(lemon)
library(this.path)

theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5), 
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))

BASE=dirname(this.path())
CN=file.path(BASE, "input/0_TCGA_Segments_dCIN.rds")
EXP=file.path(BASE, "input/Signature_Compendium_v5_Cosine-0.74_Activities_THRESH95_NAMESAPRIL21.rds")

## These are the canonical chromosome lengths. But they may be too long for SNP6.0 arrays.
LENGTHS=file.path(BASE, "input/hg19.chrom.sizes.txt")
WIGGLE = 0.05

OUTTABLE=file.path(BASE, "input")

## Load data
cn=readRDS(CN)
cn$lengths = cn$end - cn$start

exp=readRDS(EXP)
dtExp = data.table(melt(exp))

# Prepare canonical lengths
lens=fread(LENGTHS)
lens = lens[1:24,]
lens$V1 = substr(lens$V1, 4, stop = 6)

# Extract longest chromosome lengths present in data set
chrLens = aggregate(lengths ~ chromosome, data = cn, max)
chrLens$hg19 = lens$V2[ match(chrLens$chromosome, lens$V1) ]
chrLens$test = chrLens$lengths > (1-WIGGLE)*chrLens$hg19

# For all chromosomes except chromosomes 13, 14, 15, 21, 22 the lengths of the array match roughly the lengths
# of the canonical chromosomes. Because there is a bit of a discrepancy, use the array lengths
# to check the dataset for whole chromosomes.
cn$chrLength = (1-WIGGLE) * chrLens$lengths[ match(cn$chromosome, chrLens$chromosome) ]
cn$whole = cn$lengths > cn$chrLength

wholeChrCNAs = cn[ cn$segVal != 2 & cn$whole, ]
dtWhole = data.table(aggregate(whole ~ sample, data = wholeChrCNAs, sum, drop = FALSE))

# Zeros became NAs
dtWhole$whole[ is.na(dtWhole$whole) ] = 0

hist(dtWhole$whole, breaks = 20)

# Save output
saveRDS(dtWhole, file.path(OUTTABLE, "Whole_chromosome_CNAs.rds"))


#### Correlation done here but plotted with covariates in "SuppFig_Covariates.R"

## Correlate number of whole chromosomes for all signatures
## Transfer sig 1 exposure and test correlation.
allSigs = levels(dtExp$Var2)
lSpear = lapply(allSigs, function(thisSig) {
  
  dtCS1 = dtExp[ dtExp$Var2 == thisSig, ]
  dtWhole$sig1 = dtCS1$value[ match(dtWhole$sample, dtCS1$Var1) ]
  spTest = cor.test(dtWhole$whole, dtWhole$sig1, method = "spearman", use = "pairwise.complete.obs")
  out = c(thisSig, spTest$estimate, spTest$p.value)
  return(out)
})

dfSpear = data.frame(do.call(rbind, lSpear))
colnames(dfSpear) = c("Sig", "Rho", "pVal")
dfSpear$Rho = as.numeric(as.character(dfSpear$Rho))
dfSpear$pVal = as.numeric(as.character(dfSpear$pVal))
dfSpear$pAdj = p.adjust(dfSpear$pVal)

dfSpear

saveRDS(dfSpear, file.path(OUTTABLE, "Whole_chromosome_CNAs_Spearman.rds"))

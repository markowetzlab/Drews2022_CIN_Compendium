## Compare hypothesised and identified signatures

rm(list=ls(all=TRUE))

library(this.path)
library(data.table)
library(reshape2)
library(lsa)

## Paths
thisPath=dirname(this.path())
INPUT=file.path(thisPath, "input")
OUT=file.path(thisPath, "output")


# Taken from PharmacoGx (it's just a very big package and I only need a small function) 
cosinePerm = function (x, y, nperm = 1000, alternative = c("two.sided", "less", 
                                                           "greater"), include.perm = FALSE, setseed = 12345, nthread = 1) 
{
  set.seed(setseed)
  alternative <- match.arg(alternative)
  if ((length(x) != length(y))) {
    stop("x and y must be vectors of the same length")
  }
  res <- c(estimate = NA, p.value = NA)
  x <- as.numeric(x)
  y <- as.numeric(y)
  res["estimate"] <- drop(lsa::cosine(x = x, y = y))
  if (nperm > 0) {
    splitix <- parallel::splitIndices(nx = nperm, ncl = nthread)
    splitix <- splitix[sapply(splitix, length) > 0]
    mcres <- parallel::mclapply(splitix, function(x, xx, 
                                                  yy) {
      res <- sapply(x, function(x, xx, yy) {
        xx <- sample(xx)
        yy <- sample(yy)
        return(drop(lsa::cosine(x = xx, y = yy)))
      }, xx = xx, yy = yy)
      return(res)
    }, xx = x, yy = y)
    mcres <- unlist(mcres)
    switch(alternative, two.sided = {
      res["p.value"] <- 2 * (min(sum(mcres < res["estimate"]), 
                                 sum(mcres > res["estimate"]))/sum(!is.na(mcres)))
    }, less = {
      res["p.value"] <- sum(mcres < res["estimate"])/sum(!is.na(mcres))
    }, greater = {
      res["p.value"] <- sum(mcres > res["estimate"])/sum(!is.na(mcres))
    })
    if (res["p.value"] == 0) {
      res["p.value"] <- 1/(nperm + 1)
    }
  }
  res <- as.list(res)
  if (include.perm) {
    res <- c(res, list(estimate.random = mcres))
  }
  return(res)
}




# Prepare hypothesised sigs
vFeatComps = paste0(rep(c("segsize", "changepoint", "bp10MB", "bpchrarm", "osCN"), c(22, 10, 3, 5, 3)),
                    c(1:22, 1:10, 1:3, 1:5, 1:3))
dfSigs = data.frame("component" = factor(vFeatComps, levels = vFeatComps), 
                    # "feature" = factor(rep(c("segsize", "changepoint", "bp10MB", "bpchrarm", "osCN"), c(22, 10, 3, 5, 3))),
                    "HRDLST" = c(rep(0, 8), rep(1,7), rep(0,7), # Segsize
                                 c(0,0,1,1, rep(0, 6)), # Changepoint
                                 c(0,1,0), # bp10MB
                                 c(0,1,1,0,0), # bpchrarm
                                 c(1,0,1)), # osCN
                    "ECDNA" = c(rep(0,4), 1, 1, 1, rep(0, 15), # Segsize
                                c(rep(0, 8), c(1,1)), # Changepoint
                                c(0, 1, 1), # bp10MB
                                c(0, 1, 1, 0, 0), # bpchrarm
                                c(1,1,0)), # osCN
                    "CHR" = c(rep(0, 15), rep(1,7), # Segsize
                              c(0,0, 1, 1, rep(0, 6)), # Changepoint
                              c(1, 0, 0), # bp10MB
                              c(1, rep(0,4)), # bpchrarm
                              c(1,rep(0,2))), # osCN
                    "WGDearly" = c(rep(0, 15), rep(1,7), # Segsize
                                   c(rep(0,4),1,1,rep(0, 4)), # Changepoint
                                   c(1, 0, 0), # bp10MB
                                   c(1, rep(0,4)), # bpchrarm
                                   c(1,rep(0,2))),
                    "WGDlate" = c(rep(0, 4), rep(1,18), # Segsize
                                  c(rep(0,3),1,1,1,rep(0, 4)), # Changepoint
                                  c(0, 1, 1), # bp10MB
                                  c(0,1,1,0,0), # bpchrarm
                                  c(1, 1, 1))) # osCN
rownames(dfSigs) = dfSigs$component
dfSigs$component = NULL
mHypo = t(as.matrix(dfSigs))




### Load identified sigs
## Sig defs
# RENAME=TRUE
# mIdentified = mDefs = readRDS(file.path(OUT, "Out_6_Signatures_XK5TRD_normalised.rds"))
# SIGNAMES = c("W3"="HRDLST", "W5"="ecDNA", "W1"="CHR", "W4"="WGDearly", "W2"="WGDlate")

## Sig interpretation matrix
RENAME=FALSE
dfSigs = read.table(file.path(INPUT, "Out_8_InterpretationMatrix_simulations.txt"))
dfMelt = dfSigs[,c("component", "CNSig", "allSum")]
colnames(dfMelt)<-c("component", "sig", "weight")
# Order segments
dfMelt$component = factor(as.character(dfMelt$component), levels = c(paste0(rep("segsize", 22), 1:22), paste0(rep("changepoint", 10), 1:10),
                                                                     paste0(rep("bp10MB", 3), 1:3), paste0(rep("bpchrarm", 5), 1:5),
                                                                     paste0(rep("osCN", 3), 1:3)) )
# Order signatures
dfMelt$sig = factor(as.character(dfMelt$sig), levels = c("HRDLST", "ecDNA", "CHR", "WGDearly", "WGDlate"))
mIdentified=acast(dfMelt, sig~component, value.var="weight")
identical(colnames(mHypo), colnames(mIdentified))



# Rename signatures
if(RENAME){
  rownames(mIdentified) =SIGNAMES[ match(rownames(mIdentified), names(SIGNAMES)) ]
  mIdentified=mIdentified[SIGNAMES,]
}


# Do comparison
apply(mIdentified, 1, function(identified) {
  apply(mHypo, 1, function(hypo) {
    cosine(identified, hypo)
  })
})

lOut2 = lapply(1:nrow(mIdentified), function(i){

  outerName = rownames(mIdentified)[i]
  thisIdentified = mIdentified[i,]

  lOut = lapply(1:nrow(mHypo), function(j){

    innerName = rownames(mHypo)[j]
    thisHypo = mHypo[j,]

    res = cosinePerm(thisIdentified, thisHypo, nperm = 1e4, nthread = 4)
    out = c(outerName, innerName, res$estimate, res$p.value)
    return(out)
  })

  dfOut = as.data.frame(do.call(rbind, lOut))
  colnames(dfOut) = c("Identified", "Hypothesised", "Cosine", "pVal")
  dfOut$qVal = p.adjust(as.numeric(dfOut$pVal), method = "BH")
  return(dfOut)
})

dfOut = as.data.frame(do.call(rbind, lOut2))
write.table(dfOut, file.path(OUT, "SuppInfo_Genome_Simulation_Correlation_between_sigs.txt"), sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)

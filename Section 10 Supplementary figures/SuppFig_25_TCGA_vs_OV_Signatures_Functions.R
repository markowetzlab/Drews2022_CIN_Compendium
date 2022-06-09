## Functions
# Source: https://stackoverflow.com/questions/52297978/decrease-overal-legend-size-elements-and-text
addSmallLegend = function(myPlot, pointSize = 0.5, textSize = 3, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

## Taken from PharmacoGx (it's just a very big package and I only need a small function)
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

extractModelsFromFlexmix = function(thisModel) {
  
  allFeatures = names(thisModel)
  DTMODEL = list()
  for(thisFeature in allFeatures) {
    
    # Get OV data from flexmix objects
    thisOld = thisModel[[ thisFeature ]]
    
    if( sum(grepl("pois", thisOld@call)) > 0 ) {
      tableOld = data.frame( parameters(thisOld) )
      colnames(tableOld) = c("Mean")
    } else {
      tableOld = data.frame( t( data.frame( parameters(thisOld) ) ) )
      colnames(tableOld) = c("Mean", "SD")
    }
    
    tableOld$Weight = prior(thisOld)
    DTMODEL[[ thisFeature ]] = tableOld[ order(tableOld$Mean), ]
    
  }
  
  return(DTMODEL)
}

getConversionMatrix = function(MODELORIGIN, MODELTARGET, plotHeatmap = FALSE, uninfPrior = TRUE, noCN = FALSE) {
  
  # Convert to list
  if("flexmix" %in% class(MODELORIGIN[[1]])) MODELORIGIN = extractModelsFromFlexmix(MODELORIGIN)
  if("flexmix" %in% class(MODELTARGET[[1]])) MODELTARGET = extractModelsFromFlexmix(MODELTARGET)
  
  # For saving posterior probabilities
  BESTOVERLAP=list()
  # Just to make sure the rownames are the same afterwards.
  ROWNAMES=list()
  allFeatures = names(MODELTARGET)
  if(noCN) allFeatures = allFeatures[ grep("copynumber", allFeatures, invert = TRUE) ]
  for(thisFeature in allFeatures) {
    
    # Get posterior probability for mean of OV models based on TCGA models
    originDat = MODELORIGIN[[ thisFeature ]]$Mean
    targetMod = MODELTARGET[[ thisFeature ]]
    
    if( ncol(targetMod) == 2 ) {
      # Poisson model
      originDat = round(originDat)
      if(uninfPrior) {
        # Uninformative prior. As all weights would be the same, we can just drop it as we scale later anyways.
        postDatUnscaled = sapply(1:nrow(targetMod), function(x) dpois(x = originDat, lambda = targetMod[[x,"Mean"]]) )    
      } else {
        postDatUnscaled = sapply(1:nrow(targetMod), function(x) dpois(x = originDat, lambda = targetMod[[x,"Mean"]]) *
                                   targetMod[[x,"Weight"]])
      }
      
    } else {
      # Gaussian model
      if(uninfPrior) {
        # Uninformative prior. As all weights would be the same, we can just drop it as we scale later anyways.
        postDatUnscaled = sapply(1:nrow(targetMod), function(x) dnorm(x = originDat, mean = targetMod[[x,"Mean"]], 
                                                                      sd = targetMod[[x,"SD"]]) )
      } else {
        postDatUnscaled = sapply(1:nrow(targetMod), function(x) dnorm(x = originDat, mean = targetMod[[x,"Mean"]], 
                                                                      sd = targetMod[[x,"SD"]]) * targetMod[[x,"Weight"]] )
      }
    }
    postDatScaled = data.frame( postDatUnscaled / rowSums(postDatUnscaled) )
    colnames(postDatScaled) = paste0("Target_", thisFeature, 1:ncol(postDatScaled))
    
    ROWNAMES[[thisFeature]] = paste0(thisFeature, 1:nrow(postDatScaled))
    BESTOVERLAP[[thisFeature]] = postDatScaled
  }
  
  # Join entries. Not existing entries will be filled with NAs which we set to 0.
  dtConv = rbindlist(BESTOVERLAP, fill = TRUE)
  dtConv[ is.na(dtConv) ] = 0
  mConv = as.matrix(dtConv) 
  
  rownames( mConv ) = as.vector(unlist(ROWNAMES))
  
  if(plotHeatmap) show(Heatmap( t(mConv), cluster_rows = FALSE, cluster_columns = FALSE, column_title = "Conversion Matrix"))    
  return(mConv)
  
}

liftOverSigs = function(SIGS, convMatrix) {
  
  # We need: Signature by Component
  if(nrow(SIGS) > ncol(SIGS)) SIGS = t(SIGS)
  
  # Match order of columns
  SIGS = SIGS[, match(rownames(convMatrix), colnames(SIGS)) ]
  
  # Convert SIGS and rename components
  liftedSigs = t( SIGS %*% convMatrix )
  
  return(liftedSigs)
}

compareSigs = function(LIFTEDSIGS, QUERYSIGS, NPERM = 1e3, CRIT = "MAX") {
  
  # We need: Component by Signature
  if(nrow(QUERYSIGS) < ncol(QUERYSIGS)) QUERYSIGS = t(QUERYSIGS)
  
  # Match order
  QUERYSIGS = QUERYSIGS[ match(gsub("Target_", "", rownames(LIFTEDSIGS)), rownames(QUERYSIGS)), ]
  
  # Do the cosine similarity comparison by vector
  sigCos = t( apply(LIFTEDSIGS, 2, function(x) cosine(x, QUERYSIGS)) )
  
  # Estimate threshold
  sigCosSim = apply(LIFTEDSIGS, 2, function(x) {
    apply(QUERYSIGS, 2, function(y) {
      if(CRIT == "MAX") {
        max(cosinePerm(x, y, include.perm = TRUE, nperm = NPERM)$estimate.random)  
      } else {
        quantile(cosinePerm(x, y, include.perm = TRUE, nperm = NPERM)$estimate.random, probs = 0.95)  
      }
    })
  })
  
  return(list("sigCos" = sigCos, "sigCosSim" = t(sigCosSim)))
  
}

compareExposures = function(EXP1, EXP2, NPERM = 1e3, CRIT = "MAX") {
  
  # Order and filt samples
  EXP1 = EXP1[ rownames(EXP1) %in% rownames(EXP2), ]
  EXP2 = EXP2[ rownames(EXP2) %in% rownames(EXP1), ]
  EXP2 = EXP2[ match( rownames(EXP1), rownames(EXP2)), ]
  
  # Failsafe: If a TCGA signature has no exposure in OV samples, then add machine precision so the cosine function throws no error
  EXP1[ , colSums(EXP1) == 0 ] = EXP1[ , colSums(EXP1) == 0 ] + .Machine$double.eps
  EXP2[ , colSums(EXP2) == 0 ] = EXP2[ , colSums(EXP2) == 0 ] + .Machine$double.eps
  
  # Do the cosine similarity comparison by vector
  expCos = apply(EXP1, 2, function(x) cosine(x, EXP2))
  
  # Estimate threshold
  expCosSim = apply(EXP1, 2, function(x) {
    apply(EXP2, 2, function(y) {
      if(CRIT == "MAX") {
        max(cosinePerm(x, y, include.perm = TRUE, nperm = NPERM)$estimate.random)  
      } else {
        quantile(cosinePerm(x, y, include.perm = TRUE, nperm = NPERM)$estimate.random, probs = 0.95)  
      }
    })
  })
  
  return(list("expCos" = expCos, "expCosSim" = expCosSim))
  
}

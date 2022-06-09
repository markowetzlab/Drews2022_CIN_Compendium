idSmoothingTargets = function(dfAllSegs, WIGGLE, colNameSegVal, colNameChr, IGNOREDELS = TRUE) {
  
  ### Check column name
  testSegVal = dfAllSegs[[colNameSegVal]][1]
  testChr = dfAllSegs[[colNameChr]][1]
  if(! is.numeric(testSegVal)) { stop("Segment Value column has no numeric value in it. Supplied correct column name? Forgot conversion?")}
  if(is.null(testSegVal)) { stop("Chromosome column has no numeric value in it. Supplied correct column name?")}
  
  # take differences to segment down below
  dfAllSegs$diffs = c( abs( dfAllSegs[[colNameSegVal]][1:(nrow(dfAllSegs)-1)] - dfAllSegs[[colNameSegVal]][2:nrow(dfAllSegs)] ), WIGGLE+1)
  # set TRUE if difference to next segment is smaller than the user supplied cutoff
  dfAllSegs$smooth = dfAllSegs$diffs <= WIGGLE
  # set all segments which are last in a chromosome to FALSE. This also prevents leaking to other samples and cohorts.
  dfAllSegs$smooth[ cumsum( rle(as.character(dfAllSegs[[colNameChr]]))$lengths ) ] = FALSE
  
  # Ignore deletions if wished
  if(IGNOREDELS) { dfAllSegs$smooth[ dfAllSegs[[colNameSegVal]] == 0 ] = FALSE }
  
  return( dfAllSegs )
}


smoothSegments = function(lRaw, CORES, SMOOTHINGFACTOR, colNameMerge, colNameChr, colNameStart, colNameEnd,
                          IGNOREDELS = TRUE, asDf = FALSE) {
  
  ### Check column names
  test = lRaw[[1]]
  testMerge = test[[colNameMerge]][1]
  testChr = test[[colNameChr]][1]
  testStart = test[[colNameStart]][1]
  testEnd = test[[colNameEnd]][1]
  if(! is.numeric(testMerge)) { stop("Merge column has no numeric value in it. Supplied correct column name?")}
  if(is.null(testChr)) { stop("Chromosome column has no numeric value in it. Supplied correct column name?")}
  if(! is.numeric(testStart)) { stop("Start column has no numeric value in it. Supplied correct column name?")}
  if(! is.numeric(testEnd)) { stop("End column has no numeric value in it. Supplied correct column name?")}
  
  # Add diff column to names we want to keep when merging (comes from function "idSmoothingTargets").
  colNameMerge = c(colNameMerge, "diffs")
  
  registerDoMC(CORES)
  lSmooth = foreach(thisSample = lRaw, .final = function(x) setNames(x, names(lRaw)) ) %dopar% {
    
    thisOut = thisSample
    stillSmoothing = sum(thisOut$smooth)
    while( stillSmoothing > 0 ) {
      # For the while loop:
      # Read lines from thisSample and change in thisOut. Hence for a new iteration I need to sync the two.
      thisSample = thisOut
      
      rleRaw = rle(thisSample$smooth)
      # This takes the indeces of the FALSE chains and adds 1. This should give you the next segment which is TRUE.
      # Two challenges:
      # 1) Last segment always FALSE (see above), hence removal of the last number as this would indicate to a segment outside the df.
      # 2) If it starts with a TRUE segment, this would not be found when looking at the FALSE chains. Hence, adding index 1 manually if chain starts with TRUE.
      indRaw = cumsum(rleRaw$lengths)[ ! rleRaw$values ] + 1
      indRaw = indRaw[ -length(indRaw) ]
      if( rleRaw$values[1] ) { indRaw = c(1, indRaw) }
      
      # loop over start indices of TRUE chains.
      for(i in indRaw) {
        # detect length of segments to smooth. add 1 as the last segment has a FALSE value in it but still belongs to this chain.
        endOfStreak = i + rle(thisSample$smooth[i:nrow(thisSample)])$lengths[1]
        # extract reads
        dfMerge = thisSample[i:endOfStreak,]
        
        # too stupid to make this work with data.table
        newElement = as.data.frame( dfMerge[1,] )
        # Get new end and check first wether valid number.
        newEnd = dfMerge[nrow(dfMerge),][[colNameEnd]]
        if(! is.null(newEnd)) {
          newElement[[colNameEnd]] = newEnd
        } else {
          stop("New end coordinate is null. Supplied correct column name?")
        }
        ## Column "segVal" will be dealt with in a minute. Column "diffs" later when running again idSmoothingTargets.
        
        # Merge cn specifically by taking the length of the elements into consideration
        widthWeights = dfMerge[[colNameEnd]] - dfMerge[[colNameStart]]
        newElement[[colNameMerge[1]]] = weighted.mean(dfMerge[[colNameMerge[1]]], widthWeights)
        # Replace all to merge segments with the new merged segment. Later delete duplicated.
        thisOut[i:endOfStreak,] = newElement
      }
      
      # as we have replaced all segments with the new mean segment, we need to remove the duplicates
      thisOut = thisOut[ ! duplicated(thisOut), ]
      # again detect segments which needs smoothing
      thisOut = idSmoothingTargets(thisOut, SMOOTHINGFACTOR, colNameSegVal = colNameMerge[[1]], colNameChr = colNameChr,
                                   IGNOREDELS = IGNOREDELS)
      stillSmoothing = sum(thisOut$smooth)
    }
    
    # after smoothing is finished, change name of cohort
    thisOut$smooth = NULL
    thisOut$diffs = NULL
    return( thisOut )
  }
  
  if( isTRUE(asDf) ) {
    dfSmooth = setDT( rbindlist( lSmooth ) )
    return( dfSmooth )
  } else {
    return( lSmooth )
  }
  
}



## Load and extract CN profiles
calcSigsFromCNProfiles = function(DAT, OUT, OUTFILEADDON, WIGGLE, IGNOREDELS, CORES, REMOVEQUIET, PREPATH, RMNORM, INPUTMODELS,
                                  UNINFPRIOR, SIGNATUREFILE, SAVEALLFILES=FALSE) {
  dfAllSegs = readRDS(DAT)
  colnames(dfAllSegs) == c("chromosome", "start", "end", "segVal", "sample")
  
  
  ## Smoothing normal segments and merging normal and deleted segments
  # Explicit conversion to numeric. Just to be on the safe site.
  dfAllSegs$start = as.numeric( dfAllSegs$start )
  dfAllSegs$end = as.numeric( dfAllSegs$end )
  dfAllSegs$segVal = as.numeric( dfAllSegs$segVal )
  
  if(SAVEALLFILES) {
    write.table(dfAllSegs, file.path(OUT, paste0("0_TCGA_PCAWG_absolute_CN", OUTFILEADDON, ".txt")), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    saveRDS(dfAllSegs, file.path(OUT, paste0("0_TCGA_PCAWG_absolute_CN", OUTFILEADDON, ".rds")) )
  }
  
  # Set everything very close to 2 to 2
  dfAllSegs$segVal[ dfAllSegs$segVal > (2-WIGGLE) & dfAllSegs$segVal < (2+WIGGLE) ] = 2
  # Merge segments only when two normal follow each other -> SMOOTHINGFACTOR = 0
  dfAllSegs = idSmoothingTargets(dfAllSegs, WIGGLE = 0, colNameSegVal = "segVal", colNameChr = "chromosome", IGNOREDELS = IGNOREDELS)
  # Split by sample name
  lRaw = split(dfAllSegs, dfAllSegs$sample)
  
  # Smooth segments by taking the weighted average of the segVal and their lengths
  lSmooth = smoothSegments(lRaw, CORES, SMOOTHINGFACTOR = 0, colNameMerge = "segVal", colNameChr = "chromosome",
                           colNameStart = "start", colNameEnd = "end", IGNOREDELS = IGNOREDELS, asDf = FALSE)
  dtSmooth = rbindlist(lSmooth)
  
  # Write txt and rds
  if(SAVEALLFILES) {
    write.table(dtSmooth, file.path(OUT, paste0("1_TCGA_PCAWG_absolute_CN_smoothed", OUTFILEADDON, ".txt")), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    saveRDS(dtSmooth, file.path(OUT, paste0("1_TCGA_PCAWG_absolute_CN_smoothed", OUTFILEADDON, ".rds")) )
  }
  
  
  ## Avoid measurement errors
  dtSmooth$segVal[ dtSmooth$segVal < 0 ] = 0
  dtSmooth$sample = factor(dtSmooth$sample)
  
  
  
  ## Remove normal samples (should be removed beforehand but for some data sets this is done here)
  if(REMOVEQUIET) {
    # Identify CNAs per sample
    dtCNAs = dtSmooth[ dtSmooth$segVal != 2, ]
    quietSamples = names(table(dtCNAs$sample))[ table(dtCNAs$sample) < 20 ]
    
    # Remove quiet samples
    dtSmooth = dtSmooth[ ! dtSmooth$sample %in% quietSamples, ]
    dtSmooth$sample = factor(dtSmooth$sample)
  }
  
  ## Extract features
  dfBR = data.frame(dtSmooth)
  lBR = split( dfBR, dfBR$sample )
  brECNF = extractCopynumberFeatures( lBR, cores = CORES, prePath=PREPATH, rmNorm = RMNORM )
  
  if(SAVEALLFILES) {
    saveRDS(brECNF, file.path(OUT, paste0("2_TCGA_PCAWG_ECNF", OUTFILEADDON, ".rds")) )
  }
  
  ## Create input matrix
  allModels = readRDS(INPUTMODELS)
  allFeatures = names(allModels)
  
  lMats = lapply(allFeatures, function(thisFeature) {
    
    print(thisFeature)
    thisEcnf = brECNF[[ thisFeature ]]
    thisModel = allModels[[ thisFeature ]]
    
    dat = as.numeric( thisEcnf[,2] )
    # We want a posterior, hence likelihood (density) times prior (weight)
    if( ncol(thisModel) == 2 ) {
      # Poisson model
      print("Poisson-based posterior")
      if(UNINFPRIOR){
        postDatUnscaled = sapply(1:nrow(thisModel), function(x) dpois(x = dat, lambda = thisModel[[x,"Mean"]]) )
      } else {
        postDatUnscaled = sapply(1:nrow(thisModel), function(x) dpois(x = dat, lambda = thisModel[[x,"Mean"]]) * thisModel[[x, "Weight"]] )
      }
      
    } else {
      # Gaussian model
      print("Gauss-based posterior")
      if(UNINFPRIOR){
        postDatUnscaled = sapply(1:nrow(thisModel), function(x) dnorm(x = dat, mean = thisModel[[x,"Mean"]], sd = thisModel[[x,"SD"]]) )
      } else {
        postDatUnscaled = sapply(1:nrow(thisModel), function(x) dnorm(x = dat, mean = thisModel[[x,"Mean"]], sd = thisModel[[x,"SD"]]) * thisModel[[x, "Weight"]] )
      }
    }
    
    # Normalise densities to probabilities
    postDatScaled = data.frame( postDatUnscaled / rowSums(postDatUnscaled) )
    postDatScaled$Sample = thisEcnf[,1]
    matSxC = aggregate(. ~ Sample, postDatScaled, sum)
    rownames(matSxC) = matSxC$Sample
    matSxC$Sample = NULL
    matSxC = as.matrix(matSxC)
    
    # Should be sorted but just to be sure
    matSxC = matSxC[ , order(thisModel[,"Mean"]) ]
    colnames(matSxC) = paste0( thisFeature, 1:ncol(matSxC) )
    
    return(matSxC)
    
  } )
  
  allMats = do.call(cbind, lMats)
  
  if(SAVEALLFILES) {
    saveRDS(object = allMats, file = file.path(OUT, paste0("3_TCGA_PCAWG_SxC_Matrix", OUTFILEADDON, ".rds")))
    write.table(x = allMats, file = file.path(OUT, paste0("3_TCGA_PCAWG_SxC_Matrix", OUTFILEADDON, ".txt")),
                quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
  }
  
  allMats = t( allMats )
  # Definitely should be saved because signatures are derived from those files.
  write.table(x = allMats, file = file.path(OUT, paste0("3_TCGA_PCAWG_CxS_Matrix", OUTFILEADDON, ".txt")),
              quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
  saveRDS(object = allMats, file = file.path(OUT, paste0("3_TCGA_PCAWG_CxS_Matrix", OUTFILEADDON, ".rds")))
  
  
  ## Call exposures
  V = allMats
  if(grepl(pattern = "rds", x = SIGNATUREFILE, ignore.case = TRUE)) {
    W = readRDS(SIGNATUREFILE)
  } else {
    W = read.table(SIGNATUREFILE, header = TRUE, sep = "\t", row.names = 1)
  }
  
  
  # Sanity check signature matrix
  if(nrow(W) < ncol(W)) W = t(W)
  # Check order of components and fix if necessary
  if(! identical(rownames(W), rownames(V)) ) {
    W = W[ match(rownames(V), rownames(W)), ]
  }
  
  ### YAPSA needs:
  ## Full matrix V        mutCatalogue        components (rows) by samples (cols)    <= HAVE
  ## Left matrix W        sigCatalogue        components (rows) by signature (cols)   <= HAVE
  ## Right matrix H       expCatalogue        signature (rows) by samples (cols)     <= WANT
  
  # in_mutation_catalogue_df => NxM => N - Features, M - Samples => Component by Sample matrix
  # in_signatures_df => NxL => N - Features, L - Signatures => Component by Signature matrix
  Hraw = as.matrix( LCD( in_mutation_catalogue_df = V, in_signatures_df = W, in_per_sample_cutoff = 0 ) )
  if(SAVEALLFILES) {
    saveRDS(object = Hraw, file = file.path(OUT, paste0("4_TCGA_PCAWG_Exposures_to_TCGA_Signatures_raw", OUTFILEADDON, ".rds")))
  }
  
  H = t( apply(Hraw, 2, function(x) x/sum(x)) )
  if(SAVEALLFILES) {
    write.table(x = H, file = file.path(OUT, paste0("4_TCGA_PCAWG_Exposures_to_TCGA_Signatures", OUTFILEADDON, ".txt")), quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
  }
  saveRDS(object = H, file = file.path(OUT, paste0("4_TCGA_PCAWG_Exposures_to_TCGA_Signatures", OUTFILEADDON, ".rds")))
  
  return(H)
}

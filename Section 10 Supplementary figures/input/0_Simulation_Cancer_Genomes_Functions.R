## Load and prepare files
loadChromSizes = function(fileChrSizes, rmChr = FALSE) {
  # Load chromosome sizes
  dfChrSizes = read.table(fileChrSizes,sep="\t",stringsAsFactors = F)[1:24,]
  colnames(dfChrSizes) = c("chr", "size")
  # Remove y
  dfChrSizes = dfChrSizes[ dfChrSizes$chr != "chrY", ]
  
  # Remove "chr" prefix if wished
  if(rmChr) {
    # Remove chr and reinstate factors
    dfChrSizes$chr = gsub("chr", "", dfChrSizes$chr)
    chrOrder = c(1:22, "X")
  } else {
    chrOrder = paste0("chr", c(1:22, "X"))
  }
  
  # Order
  dfChrSizes = dfChrSizes[ match(chrOrder, dfChrSizes$chr), ]
  
  return(dfChrSizes)  
}

# Load centromeric regions
loadCentromeres = function(fileCentromeres, rmChr = FALSE) {
  
  dfGaps = read.table(fileCentromeres, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  dfCentromeres = dfGaps[ dfGaps[,8] == "centromere", ]
  dfCentromeres = dfCentromeres[, c(2, 3, 4)]
  colnames(dfCentromeres) = c("Chromosome", "CentStart", "CentEnd")
  
  # Remove "chr" prefix if wished
  if(rmChr) {
    # Remove chr and reinstate factors
    dfCentromeres$Chromosome = gsub("chr", "", dfCentromeres$Chromosome)
    chrOrder = c(1:22, "X")
  } else {
    chrOrder = paste0("chr", c(1:22, "X"))
  }
  
  # Order
  dfCentromeres = dfCentromeres[ match(chrOrder, dfCentromeres$Chromosome), ]
  
  return(dfCentromeres)
}




## Function for sampling from one range with gaps
sampleIntFromGaps = function(n = 1, startVector, endVector, sizeVector) {
  
  # Initialise values
  nRanges = length(sizeVector)
  probRanges = sizeVector/sum(sizeVector)
  
  # First, decide which range to chose
  whichRanges = sample(x = 1:nRanges, size = n, prob = probRanges, replace = TRUE)
  
  # Second, sample from this range
  if(n == 1) {
    vOut = sample(startVector[whichRanges]:endVector[whichRanges], size = 1)
  } else {
    vOut = sapply(1:n, function(i) {
      sample(startVector[whichRanges[i]]:endVector[whichRanges[i]], size = 1)
    })
  }
  
  return(vOut)
}

## Given a chromosome, its centromere, existing CNAs and the size of a new CNA returns
## all possible positions where this CNA could be placed without overlapping.
identifyPositions = function(segsize, THISSIZE, CENTSTART, CENTEND, EXISTCNAS, OVERLAP = FALSE) {
  
  # Output: three vectors: start, end positions of all possible regions and their size
  
  # Create iranges object with sets of three segments:
  #   1) Pre-centromeric regions (that is the region where the new CNA could end in the centromere which is not defined)
  #   2) Centromere (New CNA could start here which is not defined)
  #   3) Existing CNAs
  
  ## Catch small chromosome arms which are smaller than target segment size
  PRECENTSTART = ifelse(CENTSTART-segsize > 0, CENTSTART-segsize, 0)
  PRECENTEND = ifelse(CENTEND-segsize > 0, CENTEND-segsize, CENTSTART)
  
  # Summarise regions where CNAs can or cannot be places
  if(! OVERLAP) {
    # No overlap allowed
    irCovered = IRanges(start = c(PRECENTSTART, CENTSTART, as.numeric(EXISTCNAS$start)-segsize), 
                        end = c(PRECENTEND, CENTEND, as.numeric(EXISTCNAS$end)))
  } else {
    # Overlap allowed
    irCovered = IRanges(start = c(PRECENTSTART, CENTSTART, as.numeric(EXISTCNAS$start)-segsize, as.numeric(EXISTCNAS$end)-segsize),
                        end = c(PRECENTEND, CENTEND, as.numeric(EXISTCNAS$start), as.numeric(EXISTCNAS$end)))
  }
  
  # Edge case: If coordinates are negative (-> 3) )
  start(irCovered)[ start(irCovered) < 0 ] = 0
  
  # Reduce to remove overlaps (also sorts the elements)
  irCovered = reduce(irCovered)
  
  # Then show gaps - that's where CNAs can be placed
  irGaps = gaps(irCovered)
  
  # Add beginning of chromosome (if not already covered)
  if(start(irCovered)[1] != 0) {
    startVector = c(0, start(irGaps))
    endVector = c(start(irCovered)[1]-1, end(irGaps))
  } else {
    startVector = start(irGaps)
    endVector = end(irGaps)
  }
  
  # Add end of chromosome
  lastElemEnd = end(irCovered)[length(end(irCovered))]
  realEnd = THISSIZE - segsize
  # If end of chromosome is not covered by an element, then add it manually
  if(lastElemEnd < realEnd) {
    # Add last segment manually
    startVector = append(startVector, lastElemEnd+1)
    endVector = append(endVector, realEnd)
  }
  
  # Get sizes for easier sampling later
  sizePossible = endVector - startVector
  
  return(list("start" = startVector, "end" = endVector, "size" = sizePossible))
}


# Fill chromosome with normal segments between CNAs created by mutational processes
fillChromWithNormals = function(DTCHROM, THISCHROM, THISSIZE, NORMALCP, SIGNATURE = "BG", OUTGR = TRUE) {
  
  # Make normal and aberrant granges objects
  grNormals = GRanges(THISCHROM, IRanges(start = 0, end = THISSIZE))
  grCNAs = makeGRangesFromDataFrame(DTCHROM, keep.extra.columns=TRUE)
  
  # Merge and thereby fill gaps
  grGaps = setdiff(grNormals, grCNAs)
  
  # Return if chromosome is already filled (no gaps)
  if(length(grGaps) == 0) { return(DTCHROM) }
  
  # Define background
  grGaps$segVal = NORMALCP
  grGaps$Signature = SIGNATURE
  
  # Merge and order segments
  grComplete = c(grGaps, grCNAs)
  grComplete = grComplete[ order(grComplete) ]
  names(grComplete) = NULL
  
  if(OUTGR) {
    return(grComplete)
  } else {
    # Convert back to data frame and return
    dfCNAs = as.data.frame(grComplete)
    dfCNAs$width = NULL
    dfCNAs$strand = NULL
    colnames(dfCNAs) = colnames(DTCHROM)
    
    return(dfCNAs)  
  }
}

resolveOverlaps = function(dfCNAs) {
  
  ## Identify overlaps
  # Create genomic ranges object
  dfCNAs$start = as.numeric(dfCNAs$start)
  dfCNAs$end = as.numeric(dfCNAs$end)
  grOverlap = makeGRangesFromDataFrame(dfCNAs, keep.extra.columns = TRUE)
  
  # Check if all CNAs are not overlapping
  if(! isDisjoint(grOverlap)) {
    # All overlapping
    grNotOv = disjoin(grOverlap, with.revmap = TRUE)
    
    # Loop over all entries and decide which CN state should be put in
    grNotOv$segVal = NA
    for(ind in 1:length(grNotOv$revmap)) {
      
      # grNotOv - disjoint / not overlapping - new
      newInd = ind
      # grOverlap - not-disjoint / overlapping - old
      oldInd = grNotOv$revmap[[ ind ]]
      
      if(length(oldInd) == 1) {
        # One to one mapping from overlapping and non-overlapping Genomic Ranges overlap
        grNotOv$segVal[ newInd ] = grOverlap$segVal[ oldInd ]
      } else {
        # Overlapping element with multiple indeces
        # Since overlaps are resulted for each CNA and the CNAs are sorted
        # this leads to the fact that the old index is always the second or last entry
        grNotOv$segVal[ newInd ] = grOverlap$segVal[ oldInd[length(oldInd)] ]
      }
    }
    
    # Just to check if still disjoint afterwards
    if(! isDisjoint(grNotOv)) { print("[Warning]: CNAs still disjoint despite resvoling mechanism.") }
  } else {
    # Already disjointed
    return(dfCNAs)
  }
  
  # Convert grNotOv to dfCNAs and return
  dfCNARes = as.data.frame(grNotOv)
  dfCNARes$strand = NULL
  dfCNARes$revmap = NULL
  dfCNARes$width = NULL
  colnames(dfCNARes) = colnames(dfCNAs)
  
  return(dfCNARes)
}


# For simulating HRD (via LSTs), focal amplifications (via ecDNAs) and tolerance of WGD (via losses)
simulateMutProcess = function(THISCHROM, THISSIZE, LISTVARS, THISSIG, DTCHROM = NULL) {
  
  # Initiate process by determining number of CNAs
  nCNAs = rpois(1, LISTVARS$NUMCNASMEAN)
                
  # Failsafe to catch zeros
  if(nCNAs == 0) { nCNAs = 1 }
  
  # Calculate size of centromere and practical size of chromosomes
  centromSize = LISTVARS$CENTEND - LISTVARS$CENTSTART
  realSize = THISSIZE - centromSize
  
  # Check number of chromosomes and largest segment size work 
  if(nCNAs * LISTVARS$SEGSIZEHIGH > realSize) {
    # Case number 1: Just reduce number of CNAs
    if(LISTVARS$SEGSIZEHIGH < realSize) {
      # Is at least 1
      nCNAs = floor(realSize / LISTVARS$SEGSIZEHIGH)
    }
    # Case number 2: Segment size larger than actual chromosome size
    if(LISTVARS$SEGSIZEHIGH > realSize) {
      LISTVARS$SEGSIZEHIGH = 2/3*realSize
      # Adjust lowest end if we now created a segment size conundrum
      if(LISTVARS$SEGSIZEHIGH < LISTVARS$SEGSIZELOW){
        LISTVARS$SEGSIZELOW = 2/3*LISTVARS$SEGSIZEHIGH
      }
      nCNAs = 1
    }
  }
  
  # Communicate if wished
  if(LISTVARS$VERBOSE) {
    print(paste("uuid:", LISTVARS$uuid))
    print(paste("nCNAs:", nCNAs))
    print(paste("thisChrom:", THISCHROM))
    print(paste("thisSize:", THISSIZE))
    print(paste("realSize:", realSize))
    print(paste("thisSig:", THISSIG))
  }
  
  if(is.null(DTCHROM)) {
    # Initiate number of CNAs for this sample (no sample for now to keep it numeric)
    dfCNAs = setNames(data.frame(matrix(ncol = 4, nrow = 0)), 
                      c("chromosome", "start", "end", "segVal"))
  } else {
    # Already data frame submitted via parameter. Transfer to variable.
    dfCNAs = DTCHROM
    
    # Convert to numeric to be sure
    dfCNAs$start = as.numeric(dfCNAs$start)
    dfCNAs$end = as.numeric(dfCNAs$end)
    dfCNAs$segVal = as.numeric(dfCNAs$segVal)
  }
  
  ## Sample segment size and changepoints for the CNAS
  allSegsizes = round(runif(nCNAs, LISTVARS$SEGSIZELOW, LISTVARS$SEGSIZEHIGH))
  # Start with the largest CNA
  allSegsizes = sort(allSegsizes, decreasing = TRUE)
  # Draw change points to CNAs
  allCpsRaw = signif(runif(nCNAs, LISTVARS$CPLOW, LISTVARS$CPHIGH), 6)
  
  # If wished turn sign of a few CPs (mostly for LST)
  if(! is.null(LISTVARS$PROBSIGNPOS)) {
    # Identify a few CNAs randomly to switch sign
    whichTurn = as.logical(rbinom(length(allCpsRaw), size = 1, prob = LISTVARS$PROBSIGNPOS))
    allCpsRaw[ whichTurn ] = -allCpsRaw[ whichTurn ]
  }
  
  # Add background copy number to create final CNA copy number
  allCps = LISTVARS$NORMALCP + allCpsRaw
  
  # Loop over each CNA and first identify all possible locations where it can fit, 
  # then sample a new position from all possible locations
  for(i in 1:nCNAs) {
    
    if(LISTVARS$VERBOSE) { print(i) }
    
    ## Function to identify possible positions to fit a segment of a given size
    lPositions = identifyPositions(allSegsizes[i], THISSIZE, LISTVARS$CENTSTART, LISTVARS$CENTEND, dfCNAs, OVERLAP = LISTVARS$OVERLAP)
    
    # Check if enough space is there for new sample
    if(length(lPositions$start) == 0 || length(lPositions$end) == 0 || length(lPositions$size) == 0) {
      if(LISTVARS$VERBOSE) {print(paste("Skip", i))}
      next
    }
    
    # Sample from all possible ranges (sampling from range with gaps)
    thisPos = sampleIntFromGaps(n = 1, lPositions$start, lPositions$end, lPositions$size)
    
    # Add to final data frame
    dfCNAs[nrow(dfCNAs)+1,] = c(THISCHROM, thisPos, thisPos+allSegsizes[i], allCps[i])
    
    # Sort CNAs
    dfCNAs = dfCNAs[order(as.numeric(dfCNAs$start)),]
    
    # If overlaps are allowed, check for overlapping CNA regions and if so resolve them
    if(nrow(dfCNAs) > 1) {
      overlaps = sum(as.numeric(dfCNAs$end[1:(nrow(dfCNAs)-1)]) > as.numeric(dfCNAs$start[2:nrow(dfCNAs)]))
      if(LISTVARS$OVERLAP && overlaps > 0) { dfCNAs = resolveOverlaps(dfCNAs) }
    }
  }
  
  return(dfCNAs)
}

# For simulating mitotic errors
simulateChrSegregation = function(THISCHROM, THISSIZE, LISTVARS, DTCHROM = NULL)  {
  
  # Create dfCNA
  if(is.null(DTCHROM)) {
    # Initiate number of CNAs for this sample (no sample for now to keep it numeric)
    dfCNAs = setNames(data.frame(matrix(ncol = 4, nrow = 0)), 
                      c("chromosome", "start", "end", "segVal"))
  } else {
    # Already data frame submitted via parameter. Transfer to variable.
    dfCNAs = DTCHROM
    
    # Convert to numeric to be sure
    dfCNAs$start = as.numeric(dfCNAs$start)
    dfCNAs$end = as.numeric(dfCNAs$end)
    dfCNAs$segVal = as.numeric(dfCNAs$segVal)
  }
  
  # Decide if chromosome undergoes missegregation and then whether its lost or gained
  if(as.logical(rbinom(1, 1, LISTVARS$PROB))) {
    
    # Decide if lost or gained
    if(as.logical(rbinom(1, 1, 0.5))) {
      # A bit narrower than usually because whole-chromosome changes are mostly captured at integer states
      thisCP = LISTVARS$NORMALCP - signif(runif(1, LISTVARS$CPLOW, LISTVARS$CPHIGH), 6)
    } else {
      thisCP = LISTVARS$NORMALCP + signif(runif(1, LISTVARS$CPLOW, LISTVARS$CPHIGH), 6)
    }
    
  } else {
    # No chromosome missegregation is taking place. Assign background copy number.
    thisCP = LISTVARS$NORMALCP  
  }
  
  # Create data frame with CNA(s) depending on whether previous CNAs exist
  if(is.null(DTCHROM) || nrow(dfCNAs) == 0) {
    # Assemble CNA and return - simple case: just one entry
    dfCNAs[nrow(dfCNAs)+1,] = c(THISCHROM, 1, THISSIZE, thisCP)
  } else {
    # More complex case with pre-existing CNAs that need to be amended and the space between them filled
    # Adjust copy number of already existing CNAS
    dfCNAs$segVal = dfCNAs$segVal + ( thisCP - LISTVARS$NORMALCP )
    thisSig = dfCNAs$Signature[1]
    dfCNAs = fillChromWithNormals(DTCHROM = dfCNAs, THISCHROM, THISSIZE, NORMALCP = thisCP, SIGNATURE = thisSig, OUTGR = FALSE)
    
  }
  
  return(dfCNAs)  
}


rescueHomdels = function(DTCHROM, RESCUECNLOW = 0.51, RESCUECNHIGH = 0.55, MAXSIZE = 1e7) {
  
  # Return if no present
  if(sum(DTCHROM$segVal < 0) == 0) { return(DTCHROM) }
  
  # Failsafe: Check that all columns are numeric
  DTCHROM$start = as.numeric(DTCHROM$start)
  DTCHROM$end = as.numeric(DTCHROM$end)
  DTCHROM$segVal = as.numeric(DTCHROM$segVal)
  
  # All homdels with size smaller MAXSIZE set to zero
  DTCHROM$width = DTCHROM$end - DTCHROM$start
  DTCHROM$segVal[ DTCHROM$segVal < 0 & DTCHROM$width < MAXSIZE ] = 0
  
  # Rescue homdels longer than MAXSIZE
  numRescue = nrow(DTCHROM[ DTCHROM$segVal < 0 & DTCHROM$width >= MAXSIZE, ])
  DTCHROM$segVal[ DTCHROM$segVal < 0 & DTCHROM$width >= MAXSIZE ] = signif(runif(numRescue, RESCUECNLOW, RESCUECNHIGH), 6)
  
  DTCHROM$width = NULL
  
  return(DTCHROM)
}



callSignature = function(THISCHROM, THISSIZE, THISSIG, LISTVARS, DTCHROM = NULL) {

  if(THISSIG == "LST") { 
    
    ## Set the right variables for this mutational process
    if(LISTVARS$VERBOSE) { print("LST") }
    
    if(is.null(LISTVARS[["NUMCNASMEAN"]])) { LISTVARS[["NUMCNASMEAN"]] = 7 }
    if(is.null(LISTVARS[["SEGSIZELOW"]])) { LISTVARS[["SEGSIZELOW"]] = 5e6 }
    if(is.null(LISTVARS[["SEGSIZEHIGH"]])) { LISTVARS[["SEGSIZEHIGH"]] = 3e7 }
    if(is.null(LISTVARS[["CPLOW"]])) { LISTVARS[["CPLOW"]] = -1.25 }
    if(is.null(LISTVARS[["CPHIGH"]])) { LISTVARS[["CPHIGH"]] = -0.85 }
    if(is.null(LISTVARS[["PROBSIGNPOS"]])) { LISTVARS[["PROBSIGNPOS"]] = 0.25 }
    
    
    ## Call the mutational process with or without previous segments
    if(is.null(DTCHROM)) {
      dtChrom = simulateMutProcess(THISCHROM, THISSIZE, LISTVARS, THISSIG)   
    } else {
      dtChrom = simulateMutProcess(THISCHROM, THISSIZE, LISTVARS, THISSIG, DTCHROM = DTCHROM) 
    }
    
  }
  
  if(THISSIG == "ECDNA") { 
    
    ## Set the right variables for this mutational process
    if(LISTVARS$VERBOSE) { print("ECDNA") }
    
    if(is.null(LISTVARS[["NUMCNASMEAN"]])) { LISTVARS[["NUMCNASMEAN"]] = 3 }
    if(is.null(LISTVARS[["SEGSIZELOW"]])) { LISTVARS[["SEGSIZELOW"]] = 1e6 }
    if(is.null(LISTVARS[["SEGSIZEHIGH"]])) { LISTVARS[["SEGSIZEHIGH"]] = 3e6 }
    if(is.null(LISTVARS[["CPLOW"]])) { LISTVARS[["CPLOW"]] = 7 }
    if(is.null(LISTVARS[["CPHIGH"]])) { LISTVARS[["CPHIGH"]] = 50 }
    
    
    ## Call the mutational process with or without previous segments
    if(is.null(DTCHROM)) {
      dtChrom = simulateMutProcess(THISCHROM, THISSIZE, LISTVARS, THISSIG)   
    } else {
      dtChrom = simulateMutProcess(THISCHROM, THISSIZE, LISTVARS, THISSIG, DTCHROM = DTCHROM) 
    }
    
  }
  
  if(THISSIG == "WGD") { 
    
    ## Set the right variables for this mutational process
    if(LISTVARS$VERBOSE) { print("WGD") }
    
    dtChrom = setNames(data.frame(matrix(ncol = 4, nrow = 0)), 
                       c("chromosome", "start", "end", "segVal"))
    dtChrom[nrow(dtChrom)+1,] = c(THISCHROM, 1, THISSIZE, LISTVARS$NORMALCP)
    
  }
  
  
  if(THISSIG == "CHR") { 
    
    ## Set the right variables for this mutational process
    if(is.null(LISTVARS[["CPLOW"]])) { LISTVARS[["CPLOW"]] = 0.9 }
    if(is.null(LISTVARS[["CPHIGH"]])) { LISTVARS[["CPHIGH"]] = 1.1 }
    if(is.null(LISTVARS[["PROB"]])) { LISTVARS[["PROB"]] = 0.2 }
    
    
    ## Call the mutational process with or without previous segments
    if(is.null(DTCHROM)) {
      dtChrom = simulateChrSegregation(THISCHROM, THISSIZE, LISTVARS)
    } else {
      dtChrom = simulateChrSegregation(THISCHROM, THISSIZE, LISTVARS, DTCHROM = DTCHROM)
    }
    
  }
  
  if(THISSIG == "PASS") {
    
    ## Passthrough diploid chromosome (for other copy number states see "WGD" function)
    dtChrom = setNames(data.frame(matrix(ncol = 4, nrow = 0)), 
                      c("chromosome", "start", "end", "segVal"))
    dtChrom[nrow(dtChrom)+1,] = c(THISCHROM, 1, THISSIZE, 2)
    
  }
  
  return(dtChrom)
}


simulateAGenome = function(WHICHSIG, CHRSIZES, CENTSIZES, NORMALCP = 2, OUTGR = TRUE, LISTVARS = NULL, VERBOSE = FALSE) {
  
  # Initiate sample by assigning an ID
  thisUUID = stri_rand_strings(1, 6, pattern = "[A-Z0-9]")
  
  # Initiate list of variables for simulation functions
  if(is.null(LISTVARS)) {
    LISTVARS = list()
  }
  LISTVARS[["uuid"]] = thisUUID
  LISTVARS[["VERBOSE"]] = VERBOSE
  LISTVARS[["NORMALCP"]] = NORMALCP
  
  # Check if overlaps of CNAs are allowed or not
  if(is.null(LISTVARS[["OVERLAP"]])) { LISTVARS[["OVERLAP"]] = FALSE }
  
  ## If a '|' is present then assume that WHICHSIG has been produced by sampleMutActivity function and therefore is valid
  multiSig = ifelse(sum(grepl('\\|', WHICHSIG)) == 0, FALSE, TRUE)
  # Checking WHICHSIG if no multiple signatures are present and therefore cannot be assumed to be produced by sampleMutActivity function
  if(! multiSig) {
    # Failsafe for when WHICHSIG is shorter than number of chromosomes
    if(length(WHICHSIG) == 1) {
      WHICHSIG = rep(WHICHSIG[1], nrow(CHRSIZES))
    }
    if(length(WHICHSIG) < nrow(CHRSIZES)) {
      print("[Warning]: Length of WHICHSIG is smaller than number of chromosomes but larger than 1. Take first element and repeat it.")
      WHICHSIG = rep(WHICHSIG[1], nrow(CHRSIZES))
    }
  }
  
  # Identify and loop over chromosomes
  allChroms = unique(CHRSIZES[,1])
  lAllChroms = lapply(1:length(allChroms), function(i) {

    thisChrom = allChroms[i]
    thisSize = CHRSIZES[,2][ CHRSIZES[,1] == thisChrom ]
    thisSig = WHICHSIG[i]
    
    # Add variables for centromeric regions
    LISTVARS[["CENTSTART"]] = CENTSIZES$CentStart[ CENTSIZES$Chromosome == thisChrom ]
    LISTVARS[["CENTEND"]] = CENTSIZES$CentEnd[ CENTSIZES$Chromosome == thisChrom ]
    
    # Checking whether this chromosome will be affected by multiple signatures,
    # If not, then simple production of CNAs
    thisMultiSig = ifelse(sum(grepl('\\|', thisSig)) == 0, FALSE, TRUE)
    if(! thisMultiSig) {
      # Single signature used for creating CNAs on a chromosome
      # Call signature creation function 
      dtChrom = callSignature(thisChrom, thisSize, thisSig, LISTVARS)
      dtChrom$Signature = thisSig
    } else {
      # Multiple signatures present - create background chrom first, then loop over each signature
      theseSigs = strsplit(thisSig, "\\|")[[1]]
      
      # If CHR present, move it to the end
      # Caveat: A signature can only be present once, so just remove entry and add it to the end
      if("CHR" %in% theseSigs) {
        theseSigs = c(theseSigs[ theseSigs != "CHR" ], "CHR")
      }
      
      # Create first dtChrom (just empty data frame)
      dtChrom = setNames(data.frame(matrix(ncol = 4, nrow = 0)), 
                         c("chromosome", "start", "end", "segVal"))
      
      # Then call repeatedly callSignature, every time with the output of the previous call
      for(thatSig in theseSigs) {
        dtChrom = callSignature(thisChrom, thisSize, thatSig, LISTVARS, DTCHROM = dtChrom)
      }
      
      # Failsafe: Ensure copy number states are non-negative
      # With multiple mutational processes it can happen that segVals become negative
      # CNAs smaller than 10MB set to 0, larger than that set to around 0.5
      dtChrom = rescueHomdels(DTCHROM = dtChrom, RESCUECNLOW = 0.51, RESCUECNHIGH = 0.55, MAXSIZE = 1e7)
      
      # Mention mutational processes
      dtChrom$Signature = thisSig
    }
    
    ## Fill up chromosome with normal segments
    dtChrom = fillChromWithNormals(DTCHROM = dtChrom, THISCHROM = thisChrom, THISSIZE = thisSize, NORMALCP = LISTVARS$NORMALCP, OUTGR = FALSE)
    
    return(dtChrom)
  })
  
  # Collate results from lists
  dfCNAs = do.call(rbind, lAllChroms)
  
  # Add sample
  dfCNAs$sample = thisUUID
  
  # Confirm correct order of columns
  correctOrder = c("chromosome", "start", "end", "segVal", "sample", "Signature")
  dfCNAs = dfCNAs[ , correctOrder]
  
  # Return object
  if(OUTGR) {
    grComplete = makeGRangesFromDataFrame(dfCNAs, keep.extra.columns=TRUE)
    return(grComplete)
  } else {
    return(dfCNAs)  
  }
}

## Looping function for generating multiple samples
simulateMultipleGenomes = function(N, WHICHSIG, CHRSIZES, CENTSIZES, AFFECTEDCHR = NULL, SEED = 1, NORMALCP = 2, OUTGR = TRUE, 
                                   LISTVARS = NULL, VERBOSE = FALSE, PROGRESS = TRUE, PNUM = 10) {
  
  set.seed(SEED)
  # If WHICHSIG is a vector (class character)
  if(class(WHICHSIG) == "character") {
    
    lLST = lapply(1:N, function(i) {
      # Report progress
      if(PROGRESS) {
        if((i %% PNUM) == 0) { print(i) }  
      }
      
      # If number of affected chromosomes should be randomly determined
      if(! is.null(AFFECTEDCHR)) {
        
        # Determine number of affected chromosomes while avoiding quiet genomes
        NUMCHR = 0
        while(NUMCHR == 0) { NUMCHR = rpois(1, AFFECTEDCHR) }
        
        # Determine which chromosomes should be affected
        NEWSIG = rep("PASS", nrow(CHRSIZES))
        NEWSIG[ sample(1:nrow(CHRSIZES), NUMCHR, replace = FALSE) ] = WHICHSIG[1]
        WHICHSIG = NEWSIG
      }
      
      return( simulateAGenome(WHICHSIG, CHRSIZES, CENTSIZES, NORMALCP = NORMALCP, 
                              OUTGR = OUTGR, LISTVARS = LISTVARS, VERBOSE = VERBOSE) )
    })
    
  } else {
  
    # If WHICHSIG is a list then ignore N
    lLST = lapply(1:length(WHICHSIG), function(i) {
      # Report progress
      if(PROGRESS) {
        if((i %% PNUM) == 0) { print(i) }  
      }
      
      THESESIGS = WHICHSIG[[i]]
      return( simulateAGenome(WHICHSIG = THESESIGS, CHRSIZES, CENTSIZES, NORMALCP = NORMALCP, 
                              OUTGR = OUTGR, LISTVARS = LISTVARS, VERBOSE = VERBOSE) )
    })
  }
  
  # Identify names and rename list items
  names(lLST) = sapply(lLST, function(x) unique(x$sample))
  
  return(lLST)
}


## General handling of lists of GenomicRanges objects
convertGRListToDF = function(GRLIST) {
  
  lConvert = lapply(GRLIST, function(thisGR) {
    
    # Convert
    dfOut = as.data.frame(thisGR)
    dfOut$width = NULL
    dfOut$strand = NULL
    
    # Check all necessary columns are numerical
    dfOut$start = as.numeric(dfOut$start)
    dfOut$end = as.numeric(dfOut$end)
    dfOut$segVal = as.numeric(dfOut$segVal)
    
    # Rename to be in line with existing workflows
    colnames(dfOut) = c("chromosome", "start", "end", "segVal", "sample", "signature")
    
    return(dfOut)  
  })
  
  return(lConvert)
}


## Plotting functions
convertForPlot = function(thisSample, dfChrSizes) {
  
  # Prepare vector of chromosome sizes (ordering implicitly required)
  vChrSizes = dfChrSizes[,2]
  names(vChrSizes) = dfChrSizes[,1]
  toAdd = c(0, head(cumsum(as.numeric(vChrSizes)), -1))
  names(toAdd) = names(vChrSizes)
  
  thisSample$WGstart = thisSample$start + toAdd[ match(thisSample$chromosome, names(toAdd)) ]
  thisSample$WGend = thisSample$end + toAdd[ match(thisSample$chromosome, names(toAdd)) ]
  
  return(thisSample)
}


chrForPlot = function(dfChrSizes) { 
  vChrSizes = dfChrSizes[,2]
  names(vChrSizes) = dfChrSizes[,1]
  
  vlines = cumsum(as.numeric(vChrSizes))
  textpos = as.numeric(vChrSizes)/2 + c(0, head(cumsum(as.numeric(vChrSizes)), -1))
  
  dfChrDetails = data.frame(chrNames=names(vChrSizes), vlines=vlines, textpos=textpos)
  
  return(dfChrDetails)
}


plotCNProfile = function(dfWG, chrPlotDetails, maxVal = 5, smallPlot = TRUE, name) {
  
  # Prepare arrows
  OUTLIERS=FALSE
  if( sum(dfWG$segVal>maxVal) > 0 ) {
    outliers = dfWG[ dfWG$segVal>maxVal,]
    outliers$lineend = c('butt')
    outliers$linejoin = c('mitre')
    OUTLIERS=TRUE
  }
  
  # Main plot
  p = ggplot(dfWG, aes(x=WGstart))+ ylim(0,5) + 
    theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
    geom_segment(xend=dfWG$WGend, y=dfWG$segVal, yend=dfWG$segVal, size = 1.5 )
  
  # Now add chr lines, change styles and names
  p = p + geom_vline(xintercept = chrPlotDetails$vlines, colour = "gray80", linetype = 2) +
    ylab("Copy number") + ggtitle(name)
  
  # Add red arrows indicating outlier segments
  if(isTRUE(OUTLIERS)) {
    p = p + geom_segment(data = outliers, aes(y = maxVal*0.95, yend = maxVal, xend = WGstart), colour = "red",
                         lineend = outliers$lineend, linejoin = outliers$linejoin, arrow = arrow(length = unit(0.1,"cm")) )
  }
  
  # If TRUE, only every second chr label will be printed
  if(isTRUE(smallPlot)) {
    p = p + geom_text(data = chrPlotDetails[c(TRUE, FALSE),], aes(y=maxVal, x=textpos, label=chrNames), 
                      size=3,vjust=-0.75, hjust=0.5, colour = "gray60")
  } else {
    p = p + geom_text(data = chrPlotDetails, aes(y=maxVal, x=textpos, label=chrNames), 
                      size=3,vjust=-0.75, hjust=0.5, colour = "gray60")
  }
  
  return(p)
}

plotAllSamples = function(lDFSim21, dfChrSizes) {
  
  lAllPlots = lapply(lDFSim21, function(thisSample) {
    
    # Add corrected coordinates for plotting
    thisConv = convertForPlot(thisSample, dfChrSizes)
    
    # Get coordinates for chromosome labels and vertical lines
    dfChrDetails = chrForPlot(dfChrSizes)
    
    # Plot
    sigsUsed = names(table(thisConv$Signature))
    thatName = thisConv$sample[1]
    pWG = plotCNProfile(thisConv, dfChrDetails, name = paste(sigsUsed, "-", thatName))
    
    return(pWG)
  })
}

convAndMelt = function(lLST, NAME) {
  lLSTdf = convertGRListToDF(lLST)
  lstECNF = extractCopynumberFeatures(lLSTdf, cores = 6, prePath=prePath, rmNorm = TRUE)
  lstECNF$copynumber = NULL
  dfLSTECNF = melt(lstECNF)
  dfLSTECNF$sig = NAME
  return(dfLSTECNF)
}


plotFeaturesPerSig = function(dfPlot, SAVE = TRUE) {
  
  # Segment size
  pSeg = ggplot(dfPlot[dfPlot$L1 == "segsize",], aes(x = value, colour = sig)) + geom_density() + 
    scale_x_log10() + labs(x = "Segment size", colour = "Signature") + theme(legend.position = "none") +
    coord_capped_cart(left = "both", bottom = "both") +
    scale_colour_brewer(palette = "Dark2")
  
  # Changepoint
  pChange = ggplot(dfPlot[dfPlot$L1 == "changepoint",], aes(x = value, colour = sig)) + geom_density() +
    labs(x = "Changepoint", colour = "Signature") + theme(legend.position = "none") +
    coord_capped_cart(left = "both", bottom = "both") +
    scale_colour_brewer(palette = "Dark2")
  
  # Osc. chains
  pOsc = ggplot(dfPlot[dfPlot$L1 == "osCN",], aes(x = value, fill = sig)) + geom_histogram(bins = 10) + 
    labs(x = "Osc. chains", fill = "Signature") +
    coord_capped_cart(left = "both", bottom = "both") +
    scale_fill_brewer(palette = "Dark2")
  
  # bp10MB
  pbp10 = ggplot(dfPlot[dfPlot$L1 == "bp10MB",], aes(x = value, fill = sig)) + geom_histogram(bins = 20) + 
    labs(x = "BP 10MB", fill = "Signature") + theme(legend.position = "none") +
    coord_capped_cart(left = "both", bottom = "both") +
    scale_fill_brewer(palette = "Dark2")
  
  # bp chr arm
  pbparm = ggplot(dfPlot[dfPlot$L1 == "bpchrarm",], aes(x = value, fill = sig)) + geom_histogram(bins = 30) + 
    labs(x = "BP per chromosome arm", fill = "Signature") + theme(legend.position = "none") +
    coord_capped_cart(left = "both", bottom = "both") +
    scale_fill_brewer(palette = "Dark2")
  
  
  pFeatDist = pSeg / pChange / pOsc / pbp10 / pbparm
  
  if(SAVE) {
    ggsave(file.path(thisPath, "Out_3_Step_1_ECNF_LST_ECDNA_WGD_CHR_100each.png"), pFeatDist, 
           device = "png", width = 12, height = 12, units = "cm", bg = "white")
    
    ggsave(file.path(thisPath, "Out_3_Step_1_ECNF_LST_ECDNA_WGD_CHR_100each.svg"), pFeatDist, 
           device = "png", width = 12, height = 12, units = "cm", bg = "white")
  }
  
  return(pFeatDist)
}


calcInputMatrix = function(brECNF, INPUTMODELS, UNINFPRIOR, NAME, SAVE = TRUE, THISPATH) {
  
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
  
  # Save output
  if(SAVE) {
    saveRDS(brECNF, file.path(THISPATH, paste0(NAME, ".rds")))
    write.table(x = allMats, file = file.path(THISPATH, paste0(NAME, ".txt")), 
                quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
    
    # Transpose and save matrix again
    allMats = t( allMats )
    newName = gsub("SxC", "CxS", NAME)
    write.table(x = allMats, file = file.path(THISPATH, paste0(newName, ".txt")),
                quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
  }
  
  return(t(allMats))
}



sampleMutActivity = function(N = 1, SIGS = c("LST", "ECDNA", "CHR"), VALS = c(9, 1, 2), REST = "WGD", dfChrSizes, ALLOWMULTI = FALSE, RESCUE = TRUE) {
  
  # Order by frequency - increasing
  newOrder = order(VALS)
  VALS = VALS[newOrder]
  SIGS = SIGS[newOrder]
  
  lMuts = lapply(1:N, function(i) {
    
    # Prepare chromosome vector
    numChr = nrow(dfChrSizes)
    vMuts = rep(REST, numChr)
    
    # Loop over mut processes and place them across the chromosomes
    # More dominant process can overwrite previous processes
    for(i in 1:length(SIGS)) {
      
      thisMut = SIGS[i]
      thisFreq = VALS[i]
      
      # Catch NA and replace with standard values
      if(is.na(thisFreq)) {
        thisFreq = 2
      }
      # print(thisMut)
      
      ## Determine which chromosomes should be affected
      numAffect = rpois(1, thisFreq)
      # Catch edge cases
      if(numAffect == 0) { numAffect = 1}
      if(numAffect >= numChr) { numAffect = numChr - 1 }
      
      affectedChrs = sample(1:numChr, numAffect)
      
      if(! ALLOWMULTI) {
        vMuts[ affectedChrs ] = thisMut  
      } else {
        ## Remove REST because we now replace it with a mutational process
        vMuts[ affectedChrs ] = gsub(REST, "", vMuts[ affectedChrs ])
        
        ## Add with '|'
        vMuts[ affectedChrs ] = sapply(vMuts[ affectedChrs ], function(thisChrom) paste(thisChrom, thisMut, sep = "|"))
      }
      
      ## Remove first character if it is a '|'
      vMuts = gsub("^\\|", "", vMuts)
    }
    
    # If you want to have every mutational process be present at least once then rescue those that are currently not present
    # Only necessary if ALLOWMULTI is FALSE.
    if(RESCUE && ! ALLOWMULTI) {
      # Rescue mut processes => Every supplied process will be present at least once
      notPresent = SIGS[ ! SIGS %in% names(table(vMuts)) ]
      # Very rarely this step overwrites another mut process that is only present once
      while(length(notPresent) > 0 ){
        vMuts[ sample(1:numChr, length(notPresent)) ] = notPresent
        notPresent = SIGS[ ! SIGS %in% names(table(vMuts)) ]
      }
    }
    
    return(vMuts)
  })
  
  if(N==1) { return(lMuts[[1]]) }
  return(lMuts)
}

plotInputMatrix = function(allMats, lDFSim21, NAME, THISPATH, SAVE = TRUE) {
  
  ## Plot input matrix
  dfAllMats = melt(allMats)
  colnames(dfAllMats) = c("Sample", "Component", "Posterior")
  orderForLater = names(lDFSim21)
  dfAllMats$Sample = factor(dfAllMats$Sample, levels = orderForLater)
  
  pInput = ggplot(dfAllMats, aes(x = Component, y = Sample, fill = Posterior)) +
    geom_raster() + scale_fill_viridis_c() + labs(x = "Feature component", y = "Samples") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom") +
    coord_capped_cart(left = "both", bottom = "both")
  
  if(SAVE) {
    ggsave(file.path(THISPATH, paste0(NAME, ".png")), pInput, 
           device = "png", width = 12, height = 12, units = "cm", bg = "white")
    
    ggsave(file.path(THISPATH, paste0(NAME, ".svg")), pInput, 
           device = "png", width = 12, height = 12, units = "cm", bg = "white")
  }
  
  return(pInput)
}

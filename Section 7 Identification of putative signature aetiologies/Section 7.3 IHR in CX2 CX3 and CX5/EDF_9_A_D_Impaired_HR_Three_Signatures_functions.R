#### Functions for extended data figure xx comparing the three signatures to HRD scores

prepareWangDataset = function(WANG, cnState, sampleNames) {
  
  require(data.table)
  
  # Wang et al.
  wang = fread(WANG)
  wangBrca = wang[ wang$BRCA.status != "WT", ]
  wangBrca = wangBrca[ wangBrca$Case_ID %in% sampleNames, ]
  
  # Check LOH status for somatic genes
  wangBrca$origin = "Somatic"
  wangBrca$origin[ grepl("Germline", wangBrca$BRCA.status) ] = "Germline"
  wangBrca$Gene = sapply(strsplit(wangBrca$BRCA.status, "\\."), function(x) x[[1]])
  
  wangBrca$LOH = "No"
  for(i in 1:nrow(wangBrca)) {
    thisRow = wangBrca[i,]
    thisCN = cnState$consequence[ cnState$gene == thisRow$Gene & 
                                    substr(cnState$sample,1,12) == thisRow$Case_ID ]
    # Avoid no find
    if(length(thisCN)==0) next
    # Avoid potential AMPs and classify both LOH and DEL as "LOH" for this plot
    if(thisCN %in% c("LOH", "DEL")) wangBrca[i, "LOH"] = "LOH"
  }
  
  return(wangBrca)   
}


plotPanelA = function(dtExp, SIGNATURE = "CX3", CANCERS = "OV|BRCA", MAXXAXIS = 0.6, 
                      REMOVEYLABELS = FALSE, RESCALE = FALSE, VERBOSE = FALSE) {
  
  if(is.null(CANCERS)) {
    dtOVBRCA = dtExp[ dtExp$Signature == SIGNATURE, ]
    ## Replace split in condition based on cancer types
    dtOVBRCA$Status[ dtOVBRCA$Status == "WT BRCA1/2 TCGA" ] = "WT BRCA1/2"
  } else {
    dtOVBRCA = dtExp[ grepl(CANCERS, dtExp$Cancer) & dtExp$Signature == SIGNATURE, ]  
  }
  
  
  ## Rescale signature activities
  if(RESCALE) {  
    dtOVBRCA$Exposure = (dtOVBRCA$Exposure - min(dtOVBRCA$Exposure)) / 
      (max(dtOVBRCA$Exposure) - min(dtOVBRCA$Exposure))  
    NAMEY = paste("Rescaled", SIGNATURE, "activity")
  } else {
    NAMEY = paste(SIGNATURE, "activity")
  }
  
  ## Again simple stats for paper
  if(VERBOSE) {
    print(table(dtOVBRCA$Status))
    print(aggregate(Exposure ~ Status, data = dtOVBRCA, median))
    print(aggregate(Exposure ~ Status, data = dtOVBRCA, mean))
  }
  
  pB = ggplot(dtOVBRCA, aes(x = Status, y = Exposure)) +
    rasterise(geom_jitter(width = 0.2, size = 0.4, height = 0, alpha = 0.4, aes(colour = Cancer)), dpi = 600) + 
    geom_boxplot(outlier.alpha = NA, outlier.shape = NA, outlier.size = NA, alpha = 0.5) + 
    theme(legend.position = "none") + labs(y = NAMEY, x = "Genomic status of\n OV and BRCA samples") +
    coord_capped_flip(left = "both", bottom = "both", ylim = c(0, MAXXAXIS)) + scale_colour_manual(values = vCols)
  
  if(REMOVEYLABELS) { pB = pB + theme(axis.text.y = element_blank(), axis.title.y = element_blank()) }
  
  return(pB)
  
}

testSigPanelA = function(dtExp, RESCALE = TRUE, SIGNATURE = "CX3", CANCERS = "OV|BRCA", BACKGROUND = "WT BRCA1/2") {
  
  if(is.null(CANCERS)) {
    dtSig = dtExp[ dtExp$Signature == SIGNATURE, ]
    ## Replace split in condition based on cancer types
    dtSig$Status[ dtSig$Status == "WT BRCA1/2 TCGA" ] = "WT BRCA1/2"
  } else {
    dtSig = dtExp[ grepl(CANCERS, dtExp$Cancer) & dtExp$Signature == SIGNATURE, ]  
  }
  
  
  if(RESCALE) {  
    dtSig$Exposure = (dtSig$Exposure - min(dtSig$Exposure)) / 
      (max(dtSig$Exposure) - min(dtSig$Exposure))
  }
  
  ## Extract background values
  x = dtSig$Exposure[ dtSig$Status == BACKGROUND ]
  allLevels = levels(dtSig$Status)
  lSig = lapply(allLevels, function(thisStatus) {
    
    if(thisStatus == BACKGROUND | thisStatus == "WT BRCA1/2 TCGA" ) { return(NA) }
    
    y = dtSig$Exposure[ dtSig$Status == thisStatus ]
    if(length(y) < 2) { return(NA) } 
    tTest = t.test(x, y, var.equal = FALSE)
    
    ## Prepare output
    out = c(BACKGROUND, thisStatus, length(y), signif(tTest$estimate[1], 2), 
            signif(tTest$estimate[2], 2), signif(tTest$estimate[2] - tTest$estimate[1], 2),
            signif(tTest$statistic, 4), signif(tTest$p.value, 4))
    return(out)
    
  })
  
  dtSig = data.table(do.call(rbind, lSig))
  ## Remove NA row
  dtSig = dtSig[ ! rowSums(is.na(dtSig)) == 7,]
  colnames(dtSig) = c("Background", "Status", "NStatus", "MeanBack", "MeanStatus", "MeanDiff", 
                      "t", "pVal")
  dtSig$NStatus = as.numeric(dtSig$NStatus)
  dtSig$MeanBack = as.numeric(dtSig$MeanBack)
  dtSig$MeanStatus = as.numeric(dtSig$MeanStatus)
  dtSig$MeanDiff = as.numeric(dtSig$MeanDiff)
  dtSig$t = as.numeric(dtSig$t)
  dtSig$pVal = as.numeric(dtSig$pVal)
  
  dtSig$pAdj = p.adjust(dtSig$pVal, method = "BH")
  
  return(dtSig)
  
}


## Scatter plot with linear fit
plotActivityVSFeature = function(dtExp, SIGNATURES = c("CX3", "CX5", "CX2"), RESCALE = FALSE, FEATURE = "HRD_Score", PLOTNAME = "HRD Score", MAXXAXIS = 0.6) {
  
  
  ## Make plots for three signatures
  lPlots = lapply(SIGNATURES, function(thisSig) {
    
    dtSig = dtExp[ dtExp$Signature == thisSig, ]
    
    if(RESCALE) {  
      dtSig$Exposure = (dtSig$Exposure - min(dtSig$Exposure)) / (max(dtSig$Exposure) - min(dtSig$Exposure))  
      NAMEX = paste("Rescaled", thisSig, "activity")
    } else {
      NAMEX = paste(thisSig, "activity")
    }
    
    # Linear model for showing values inside plot
    fm = as.formula(paste("Exposure ~", FEATURE))
    lmFeat = lm(fm, data = dtSig)
    
    YPOS = max(dtSig[,..FEATURE], na.rm = TRUE)/2
    
    pPlot = ggplot(dtSig, aes_string(x = "Exposure", y = FEATURE)) + 
      # rasterise(geom_point(size = 0.5, alpha = 0.35), dpi = 600) +
      geom_point(size = 0.5, alpha = 0.35) +
      geom_smooth(method='lm', colour = "red") + labs(x = NAMEX, y = PLOTNAME) +
      annotate("text", x = 0.05, y = YPOS, size = 3,
               label = paste0("Intercept=", signif(coefficients(lmFeat)[1],2), 
                              "\nSlope=", signif(coefficients(lmFeat)[2],2),
                              "\np=", signif(summary(lmFeat)$coefficients[,4][2], digits = 2))) +
      coord_cartesian(xlim = c(0, MAXXAXIS))
    
  })
  
  ## Amend plots manually
  lPlots[[2]] = lPlots[[2]] + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
  lPlots[[3]] = lPlots[[3]] + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
  
  ## Merge and output
  # pOut = lPlots[[1]] + lPlots[[2]] + lPlots[[3]] + plot_layout(nrow = 1)
  return(lPlots)
  
} 

## Hex bin plot
plotActivityVSFeatureHex = function(dtExp, SIGNATURES = c("CX3", "CX5", "CX2"),  RESCALE = FALSE, FEATURE = "tp53_score", PLOTNAME = "TP53 Inactivation Score", 
                                    BREAKS = c(5, 10, 20, 30, 1000), BREAKLABELS = c("1-5", "5-10", "10-20", "20-30", "30+"),
                                    ADDLM = FALSE) {
  
  
  ## Make plots for three signatures
  lPlots = lapply(SIGNATURES, function(thisSig) {
    
    dtSig = dtExp[ dtExp$Signature == thisSig, ]
    
    if(RESCALE) {  
      ## Rescale signature activities
      dtSig$Exposure = (dtSig$Exposure - min(dtSig$Exposure)) / (max(dtSig$Exposure) - min(dtSig$Exposure))  
      NAMEX = paste("Rescaled", thisSig, "activity")
      
      ## Rescale feature (sometimes like PARPi7 not between [0,1])
      dtSig[[FEATURE]] = (dtSig[[FEATURE]] - min(dtSig[[FEATURE]], na.rm = TRUE)) / 
        (max(dtSig[[FEATURE]], na.rm = TRUE) - min(dtSig[[FEATURE]], na.rm = TRUE))  
      
    } else {
      NAMEX = paste(thisSig, "activity")
    }
    
    
    #### Prepare hex bin plot
    #### Taken from: https://stackoverflow.com/questions/46636473/consistent-hexagon-sizes-and-legend-for-manually-assignment-of-colors
    xbins = 25
    ## Identify bounds of created data
    minVal = min(dtSig$Exposure, dtSig[[FEATURE]], na.rm = TRUE)
    maxVal = max(dtSig$Exposure, dtSig[[FEATURE]], na.rm = TRUE)
    maxRange = c(minVal, maxVal)
    ## Add a buffer so that the last hex with data doesn't overlap margins
    buffer = (maxRange[2] - maxRange[1]) / (xbins / 2)
    ## Get our data into shape
    dtSig$factor = as.factor(1)
    bindata = dtSig[,c("Exposure", ..FEATURE, "factor")]
    
    ## Convert to hexbin object: object with number of counts in each cell with at least one count
    h = hexbin(bindata, xbins = xbins, IDs = TRUE, xbnds = maxRange, ybnds = maxRange)
    
    
    ### Convert back into an accessible data.frame
    ## Hextapply applies to each hex. Get ID of bin (names) and their respective point count (value)
    counts = hexTapply(h, bindata$factor, table)
    ## Convert to matrix with one row
    counts = t(simplify2array(counts))
    ## Melt into data frame: factor 1 (column 1), ID of hexbin (column 2), number of points in this bin (column 3)
    counts = melt(counts)
    colnames (counts)  <- c ("factor", "ID", "counts")
    counts$factor = as.factor(counts$factor)
    
    
    ## Get x and y coordinates from hexbin object and store in data.frame
    hexdf = data.frame(hcell2xy(h), ID = h@cell)
    ## Master object for plotting: ID of bin with number of points and x/y position
    hexdf <- merge(counts, hexdf)
    
    
    ## Create colour scale (make some more for more colour options)
    # clrs = brewer.pal(length(BREAKS) + 3, "Blues")
    clrs = brewer.pal(length(BREAKS)+1, "Blues")
    ## Chose darker colours (needs to be one larger than number of breaks for "Inf" and "0" category)
    # clrs = clrs[3:length(clrs)]
    
    
    
    ### Create new bin variable
    ## Add "0" to breaks for correct starting bin
    all_breaks <- c(0, BREAKS)
    ## Get number (or ID) for each break
    breaks_n <- 1:length(all_breaks)
    ## Function for identifying the break category each bin ends up in
    get_break_n <- function(n) {
      break_idx <- max(which((all_breaks - n) < 0))
      breaks_n[break_idx]
    }
    ## Apply function to master data frame
    hexdf$bin <- sapply(hexdf$counts, get_break_n)
    
    
    ## Create final plot
    pPlot = ggplot(hexdf, aes(x=x, y=y)) +
      geom_hex(stat="identity", aes(fill=bin)) +
      scale_fill_gradientn(colors=rev(clrs[-1]),
                           guide="legend",
                           labels=BREAKLABELS,
                           name="Count") +
      labs(x = NAMEX, y = PLOTNAME) +
      coord_fixed(xlim = c(0, (maxRange[2]+buffer)),
                  ylim = c(0, (maxRange[2]+buffer))) +
      theme(aspect.ratio=1)
      
    
    ## Add linear model if wished
    if(ADDLM) {
      
      # Linear model for showing values inside plot
      fm = as.formula(paste("Exposure ~", FEATURE))
      lmFeat = lm(fm, data = dtSig)
      
      YPOS = max(dtSig[,..FEATURE], na.rm = TRUE)/2
      
      pPlot = pPlot + geom_smooth(method='lm', colour = "red") +
        annotate("text", x = 0.05, y = YPOS, size = 3,
                 label = paste0("Intercept=", signif(coefficients(lmFeat)[1],2),
                                "\nSlope=", signif(coefficients(lmFeat)[2],2),
                                "\np=", signif(summary(lmFeat)$coefficients[,4][2], digits = 2)))
      
    }
    
    return(pPlot)
    
  })
  
  ## Amend plots manually
  lPlots[[2]] = lPlots[[2]] + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
  lPlots[[3]] = lPlots[[3]] + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
  
  ## Merge and output
  # pOut = lPlots[[1]] + lPlots[[2]] + lPlots[[3]] + plot_layout(nrow = 1)
  return(lPlots)
  
} 





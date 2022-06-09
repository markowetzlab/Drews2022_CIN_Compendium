createRandString<- function() {
    digits = 0:9
    v = c(sample(LETTERS, 2, replace = TRUE),
          sample(digits, 2, replace = TRUE),
          sample(LETTERS, 2, replace = TRUE))
    return(paste0(v, collapse = ""))
}

# # Works only for Mclust objects
# calcEntau = function(thisModel) {
#     matZ = thisModel$z
#     logMatZ = log(matZ)
#     logMatZ[ is.infinite(logMatZ) ] = -1000
#     entau = -sum( matZ*logMatZ )
#     return(entau)
# }
# 
# # Works only for Mclust objects
# iclbic = function(thisModel) {
#     entau = calcEntau(thisModel)
#     iclbic = -2*thisModel$loglik + 2*entau + thisModel$df*log(thisModel$n)
#     return(iclbic)
# }

fitComponent2 = function(dat,dist="norm",seed=77777,model_selection="BIC",min_prior=0.001,niter=1000,nrep=1,min_comp=2,max_comp=10) {
 
    thisRange = min_comp:max_comp
    set.seed(seed)
    
    if(dist=="norm") {
        if( model_selection == "BIC" ) {
            # "V" = mixture of Gaussians with variable variance
            # "E" = mixture of Gaussians with identical variance <- not really useful for our data. Leaving it out.
            return( Mclust(dat, modelNames = "V", G = thisRange) )
        } else {
            ## ICL
            # mclustICL to loop over range of K and determine K based on ICL.
            # Problem is that this functions outputs a "mclustICL" not an "Mclust" object. This is shit.
            # Hence, I need to rerun Mclust function with the optimal K to have a "Mclust" object.
            iclModel = mclustICL(dat, G = thisRange)
            modelInfo = strsplit( names(summary( iclModel ))[1], "," )[[1]]
            thisModelName = modelInfo[1]
            thisG = as.numeric( modelInfo[2] )
            print(paste("ICL-based model decision determined", thisG, "components."))
            return( Mclust(dat, modelNames = thisModelName, G = thisG) )
        }
    } else {
        # Poisson mixture model
        control<-new("FLXcontrol")
        control@minprior<-min_prior
        control@iter.max<-niter
        
        if(min_comp==max_comp) {
            fit = flexmix::flexmix(dat ~ 1, model=flexmix::FLXMCmvpois(), k=min_comp, control=control)
        } else {
            fit = flexmix::stepFlexmix(dat ~ 1,model = flexmix::FLXMCmvpois(), k=min_comp:max_comp, nrep=nrep, control=control)
            fit = flexmix::getModel(fit, which = model_selection)
        }
        
        return( fit )
    }
    
}

fitComponent<-function(dat,dist="norm",seed=77777,model_selection="BIC",min_prior=0.001,niter=1000,nrep=1,min_comp=2,max_comp=10) {
    
    control<-new("FLXcontrol")
    control@minprior<-min_prior
    control@iter.max<-niter
    
    set.seed(seed)
    if(dist=="norm") {
        fit = Mclust(dat, G = min_comp:max_comp)
        
    } else if(dist=="pois") {
        if(min_comp==max_comp) {
            fit = flexmix::flexmix(dat ~ 1, model=flexmix::FLXMCmvpois(), k=min_comp, control=control)
        } else {
            fit = flexmix::stepFlexmix(dat ~ 1,model = flexmix::FLXMCmvpois(), k=min_comp:max_comp, nrep=nrep, control=control)
            fit = flexmix::getModel(fit, which = model_selection)
        }
    }
    
    return( fit )
}


fitComponentGM<-function(dat,dist="norm",seed=77777,model_selection="BIC",min_prior=0.001,niter=1000,nrep=1,min_comp=2,max_comp=10)
{
    control<-new("FLXcontrol")
    control@minprior<-min_prior
    control@iter.max<-niter
    
    set.seed(seed)
    if(dist=="norm")
    {
        if(min_comp==max_comp)
        {
            fit<-flexmix::flexmix(dat ~ 1,model=flexmix::FLXMCnorm1(),k=min_comp,control=control)
        }else{
            fit<-flexmix::stepFlexmix(dat ~ 1,model = flexmix::FLXMCnorm1(),k=min_comp:max_comp,nrep=nrep,control=control)
            fit<-flexmix::getModel(fit,which=model_selection)
        }
        
    }else if(dist=="pois")
    {
        if(min_comp==max_comp)
        {
            fit<-flexmix::flexmix(dat ~ 1,model=flexmix::FLXMCmvpois(),k=min_comp,control=control)
        }else{
            fit<-flexmix::stepFlexmix(dat ~ 1,model = flexmix::FLXMCmvpois(),k=min_comp:max_comp,nrep=nrep,control=control)
            fit<-flexmix::getModel(fit,which=model_selection)
        }
    }
    fit
}


getParams = function(var, lenRuns) {
    if(length(var) == 1 ) {
        rep(var, lenRuns)
    } else if ( length(var) == lenRuns) {
        var        
    } else if (lenRuns > length(var)) {
        c(var, rep(var[length(var)], lenRuns-length(var)))
    } else {
        var[1:lenRuns]
    }
}


calculateSumOfPosteriors<-function(CN_feature,components,name, rowIter = 1000, cores = 1)
{
    
    if(cores > 1){
        require(foreach)
        require(doMC)

        len = dim(CN_feature)[1]
        iters = floor( len / rowIter )
        lastiter = iters[length(iters)]

        registerDoMC(cores)
        curr_posterior = foreach( i=0:iters, .combine=rbind) %dopar% {
            start = i*rowIter+1
            if(i != lastiter) { end = (i+1)*rowIter } else { end = len }
                flexmix::posterior(components,data.frame(dat=as.numeric(CN_feature[start:end,2])))
        }
    } else {
        curr_posterior<-flexmix::posterior(components,data.frame(dat=as.numeric(CN_feature[,2])))
    }

    mat<-cbind(CN_feature,curr_posterior)
    posterior_sum<-c()
    mat = mat[,c(-2)]
    posterior_sum = aggregate(. ~ ID, mat, sum)
    rownames(posterior_sum) = posterior_sum$ID
    posterior_sum$ID = NULL

    params<-flexmix::parameters(components)
    if(!is.null(nrow(params)))
    {
        posterior_sum<-posterior_sum[,order(params[1,])]
    } else {
        posterior_sum<-posterior_sum[,order(params)]
    }
    colnames(posterior_sum)<-paste0(name,1:ncol(posterior_sum))
    posterior_sum
}

getSegsize<-function(abs_profiles, rmNorm = FALSE)
{
    out<-c()
    samps<-getSampNames(abs_profiles)
    for(i in samps)
    {
        if(class(abs_profiles)=="QDNAseqCopyNumbers")
        {
            segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        }
        else
        {
            segTab<-abs_profiles[[i]]
            colnames(segTab)[4]<-"segVal"
        }
        # First tap.
        segTab$segVal = as.numeric(segTab$segVal)
        # If wished, don't consider normal segments
        if(rmNorm) { segTab = segTab[ segTab$segVal != 2, ] }
        # Avoiding potential artefact
        segTab$segVal[segTab$segVal<0]<-0
        seglen<-segTab$end-segTab$start
        seglen<-seglen[seglen>0]
        # Double tap.
        out<-rbind(out,cbind(ID=rep(i,length(seglen)),value=as.numeric(seglen)))
    }
    rownames(out)<-NULL
    data.frame(out,stringsAsFactors = F)
}

getBPnum<-function(abs_profiles,chrlen, SIZE = 10000000)
{
    out<-c()
    samps<-getSampNames(abs_profiles)
    for(i in samps)
    {
        if(class(abs_profiles)=="QDNAseqCopyNumbers")
        {
            segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        }else
        {
            segTab<-abs_profiles[[i]]
            colnames(segTab)[4]<-"segVal"
        }
        chrs<-unique(segTab$chromosome)
        allBPnum<-c()
        for(c in chrs)
        {
            currseg<-segTab[segTab$chromosome==c,]
            intervals<-seq(1,chrlen[chrlen[,1]==paste0("chr",c),2]+SIZE,SIZE)
            res <- hist(as.numeric(currseg$end[-nrow(currseg)]),breaks=intervals,plot=FALSE)$counts
            allBPnum<-c(allBPnum,res)
        }
        # Make sure it's really numeric
        out<-rbind(out,cbind(ID=rep(i,length(allBPnum)),value=as.numeric(allBPnum)))
    }
    rownames(out)<-NULL
    data.frame(out,stringsAsFactors = F)
}

getOscilation<-function(abs_profiles,chrlen)
{
    out<-c()
    samps<-getSampNames(abs_profiles)
    for(i in samps)
    {
        if(class(abs_profiles)=="QDNAseqCopyNumbers")
        {
            segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        }else
        {
            segTab<-abs_profiles[[i]]
            colnames(segTab)[4]<-"segVal"
        }
        chrs<-unique(segTab$chromosome)
        oscCounts<-c()
        for(c in chrs)
        {
            currseg<-as.numeric(segTab$segVal[segTab$chromosome==c])
            currseg<-round(as.numeric(currseg))
            if(length(currseg)>3)
            {
                prevval<-currseg[1]
                count=0
                for(j in 3:length(currseg))
                {
                    if(currseg[j]==prevval&currseg[j]!=currseg[j-1])
                    {
                        count<-count+1
                    }else{
                        oscCounts<-c(oscCounts,count)
                        count=0
                    }
                    prevval<-currseg[j-1]
                }
            }
        }
        # Make sure it's really numeric
        out<-rbind(out,cbind(ID=rep(i,length(oscCounts)),value=as.numeric(oscCounts)))
        if(length(oscCounts)==0)
        {
            out<-rbind(out,cbind(ID=i,value=0))
        }
    }
    rownames(out)<-NULL
    data.frame(out,stringsAsFactors = F)
}

getCentromereDistCounts<-function(abs_profiles,centromeres,chrlen)
{
    out<-c()
    samps<-getSampNames(abs_profiles)
    for(i in samps)
    {
        if(class(abs_profiles)=="QDNAseqCopyNumbers")
        {
            segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        }else
        {
            segTab<-abs_profiles[[i]]
            colnames(segTab)[4]<-"segVal"
        }
        chrs<-unique(segTab$chromosome)
        all_dists<-c()
        for(c in chrs)
        {
            if(nrow(segTab)>1)
            {
                starts<-as.numeric(segTab[segTab$chromosome==c,2])[-1]
                segstart<-as.numeric(segTab[segTab$chromosome==c,2])[1]
                ends<-as.numeric(segTab[segTab$chromosome==c,3])
                segend<-ends[length(ends)]
                ends<-ends[-length(ends)]
                centstart<-as.numeric(centromeres[substr(centromeres[,2],4,5)==c,3])
                centend<-as.numeric(centromeres[substr(centromeres[,2],4,5)==c,4])
                chrend<-chrlen[substr(chrlen[,1],4,5)==c,2]
                ndist<-cbind(rep(NA,length(starts)),rep(NA,length(starts)))
                ndist[starts<=centstart,1]<-(centstart-starts[starts<=centstart])/(centstart-segstart)*-1
                ndist[starts>=centend,1]<-(starts[starts>=centend]-centend)/(segend-centend)
                ndist[ends<=centstart,2]<-(centstart-ends[ends<=centstart])/(centstart-segstart)*-1
                ndist[ends>=centend,2]<-(ends[ends>=centend]-centend)/(segend-centend)
                ndist<-apply(ndist,1,min)
                
                all_dists<-rbind(all_dists,sum(ndist>0))
                all_dists<-rbind(all_dists,sum(ndist<=0))
            }
        }
        if(nrow(all_dists)>0)
        {
            # Make sure it's really numeric
            out<-rbind(out,cbind(ID=i,ct1=as.numeric(all_dists[,1])))
        }
    }
    rownames(out)<-NULL
    data.frame(out,stringsAsFactors = F)
}


getChangepointCN<-function(abs_profiles, allowedError = 0.1, rmNorm = FALSE)
{
    out<-c()
    samps<-getSampNames(abs_profiles)
    for(i in samps)
    {
        if(class(abs_profiles)=="QDNAseqCopyNumbers") {
            segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        } else {
            segTab<-abs_profiles[[i]]
            colnames(segTab)[4]<-"segVal"
        }
        segTab$segVal = as.numeric(segTab$segVal)
        segTab$segVal[segTab$segVal<0]<-0
        chrs<-unique(segTab$chromosome)
        allcp<-c()
        for(c in chrs) {
            currseg<-as.numeric(segTab$segVal[segTab$chromosome==c])
            firstSeg = abs(2 - currseg[1] )
            # As we look only at the left end of a CNA, we might miss a changepoint at the beginning of the p-arm
            # That's why we check manually but only regard this value if it is higher than an allowed error rate.
            if(firstSeg <= allowedError) {
                theseChanges = abs(currseg[-1]-currseg[-length(currseg)])
                if(rmNorm) { theseChanges = theseChanges[ currseg[-1] != 2 ] }
                allcp<-c(allcp, theseChanges)
            } else {
                theseChanges = c( firstSeg, abs(currseg[-1]-currseg[-length(currseg)]) )
                if(rmNorm) { theseChanges = theseChanges[ currseg != 2 ] }
                allcp<-c(allcp, theseChanges)
            }
            
        }
        if(length(allcp)==0) {
            allcp<-0 #if there are no changepoints
        }
        # Make sure it's really numeric
        out<-rbind(out,cbind(ID=rep(i,length(allcp)),value=as.numeric(allcp)))
    }
    rownames(out)<-NULL
    data.frame(out,stringsAsFactors = F)
}

getCN<-function(abs_profiles, rmNorm = FALSE)
{
    out<-c()
    samps<-getSampNames(abs_profiles)
    for(i in samps)
    {
        if(class(abs_profiles)=="QDNAseqCopyNumbers")
        {
            segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        }
        else
        {
            segTab<-abs_profiles[[i]]
            colnames(segTab)[4]<-"segVal"
        }
        segTab$segVal[as.numeric(segTab$segVal)<0]<-0
        # If wished, don't consider normal segments
        if(rmNorm) { segTab = segTab[ segTab$segVal != 2, ] }
        cn<-as.numeric(segTab$segVal)
        # Make sure it's really numeric.
        out<-rbind(out,cbind(ID=rep(i,length(cn)),value=as.numeric(cn)))
    }
    rownames(out)<-NULL
    data.frame(out,stringsAsFactors = F)
}

getSampNames<-function(abs_profiles)
{
    if(class(abs_profiles)=="QDNAseqCopyNumbers")
    {
        samps<-colnames(abs_profiles)
    }
    else
    {
        samps<-names(abs_profiles)
    }
    samps
}

getSegTable<-function(x)
{
    dat<-x
    sn<-Biobase::assayDataElement(dat,"segmented")
    fd <- Biobase::fData(dat)
    fd$use -> use
    fdfiltfull<-fd[use,]
    sn<-sn[use,]
    segTable<-c()
    for(c in unique(fdfiltfull$chromosome))
    {
        snfilt<-sn[fdfiltfull$chromosome==c]
        fdfilt<-fdfiltfull[fdfiltfull$chromosome==c,]
        sn.rle<-rle(snfilt)
        starts <- cumsum(c(1, sn.rle$lengths[-length(sn.rle$lengths)]))
        ends <- cumsum(sn.rle$lengths)
        lapply(1:length(sn.rle$lengths), function(s) {
            from <- fdfilt$start[starts[s]]
            to <- fdfilt$end[ends[s]]
            segValue <- sn.rle$value[s]
            c(fdfilt$chromosome[starts[s]], from, to, segValue)
        }) -> segtmp
        segTableRaw <- data.frame(matrix(unlist(segtmp), ncol=4, byrow=T),stringsAsFactors=F)
        segTable<-rbind(segTable,segTableRaw)
    }
    colnames(segTable) <- c("chromosome", "start", "end", "segVal")
    segTable
}


getPloidy<-function(abs_profiles)
{
  out<-c()
  samps<-getSampNames(abs_profiles)
  for(i in samps)
  {
    if(class(abs_profiles)=="QDNAseqCopyNumbers")
    {
      segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
    }
    else
    {
      segTab<-abs_profiles[[i]]
      colnames(segTab)[4]<-"segVal"
    }
    segLen<-(as.numeric(segTab$end)-as.numeric(segTab$start))
    ploidy<-sum((segLen/sum(segLen))*as.numeric(segTab$segVal))
    out<-c(out,ploidy)
  }
  data.frame(out,stringsAsFactors = F)
}


normaliseMatrix<-function(signature_by_sample,sig_thresh=0.01)
{
    norm_const<-colSums(signature_by_sample)
    sample_by_signature<-apply(signature_by_sample,1,function(x){x/norm_const})
    sample_by_signature<-apply(sample_by_signature,1,lower_norm,sig_thresh)
    signature_by_sample<-t(sample_by_signature)
    norm_const<-apply(signature_by_sample,1,sum)
    sample_by_signature<-apply(signature_by_sample,2,function(x){x/norm_const})
    signature_by_sample<-t(sample_by_signature)
    signature_by_sample
}

lower_norm<-function(x,sig_thresh=0.01)
{
    new_x<-x
    for(i in 1:length(x))
    {
        if(x[i]<sig_thresh)
        {
            new_x[i]<-0
        }
    }
    new_x
}

# Taken from: https://www.r-bloggers.com/identifying-the-os-from-r/
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

listNMF = function(mat, ranks = ranks, seed = seed, 
    iter = iter, nmfalg = nmfalg, cores = cores) {
    
    if( nmfalg == "snmf/r" ) {
        require(NMF)
        lNMFs = list()
        index = 1
        dispScore = 1
        while( dispScore == 1 && ! is.na(ranks[index])  ) {
            print(paste("rank:", ranks[index]))
            thisRun = NMF::nmf(mat, ranks[index], seed=seed, nrun=iter, method=nmfalg, .opt=list(paste0("vp", cores) ) )
            dispScore = dispersion(thisRun)

            # Store results and prepare next iteration
            lNMFs[[ index ]] = thisRun
            index = index + 1
        }
        return( lNMFs )
    } else {
        require(foreach)
        require(NMF)

        lNMFs = foreach( r = ranks ) %do% {
            print(paste("rank:", r))
            NMF::nmf(mat, r, seed=seed, nrun=iter, method=nmfalg, .opt=list(paste0("vp", cores) ) )
        }
        return( lNMFs )
    }
}

calcScore = function(nmfFit, score, mat=NULL) {
    require(NMF)

    if(score == "Cophenetic") {
        return( cophcor(nmfFit) )
    } else if (score == "Sparseness") {
        # [1] is "base" or matrix W
        # [2] is "coef" or matrix H
        # We are interested in the sparseness of the signatures, but as we have transposed our matrix, we want the sparseness of H.
        return( sparseness(nmfFit)[2] )
    } else if (score == "RSS") {
        return( rss(nmfFit, mat) )
    } else if (score == "Dispersion.coef") {
        return( dispersion(nmfFit) )
    }
}

getIndScore = function(lNMFs, ranks, mat) {
    len = length(lNMFs)
    scores = c("Cophenetic", "Sparseness", "RSS", "Dispersion.coef")

    lScores = lapply(lNMFs, function(thisFit) {
        thisRank = as.character( dim(thisFit@fit@H)[1] )
        lScores=list()
        sapply(scores, function(s) {
            calcScore( thisFit, score = s, mat )
        })

    })

    dfScores = as.data.frame(lScores)
    names(dfScores) = ranks
    dfScores = melt(t(dfScores))
    colnames(dfScores) = c("Rank", "Score", "Value")
    return(dfScores)
}

getScores = function(lNMFs, lNMFsran, ranks, mat, randomMat=TRUE) {

    if( isTRUE(randomMat) ) {
        dfScoresReal = getIndScore(lNMFs, ranks, mat)
        dfScoresReal$Condition = "Observed"
        dfScoresRand = getIndScore(lNMFsran, ranks, mat)
        dfScoresRand$Condition = "Random"
        dfScores = rbind(dfScoresReal, dfScoresRand)

        return(dfScores)
    } else {
        dfScoresReal = getIndScore(lNMFs, ranks, mat)
        dfScoresReal$Condition = "Observed"
        return(dfScoresReal)
    }

}

plotScores = function(dfScores) {

    ggplot(dfScores, aes(x=Rank, y=Value, colour=Condition)) + geom_line() + facet_grid(Score ~ ., scales = "free")

}

decideNumSig = function(tcga.numSigs, decisionScore = "Sparseness.coef" ) {
  dfSigs = tcga.numSigs$scores
  if( decisionScore == "Sparseness.coef" ) {
    sigsObs = dfSigs[ dfSigs$Condition == "Observed" & 
                      dfSigs$Score == "Sparseness.coef", ]
    sigsObs$Random = dfSigs[ dfSigs$Condition == "Random" & 
                      dfSigs$Score == "Sparseness.coef", "Value" ]
    sigsObs$Diff = sigsObs$Value - sigsObs$Random
    
    ## Some code for eventual in-depth checking of the potential number of signatures
    # potRank = as.numeric( sigsObs[ sigsObs$Diff < 0, "Rank" ][1] ) -1
    # # check if chosen rank is actually in the screened rank. if not, maybe mistake as small ranks tend to have very low sparse values compared to random matrices.
    # # check if it's the first rank but the next is not anymore, then chose to ignore it.
    # minSig = min(sigsObs$Rank)
    # if( potRank <= minSig ) {
    #     # the rank +1 is the one with the negative difference
    #     # the rank +2 is the one after that. If that is still negative, then take it, otherwise ignore them and chose new
    #     if(sigsObs[sigsObs$Rank == potRank + 2,"Diff"] > 0){
    #         sigsObs = sigsObs[-which( sigsObs$Rank %in% (potRank+1)),]
    #         potRank = as.numeric( sigsObs[ sigsObs$Diff < 0, "Rank" ][1] ) -1
    #     }
    # }

    # # check if rank exists for this criterion, otherwise take the largest dip in cophenetic score
    # if( is.na(potRank) ) {
    #     sigsObs = dfSigs[ dfSigs$Condition == "Observed" & 
    #                   dfSigs$Score == "Cophenetic", ]

    #     diffs = diff(sigsObs$Value)
    #     # Rank with largest dip is actually one rank too far but as the difference is always one rank shorter, it works out.
    #     potRank = which(diffs %in% min(diffs))
    # }
    
    return( as.numeric( sigsObs[ sigsObs$Diff < 0, "Rank" ][1] ) -1 )
  } else if ( decisionScore == "Dispersion.coef" ) {
    obs = dfSigs[ dfSigs$Condition == "Observed" & 
                      dfSigs$Score == "Dispersion.coef", "Value" ]
    # First TRUE is actually one rank too far but as the difference is always one rank shorter, it works out.
    return( dfSigs[ diff(obs) < 0, "Rank" ][1] )
  }
}


# Edited from Lynch 2016 (https://f1000research.com/articles/5-1253/v1)
getRinv = function(signatures.ref) {
    return( backsolve(chol(signatures.ref %*% t(signatures.ref)), 
                          diag(dim(signatures.ref)[1])) )
}

QPsig = function(tumour.ref, samplerow, signatures.ref, Rinverse){
    ## AL:
    # we normalize the observations so that they sum to 1
    obs<-as.numeric(tumour.ref[samplerow,]/sum(tumour.ref[samplerow,]))
    
    ## RMD:
    # deleted: conversion to matrix as already happened before.

    ## AL:
    # we use the factorized version of solve.QP - 
    # see the helpfile of that function for details of the required form
    # otherwise we would set Dmat = signatures.ref %*% t(signatures.ref) as indicated 
    # in the article
    
    ## RMD:
    # deleted: calculation of Rinverse as outsourced to function prior to this call
    
    ## AL:
    # we also need to define the linear part of the function that we are minimizing
    dvec <- (obs) %*% t(signatures.ref)
    
    ## AL:
    # we have one constraint that the sum of weights is equal to 1
    # we have more constraints that each weight is positive
    Amat <- cbind(rep(1,dim(Rinverse)[1]), diag(dim(Rinverse)[1]))
    bvec <- c(1,rep(0,dim(Rinverse)[1]))
    
    ## AL:
    # we now call the solve.QP function from the quadprog library
    myQP<-quadprog::solve.QP(Dmat = Rinverse, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1, 
    factorized = TRUE)
    return(myQP$solution)    
}

plotpdfpng = function( mSig1VsSig2, strNumSigs1, strNumSigs2, outputNameBase, titley, titleA, pdf=TRUE) 
{
    library("Cairo")
    if(pdf) {
        outputName=paste0(outputNameBase, "nnls_", titley, "-", strNumSigs1, "_vs_", titleA, "-", strNumSigs2, ".pdf")
        pdf(outputName, paper="a4r", width=0, height=0)
    } else {
        outputName=paste0(outputNameBase, "nnls_", titley, "-", strNumSigs1, "_vs_", titleA, "-", strNumSigs2, ".png")
        CairoPNG(file = outputName, width=1000, height=1000)
    }
    corrplot(t(mSig1VsSig2), method="color", cl.lim=c(min(mSig1VsSig2), max(mSig1VsSig2)), is.corr=FALSE,
             title = paste("argmin_x ||Ax-y||_2 with x>=0\n y =", titley, "and A =", titleA),
             mar = c(1, 1, 5, 1) )
    dev.off()

}

# Correlate signatures with same dimensions of features
corrSigs = function(listOV, listNaive, outputNameBase, smallFilter = 0) {
    library("nnls")
    library("corrplot")

    titley = deparse(substitute(listOV))
    titleA = deparse(substitute(listNaive))

    # Loop over first set of signatures. The ones you want to know whether they are made up of elements of the second set of signatures.
    listAllSigsNnls=lapply(listOV, function(sigs1) {

        strNumSigs1=dim(sigs1)[2]
        # Loop over second set of signatures
        listSig1VsSig=lapply(listNaive, function(sigs2) {

            strNumSigs2=dim(sigs2)[2]
            # Do the non-negative least squares for each signature from 1 against the matrix of signatures 2.
            # The idea is that a signature Y can be displayed as linear combination of matrix A, hence: Ax = y
            # The algorithm nnls() now tries to minimise the following equation:
            # argmin_x ||Ax-y||_2 with x >= 0
            # The vector x (or b) is called coefficients. The nnls is constraining the coefficients to be non-negative.
            mSig1VsSig2=apply(sigs1, 2, function(s1) {
                coefs=coef(nnls(A = sigs2, b = s1))
                return(coefs)
            })

            mSig1VsSig2[mSig1VsSig2 < smallFilter] = 0

                        # Print as pdf and png to the output
            plotpdfpng( mSig1VsSig2, strNumSigs1, strNumSigs2, outputNameBase, titley=titley, titleA=titleA, pdf=TRUE) 
            plotpdfpng( mSig1VsSig2, strNumSigs1, strNumSigs2, outputNameBase, titley=titley, titleA=titleA, pdf=FALSE)

            return(mSig1VsSig2)
        } )
        return(listSig1VsSig)
    } )
    return(listAllSigsNnls)
}

library("NMF")
library("quadprog")
#library("ComplexHeatmap")
#library("circlize")
# for the manual parallelisation and also used by NMF for their internal parallelisation
library("foreach")
# needed for faster nmf runs based on the developer's own recommendations:
# https://cran.r-project.org/web/packages/NMF/vignettes/NMF-vignette.pdf
# page 16 ff.
library("doParallel")
# Libraries "bigmemory" and "synchronicity" are recommended for non-Windows machines
# for remodelling data structures
library("reshape2")
# for plotting
library("ggplot2")
#library("Cairo")
# Correlation of signatures
library("nnls")
library("corrplot")


locationOfThisScript = function() # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
{
    this.file = NULL
    # This file may be 'sourced'
    for (i in -(1:sys.nframe())) {
        if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
    }
    
    if (!is.null(this.file)) return(dirname(this.file))
    
    # But it may also be called from the command line
    cmd.args = commandArgs(trailingOnly = FALSE)
    cmd.args.trailing = commandArgs(trailingOnly = TRUE)
    cmd.args = cmd.args[seq.int(from=1, length.out=length(cmd.args) - length(cmd.args.trailing))]
    res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)
    
    # If multiple --file arguments are given, R uses the last one
    res = tail(res[res != ""], 1)
    if (0 < length(res)) return(dirname(res))
    
    # Both are not the case. Maybe we are in an R GUI?
    return(NULL)
}
this_path<-locationOfThisScript()
source(paste(this_path,"helper_functions.R",sep="/"))

plotSignatures = function(tcga.sigs, file, width = 9, height = 7) {

    # columns are features and rows are signatures => output of generateSignatures
    rownames(tcga.sigs) = paste0("Sig", seq(1:nrow(tcga.sigs)))
    sigsNormFeat = t( normaliseMatrix(tcga.sigs, sig_thresh = 0) )
    sigsNormSigs = normaliseMatrix(t ( tcga.sigs ), sig_thresh = 0)

    pdf(file, width = width, height = height, onefile=FALSE)
    par(mfrow = c(1,2))
    aheatmap(sigsNormFeat, Rowv = NA, Colv = NA, main = "Normalised by feature", color = colorRampPalette(c("white", "darkblue"))(50) )
    aheatmap(sigsNormSigs, Rowv = NA, Colv = NA, main = "Normalised by signature", color = colorRampPalette(c("white", "darkblue"))(50) )
    dev.off()

}

plotExposures = function(tcga.exp, cnaPerSample, namesAndCancers, outfile, cancerCols) 
{
    tcga.expSort = tcga.exp[match(names(cnaPerSample), rownames(tcga.exp)),]

    dfCnaPerSample = data.frame("cnaPerSample" = cnaPerSample)
    namesAndCancers = namesAndCancers[ match(names(cnaPerSample), namesAndCancers[,1] ), ]
    colsToUse = cancerCols[ match(sort(unique(namesAndCancers[,3])), 
                                  cancerCols$tcga), "cols"]

    # Annotation plots
    ha = rowAnnotation(df = dfCnaPerSample)
    ht2 = Heatmap(namesAndCancers[,3], name = colnames(namesAndCancers)[3], 
                  show_row_names = FALSE, width = unit(5, "mm"), 
                  show_column_dend = FALSE, show_column_names = TRUE, col = colsToUse)

    # Heatmaps
    # Split by cancer
    ht1a = Heatmap(tcga.expSort, split = namesAndCancers[,3], name = "Exposure",
                  show_column_names = TRUE, column_names_gp = gpar(fontsize = 6),
                  row_title_rot = 0, row_title_gp = gpar(fontsize = 8),
                  show_row_names = FALSE, show_row_dend = FALSE, 
                  row_dend_reorder = FALSE, cluster_rows = FALSE,
                  col = colorRamp2(c(0, 1), c("white", "red")) )
    ht1aa = Heatmap(tcga.expSort, split = namesAndCancers[,3], name = "Exposure",
                  show_column_names = TRUE, column_names_gp = gpar(fontsize = 6),
                  row_title_rot = 0, row_title_gp = gpar(fontsize = 8),
                  show_row_names = FALSE, show_row_dend = FALSE, 
                  col = colorRamp2(c(0, 1), c("white", "red")) )
    # Split by cancer and colour maxed out over 0.5
    ht1b = Heatmap(tcga.expSort, split = namesAndCancers[,3], name = "Exposure",
                  show_column_names = TRUE, column_names_gp = gpar(fontsize = 6),
                  row_title_rot = 0, row_title_gp = gpar(fontsize = 8),
                  show_row_names = FALSE, show_row_dend = FALSE, 
                  row_dend_reorder = FALSE, cluster_rows = FALSE,
                  col = colorRamp2(c(0, 0.5, 1), c("white", "red", "red")) )
    ht1ba = Heatmap(tcga.expSort, split = namesAndCancers[,3], name = "Exposure",
                  show_column_names = TRUE, column_names_gp = gpar(fontsize = 6),
                  row_title_rot = 0, row_title_gp = gpar(fontsize = 8),
                  show_row_names = FALSE, show_row_dend = FALSE, 
                  col = colorRamp2(c(0, 0.5, 1), c("white", "red", "red")) )
    # No splitting, all samples in one big group
    ht1c = Heatmap(tcga.expSort, km = 1, name = "Exposure",
                  show_column_names = TRUE, column_names_gp = gpar(fontsize = 6),
                  row_title_rot = 0, row_title_gp = gpar(fontsize = 8),
                  show_row_names = FALSE, show_row_dend = FALSE, 
                  row_dend_reorder = FALSE, cluster_rows = FALSE,
                  col = colorRamp2(c(0, 0.5, 1), c("white", "red", "red")) )
    ht1ca = Heatmap(tcga.expSort, km = 1, name = "Exposure",
                  show_column_names = TRUE, column_names_gp = gpar(fontsize = 6),
                  row_title_rot = 0, row_title_gp = gpar(fontsize = 8),
                  show_row_names = FALSE, show_row_dend = FALSE,
                  col = colorRamp2(c(0, 0.5, 1), c("white", "red", "red")) )

    pdf(outfile, width = 7, height = 14)
    draw(ht1a+ha+ht2)
    draw(ht1aa+ha+ht2)
    draw(ht1b+ha+ht2)
    draw(ht1ba+ha+ht2)
    draw(ht1c+ha+ht2)
    draw(ht1ca+ha+ht2)
    dev.off()
}

quantifySignatures2<-function(mutCatalogue, sigCatalogue=NULL, method = "QP", threshold = 0.01, smallSampleNumbers=FALSE)
{
    if(is.null(sigCatalogue))
    {
        sigCatalogue<-readRDS("data/feat_sig_mat.rds")
    }

    if( method == "YAPSA" ) {
        
        ### YAPSA needs:
        ## Left matrix W        sigCatalogue        component (rows) by signature (cols)   <= HAVE
        ## Right matrix H       expCatalogue        signature (rows) by samples (cols)     <= WANT
        ## Full matrix V        mutCatalogue        components (rows) by samples (cols)    <= HAVE
        
        ### Assumptions: signatures < components < samples
        ## LCD(in_mutation_catalogue_df, in_signatures_df, in_per_sample_cutoff = 0)
        
        ## in_mutation_catalogue_df: A numeric data frame 'V' with 'n' rows and
        ##           'm' columns, 'n' being the number of features and 'm' being
        ##           the number of samples
        
        ## in_signatures_df: A numeric data frame 'W' with 'n' rows and 'l'
        ##           columns, 'n' being the number of features and 'l' being the
        ##           number of signatures

        # Meaning: less samples than components 
        if( isTRUE(smallSampleNumbers) ) {
            # samples in rows, so transpose matrix
            if( nrow(mutCatalogue) < ncol(mutCatalogue) ) { mutCatalogue = t( mutCatalogue ) }
        } else {
            # samples in rows, so transpose matrix
            if( nrow(mutCatalogue) > ncol(mutCatalogue) ) { mutCatalogue = t( mutCatalogue ) }
        }

        # signatures in rows, so transpose matrix
        if( nrow(sigCatalogue) < ncol(sigCatalogue) ) { sigCatalogue = t( sigCatalogue ) }

        # Signature columns have to be normalised to 1! not feature rows!
        expCatalogue = YAPSA::LCD( mutCatalogue, YAPSA:::normalize_df_per_dim( sigCatalogue,2 ) )
        expCatalogue = normaliseMatrix( expCatalogue, sig_thresh = threshold)
        return( t( expCatalogue ) )

    } else if ( method == "QP" ) {
        
        # For QPsig, each sample is calculated individually by inputting only the sample's feature vector and the signature catalogue
        # Assumption: signatures < components < samples

        # Meaning: less samples than components 
        if( isTRUE(smallSampleNumbers) ) {
            # samples in rows, so transpose matrix
            if( nrow(mutCatalogue) > ncol(mutCatalogue) ) { mutCatalogue = t( mutCatalogue ) }
        } else {
            # samples in rows, so transpose matrix
            if( nrow(mutCatalogue) < ncol(mutCatalogue) ) { mutCatalogue = t( mutCatalogue ) }
        }

        # signatures in rows, so transpose matrix
        if( nrow(sigCatalogue) > ncol(sigCatalogue) ) { sigCatalogue = t( sigCatalogue ) }
        # constraint of QPsig is normalised sigCatalogue matrix
        sigCatalogue = as.matrix( YAPSA:::normalize_df_per_dim( sigCatalogue,1 ) )
 
        expCatalogue = matrix(0, nrow = nrow(mutCatalogue), ncol = nrow(sigCatalogue) )
        Rinverse = getRinv(sigCatalogue)
        for(i in 1:nrow(mutCatalogue) ) { 
            expCatalogue[i,] = QPsig(tumour.ref = mutCatalogue, i, 
                                            signatures.ref = sigCatalogue,
                                            Rinverse = Rinverse )
        }

        rownames(expCatalogue) = rownames(mutCatalogue)
        colnames(expCatalogue) = paste0("s", seq(1:ncol(expCatalogue)) )
        # Remove signatures with too low exposure
        expCatalogue = t( normaliseMatrix( t( expCatalogue), sig_thresh = threshold) )
        return(expCatalogue)

    }

}

quantifySignatures<-function(sample_by_component,component_by_signature=NULL)
{
    if(is.null(component_by_signature))
    {
        component_by_signature<-readRDS(paste(this_path,"data/feat_sig_mat.rds",sep="/"))
    }
    signature_by_sample<-YAPSA::LCD(t(sample_by_component),
                                    YAPSA:::normalize_df_per_dim(component_by_signature,2))
    signature_by_sample<-normaliseMatrix(signature_by_sample)
    signature_by_sample
}

generateSignatures<-function(sample_by_component, nsig, iter=1000,
                             seed="nndsvd", nmfalg="snmf/r", cores=40, 
                             threshold = 0.01)
{
    modelSigs = NMF::nmf(sample_by_component, nsig, seed=seed, nrun=iter,
                         method=nmfalg, .opt = paste0("p", cores) )
    sigs = modelSigs@fit@H
    # failsafe: remove empty signatures
    sigs = sigs[ apply(sigs, 1, sum) != 0, ]
}

chooseNumberSignatures<-function(sample_by_component, outfile="numSigs.pdf", min_sig=3,
                                 max_sig=15, iter=30, cores=1, largeMatrix=FALSE,
                                 nmfalg="brunet",seed=77777, randomMat=TRUE )
{

    # Small matrices (500x40) can be easily done with largeMatrix=FALSE.
    # The estimation of the factorisation rank on larger matrices (12000x50) cannot be done via the nmfEstimateRank function due to the insane memory overhead in R. Hence the effort should be made to install additional packages ("foreach", "doParallel" and when on non-Windows "bigmemory" and "synchronicity") and largeMatrix should be set to TRUE.
    # 20 cores and 1GB memory should do for 12000x50. Runtime is around 850 seconds with peak memory usage of 11MB.

    # 30 - 50 iterations is enough to estimate the factorisation rank
    # Brunet2004; Hutchins2008


    # possible later implementation: NMF via GPU
    # https://www.rdocumentation.org/packages/nmfgpu4R/versions/0.2.5.2/topics/nmf
    # package nmfgpu4R and then algorithm AHCLS from https://arxiv.org/pdf/1407.7299.pdf
    # should be insanely fast


    # brunet works with KL divergence, lee with Euclidean distances -> lee is a factor 2 faster
    # decided to go for brunet as it results in a better factorisation from where the reconstruction of the original matrix shows less extreme errors.
    # https://cran.r-project.org/web/packages/NMF/vignettes/NMF-vignette.pdf
    # for an overview of the different NMF algorithms
    # nmfalg = "brunet"
    ranks = min_sig:max_sig

    # Orientation of the matrix matters dramatically for the speed! NMF algorithms were made for gene expression data. So many features (rows) and few samples (columns). 
    # "and X a matrix with n rows - the measured features - and p columns - the samples - with only non-negative entries" Gaujoux and Seoighe (2010)
    # However, we have genomics data with many samples (rows) and few features (columns). Hence we use the tall matrix. If the wide matrix is used, runtime and memory grows exponentially.
    mat = sample_by_component
    if(ncol(mat) > nrow(mat)) { mat=t(mat) }

    if(largeMatrix) {
        # for the nmf
        library(NMF)
        # for the manual parallelisation and also used by NMF for their internal parallelisation
        library(foreach)
        # needed for faster nmf runs based on the developer's own recommendations:
        # https://cran.r-project.org/web/packages/NMF/vignettes/NMF-vignette.pdf
        # page 16 ff.
        library(doParallel)
        # for non-Windows machines
        if( ! grepl( "windows", get_os(), ignore.case = TRUE) ) {
            library(bigmemory)
            library(synchronicity)
        }
        # for remodelling data structures
        library(reshape2)
        # for plotting
        library(ggplot2)

        lNMFs = listNMF(mat, ranks = ranks, seed = seed, 
            iter = iter, nmfalg = nmfalg, cores = cores)

        # get new rank of sigs when using "snmf/r"
        if(nmfalg == "snmf/r") {
            min_rank = ncol( basis( lNMFs[[1]] ) )
            max_rank = ncol( basis( lNMFs[[length(lNMFs)]] ) )
            ranks = min_rank:max_rank
        }

        # run random matrix only when asked for
        if( isTRUE(randomMat) ) {
             rmat = NMF::randomize(mat)
            lNMFsran = listNMF(rmat, ranks = ranks, seed = seed, 
                iter = iter, nmfalg = nmfalg, cores = cores)
            dfScores = getScores(lNMFs, lNMFsran, ranks, mat, randomMat)
        } else {
            lNMFsran=lNMFs
            dfScores = getScores(lNMFs, lNMFsran, ranks, mat, randomMat)
        }
        
        p = plotScores(dfScores)
        
        pdf(file=outfile, width=5, height=10 )
        print(p)
        dev.off()

        out = list( "nmf" = lNMFs, "nmfrandom" = lNMFsran, 
            "scores" = dfScores, "plot" = p, "randomMat" = randomMat )

        return(out)

   
    } else {
        estim.r <- NMF::nmfEstimateRank(mat, ranks,seed = seed,nrun=iter, method=nmfalg, .opt=list(paste0("vp", cores) ) )
        V.random <- NMF::randomize(mat)
        estim.r.random <- NMF::nmfEstimateRank(V.random, ranks, seed =seed,nrun=iter, method=nmfalg, .opt=list(paste0("vp", cores) ) )
    
        p<-NMF::plot(estim.r,estim.r.random, 
                what = c("cophenetic", "dispersion","sparseness", "silhouette"),
                xname="Observed",yname="Randomised",main="")

        pdf(file=outfile, width=5, height=10 )
        p
        dev.off()

        return(p)

    }

}

extractCopynumberFeatures<-function(CN_data, cores = 1, prePath="data/", allowedError = 0.1, rmNorm = FALSE)
{
    #get chromosome lengths
    chrlen<-read.table(paste0(prePath, "hg19.chrom.sizes.txt"),sep="\t",stringsAsFactors = F)[1:24,]
    
    #get centromere locations
    gaps<-read.table(paste0(prePath, "gap_hg19.txt"),sep="\t",header=F,stringsAsFactors = F)
    centromeres<-gaps[gaps[,8]=="centromere",]
    
    if(cores > 1) {
        require(foreach)
        require(doMC)

        registerDoMC(cores)

        temp_list = foreach(i=1:6) %dopar% {
            if(i == 1){
                list(segsize = getSegsize(CN_data, rmNorm = rmNorm) )
            } else if (i == 2) {
                list(bp10MB = getBPnum(CN_data,chrlen) )
            } else if (i == 3) {
                list(osCN = getOscilation(CN_data,chrlen) )
            } else if (i == 4) {
                list(bpchrarm = getCentromereDistCounts(CN_data,centromeres,chrlen) )
            } else if (i == 5) {
                list(changepoint = getChangepointCN(CN_data, allowedError, rmNorm = rmNorm) )
            } else {
                list(copynumber = getCN(CN_data, rmNorm = rmNorm) )
            }
        
        }
        
        # Another failsafe that the outcome is definitely numeric
        temp_list = unlist( temp_list, recursive = FALSE )
        outList = lapply(temp_list, function(thisDF) { 
            thisDF[,2] = as.numeric(thisDF[,2])
            return(thisDF)
        })
        return( outList )
        
    } else {  
        
        segsize<-getSegsize(CN_data)
        bp10MB<-getBPnum(CN_data,chrlen)
        osCN<-getOscilation(CN_data,chrlen)
        bpchrarm<-getCentromereDistCounts(CN_data,centromeres,chrlen)
        changepoint<-getChangepointCN(CN_data)
        copynumber<-getCN(CN_data)
        
        temp_list = list(segsize=segsize,bp10MB=bp10MB,osCN=osCN,bpchrarm=bpchrarm,changepoint=changepoint,copynumber=copynumber)
        temp_list = unlist( temp_list, recursive = FALSE )
        outList = lapply(temp_list, function(thisDF) { 
            thisDF[,2] = as.numeric(thisDF[,2])
            return(thisDF)
        })
        return( outList )
        
    }

}

fitMixtureModelsGM<-function(CN_features, seed=77777, min_comp=2, max_comp=10, min_prior=0.001, model_selection="BIC", nrep=1, niter=1000)
{

    dat<-as.numeric(CN_features[["segsize"]][,2])
    segsize_mm<-fitComponentGM(dat,seed=seed,model_selection=model_selection,
                               min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp)
    
    dat<-as.numeric(CN_features[["bp10MB"]][,2])
    bp10MB_mm<-fitComponentGM(dat,dist="pois",seed=seed,model_selection=model_selection,
                              min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp)
    
    dat<-as.numeric(CN_features[["osCN"]][,2])
    osCN_mm<-fitComponentGM(dat,dist="pois",seed=seed,model_selection=model_selection,
                            min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp)
    
    dat<-as.numeric(CN_features[["bpchrarm"]][,2])
    bpchrarm_mm<-fitComponentGM(dat,dist="pois",seed=seed,model_selection=model_selection,
                                min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp)
    
    dat<-as.numeric(CN_features[["changepoint"]][,2])
    changepoint_mm<-fitComponentGM(dat,seed=seed,model_selection=model_selection,
                                   min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp)
    
    dat<-as.numeric(CN_features[["copynumber"]][,2])
    copynumber_mm<-fitComponentGM(dat,seed=seed,model_selection=model_selection,
                                  nrep=nrep,min_comp=min_comp,max_comp=max_comp,min_prior=0.005,niter=2000)
    
    list(segsize=segsize_mm,bp10MB=bp10MB_mm,osCN=osCN_mm,bpchrarm=bpchrarm_mm,changepoint=changepoint_mm,copynumber=copynumber_mm)
    
}

fitMixtureModels<-function(CN_features, seed=77777, min_comp=2, max_comp=10, min_prior=0.001, model_selection="BIC",
                            nrep=1, niter=1000, cores = 1, featsToFit = seq(1, 6))
{

    if(cores > 1) {
        require(foreach)
        require(doMC)
        require(mclust)

        registerDoMC(cores)

        min_comp = getParams(min_comp, length(featsToFit))
        max_comp = getParams(max_comp, length(featsToFit))
        min_prior = getParams(min_prior, length(featsToFit))
        nrep = getParams(nrep, length(featsToFit))
        niter = getParams(niter, length(featsToFit))

        temp_list = foreach(i=1:6) %dopar% {

            if(i == 1 & i %in% featsToFit ){
            
                dat<-as.numeric(CN_features[["segsize"]][,2])
                list( segsize = fitComponent(dat,seed=seed,model_selection=model_selection,
                    min_prior=min_prior[i],niter=niter[i],nrep=nrep[i],min_comp=min_comp[i],max_comp=max_comp[i]) )
            
            } else if (i == 2 & i %in% featsToFit ) {
            
                dat<-as.numeric(CN_features[["bp10MB"]][,2])
                list( bp10MB = fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                    min_prior=min_prior[i],niter=niter[i],nrep=nrep[i],min_comp=min_comp[i],max_comp=max_comp[i]) )
            } else if (i == 3 & i %in% featsToFit ) {
            
                dat<-as.numeric(CN_features[["osCN"]][,2])
                list( osCN = fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                    min_prior=min_prior[i],niter=niter[i],nrep=nrep[i],min_comp=min_comp[i],max_comp=max_comp[i]) )
            
            } else if (i == 4 & i %in% featsToFit ) {
            
                dat<-as.numeric(CN_features[["bpchrarm"]][,2])
                list( bpchrarm = fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                    min_prior=min_prior[i],niter=niter[i],nrep=nrep[i],min_comp=min_comp[i],max_comp=max_comp[i]) )
            
            } else if (i == 5 & i %in% featsToFit ) {
            
                dat<-as.numeric(CN_features[["changepoint"]][,2])
                list( changepoint = fitComponent(dat,seed=seed,model_selection=model_selection,
                    min_prior=min_prior[i],niter=niter[i],nrep=nrep[i],min_comp=min_comp[i],max_comp=max_comp[i]) )

            } else if (i == 6 & i %in% featsToFit) {
            
                dat<-as.numeric(CN_features[["copynumber"]][,2])
                list( copynumber = fitComponent(dat,seed=seed,model_selection=model_selection,
                    min_prior=min_prior[i],niter=niter[i],nrep=nrep[i],min_comp=min_comp[i],max_comp=max_comp[i]) )

            }

        }
        unlist( temp_list, recursive = FALSE ) 
    } else {
        dat<-as.numeric(CN_features[["segsize"]][,2])
        segsize_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
                                 min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp)

        dat<-as.numeric(CN_features[["bp10MB"]][,2])
        bp10MB_mm<-fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                                min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp)

        dat<-as.numeric(CN_features[["osCN"]][,2])
        osCN_mm<-fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                              min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp)

        dat<-as.numeric(CN_features[["bpchrarm"]][,2])
        bpchrarm_mm<-fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                                  min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp)
        
        dat<-as.numeric(CN_features[["changepoint"]][,2])
        changepoint_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
                                     min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp)

        dat<-as.numeric(CN_features[["copynumber"]][,2])
        copynumber_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
                                nrep=nrep,min_comp=min_comp,max_comp=max_comp,min_prior=0.005,niter=2000)

        
        list(segsize=segsize_mm,bp10MB=bp10MB_mm,osCN=osCN_mm,bpchrarm=bpchrarm_mm,changepoint=changepoint_mm,copynumber=copynumber_mm)
    }
}

retrieveSumOfPosteriors = function(CN_feature, components, name) {

    if( class(components) == "flexmix" ) {
        curr_posterior = components@posterior$scaled
    } else if( class(components) == "Mclust" ) {
        curr_posterior = components$z
    } else {
        stop("Class of input is neither flexmix nor Mclust.")
    }

    # for the names only
    mat<-cbind(CN_feature,curr_posterior)
    # actual values not needed anymore
    mat = mat[,c(-2)]
    # get sum for each column for each sample
    posterior_sum = aggregate(. ~ ID, mat, sum)
    # get ID column as rownames and delete afterwards
    rownames(posterior_sum) = posterior_sum$ID
    posterior_sum$ID = NULL

    
    # sort by mean of mixture model.
    if( class(components) == "flexmix" ) {
        params<-flexmix::parameters(components)
    } else if( class(components) == "Mclust" ) {
        params = components$parameters$mean
    }
        
    if(!is.null(nrow(params)))
    {
        posterior_sum<-posterior_sum[,order(params[1,])]
    } else {
        posterior_sum<-posterior_sum[,order(params)]
    }
    # give meaningful names and return
    colnames(posterior_sum)<-paste0(name,1:ncol(posterior_sum))
    return( posterior_sum )
    
}

generateSampleByComponentMatrix2=function(CN_features, all_components, cores = 1, rowIter = 1000, 
                                          feat = "all", calcPost=FALSE) {

    if(is.null(all_components))
    {
        all_components<-readRDS("data/component_parameters.rds")
    }

    if( isFALSE(calcPost) ) {
        if(cores > 1){
            require(foreach)
            require(doMC)

            feats = c( "segsize", "bp10MB", "osCN", "bpchrarm", "changepoint", "copynumber" )
            registerDoMC(cores)

            full_mat = foreach(feat=feats, .combine=cbind) %dopar% {
                retrieveSumOfPosteriors(CN_features[[feat]],all_components[[feat]], 
                    feat)
            }
        } else {
            full_mat<-cbind(
                    retrieveSumOfPosteriors(CN_features[["segsize"]], all_components[["segsize"]],"segsize"),
                    retrieveSumOfPosteriors(CN_features[["bp10MB"]],all_components[["bp10MB"]],"bp10MB"),
                    retrieveSumOfPosteriors(CN_features[["osCN"]],all_components[["osCN"]],"osCN"),
                    retrieveSumOfPosteriors(CN_features[["bpchrarm"]],all_components[["bpchrarm"]],"bpchrarm"),
                    retrieveSumOfPosteriors(CN_features[["changepoint"]],all_components[["changepoint"]],"changepoint"),
                    retrieveSumOfPosteriors(CN_features[["copynumber"]],all_components[["copynumber"]],"copynumber"))          
        }
    } else {
        if(cores > 1){
            require(foreach)
            require(doMC)

            feats = c( "segsize", "bp10MB", "osCN", "bpchrarm", "changepoint", "copynumber" )
            registerDoMC(cores)

            full_mat = foreach(feat=feats, .combine=cbind) %dopar% {
                calculateSumOfPosteriors(CN_features[[feat]],all_components[[feat]], 
                    feat, rowIter = rowIter, cores = subcores)
            }
        } else {
        
            full_mat<-cbind(
            calculateSumOfPosteriors(CN_features[["segsize"]],all_components[["segsize"]],"segsize"),
            calculateSumOfPosteriors(CN_features[["bp10MB"]],all_components[["bp10MB"]],"bp10MB"),
            calculateSumOfPosteriors(CN_features[["osCN"]],all_components[["osCN"]],"osCN"),
            calculateSumOfPosteriors(CN_features[["bpchrarm"]],all_components[["bpchrarm"]],"bpchrarm"),
            calculateSumOfPosteriors(CN_features[["changepoint"]],all_components[["changepoint"]],"changepoint"),
            calculateSumOfPosteriors(CN_features[["copynumber"]],all_components[["copynumber"]],"copynumber"))
            
        }

    }

    rownames(full_mat)<-unique(CN_features[["segsize"]][,1])
    full_mat[is.na(full_mat)]<-0
    full_mat

}

generateSampleByComponentMatrix<-function(CN_features, all_components=NULL, cores = 1, rowIter = 1000, feat = "all")
{
    # Addendum 01/02/2019: Took out foreach loop and "subcores" in order to make it computationally less intensive.
    # The reasoning is that the splitting of each feature (done by "rowIter") is more important than the parallelisation between features.
    # With 3 cores this will run under 10 minutes for the whole TCGA data set and needs less than 10GB memory.
    if(is.null(all_components))
    {
        all_components<-readRDS("data/component_parameters.rds")
    }

    if ( feat == "all" ) {
        full_mat<-cbind(
            calculateSumOfPosteriors(CN_features[["segsize"]],all_components[["segsize"]],"segsize", rowIter = rowIter, cores = cores),
            calculateSumOfPosteriors(CN_features[["bp10MB"]],all_components[["bp10MB"]],"bp10MB", rowIter = rowIter, cores = cores),
            calculateSumOfPosteriors(CN_features[["osCN"]],all_components[["osCN"]],"osCN", rowIter = rowIter, cores = cores),
            calculateSumOfPosteriors(CN_features[["bpchrarm"]],all_components[["bpchrarm"]],"bpchrarm", rowIter = rowIter, cores = cores),
            calculateSumOfPosteriors(CN_features[["changepoint"]],all_components[["changepoint"]],"changepoint", rowIter = rowIter, cores = cores),
            calculateSumOfPosteriors(CN_features[["copynumber"]],all_components[["copynumber"]],"copynumber", rowIter = rowIter, cores = cores) )
    
        rownames(full_mat)<-unique(CN_features[["segsize"]][,1])
        full_mat[is.na(full_mat)]<-0
        full_mat

    } else {
        # One feature only
        full_mat = calculateSumOfPosteriors(CN_features[[feat]], all_components[[feat]], feat, rowIter = rowIter, cores = subcores)

        rownames(full_mat)<-unique(CN_features[[feat]][,1])
        full_mat[is.na(full_mat)]<-0
        full_mat

    }

}

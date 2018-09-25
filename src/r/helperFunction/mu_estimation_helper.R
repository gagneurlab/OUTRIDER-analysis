runMuBenchmark <- function(funs, ods, odsls, q, idx=seq_along(funs), BPPARAM=bpparam()){
    if(missing(odsls)){
        odsls <- list()
        odsls[names(funs)] <- list(NULL)
    }
    param <- NULL
    if(length(idx) < 2){
        param <- BPPARAM
        BPPARAM <- SerialParam()
    }
    system.time(
        tmpls <- bplapply(names(funs)[idx], ods=ods, funs=funs, param=param,
                        BPPARAM=BPPARAM, FUN=function(x, ods, funs, param) {
                if(is.null(param)){
                    param <- SerialParam()
                }
                if(is.null(funs[[x]])){
                    o <- autoCorrect(ods, q=q, autoCorrect=TRUE, implementation=x,
                            BPPARAM=param, debug=TRUE)
                } else {
                    o <- funs[[x]](ods, BPPARAM=param)   
                }
                if(!grepl('^ed', x)){
                    o <- fit(o, BPPARAM=param)
                }
                o <- computePvalues(o, BPPARAM=param)
                o <- computeZscores(o, peerResidual=grepl('peer', x, ignore.case=TRUE))
        })
    )
    odsls[idx] <- tmpls
    names(odsls) <- names(funs)
    return(odsls)
}

getRecallPrecisionTable <- function(odsls, index, padj=0.05){
    tTP <- sum(index != 0)
    resAll <- t(sapply(odsls, function(x){
        p <- sum(aberrant(x, padjCutoff=padj))
        tp <- sum(aberrant(x, padjCutoff=padj)[index != 0])
        return(list(tp=tp, p=p, precision=round(tp/p, 3), 
                    recall=round(tp/tTP, 3), padj=padj))
    }))
    dn <- dimnames(resAll)
    resAll <- data.table(
            matrix(unlist(resAll), ncol=length(dn[[2]]), dimnames=dn), 
            keep.rownames='Correction')
    
    return(resAll)
}


prepareMuErrorPlotData <- function(odsls, simData, met2plot, it){
    noIterNames <- grep('^((peer|pca|AE)\\d?)|(ed.*)$', perl=TRUE, value=TRUE, ignore.case=TRUE, names(odsls))
    geneIndex <- matrix(rowAnys(simData$index!=0),
            ncol=ncol(simData$index), nrow=nrow(simData$index), byrow=FALSE)
    mu <- simData$mu
    s_est <- DESeq2::estimateSizeFactorsForMatrix(k)
    trueTheta <- simData$theta
    
    dt <- data.table(
        count       = c(simData$k), 
        mu          = c(mu), 
        sdVec       = simData$sdVec,
        trueTheta   = simData$theta,
        Outlier     = c(simData$index), 
        OutlierGene = c(geneIndex))
    dtdisp <- copy(dt)
    
    for(i in noIterNames){
        dt[[i]]     <- c(normalizationFactors(odsls[[i]]))
        dtdisp[[i]] <- c(dispersions(odsls[[i]]))
    }
    
    #'
    #' Add methods with iterations
    for(midx in met2plot){
        for(i in it){
            m <- names(odsls)[midx]
            tmpods <- odsls[[m]]
            if(i == 'last'){
                i <- max(as.numeric(gsub('iter_', '', 
                        grep('iter', names(metadata(tmpods)), value=TRUE))))
            }
            name <- paste0("iter_", i)
            if(!name %in% names(metadata(tmpods))){
                next
            }
            dname <- paste0(m, '_', gsub('_([0-9])$', '_0\\1', name))
            w <- getAEData(tmpods, w=metadata(tmpods)[[name]]$w)$norm
            disp <- dispersions(odsls[[m]])
            if(is.null(disp)){
                disp <- NA_real_
            }
            dt[[dname]] <- c(w)
            dtdisp[[dname]] <- c(disp)
        }
    }
    
    
    #' 
    #' Compute statistics
    mnames <- colnames(dt)[!names(dt) %in% c('Outlier', 'OutlierGene', 'count', 'mu', 'sdVec', 'bins', 'trueTheta')]
    dt <-melt(dt, measure.vars=mnames, variable.name='Method', value.name='mu_est')
    dtdisp <-melt(dtdisp, measure.vars=mnames, variable.name='Method', value.name='disp_est')
    dt$disp_est <- dtdisp$disp_est
    dt[,logError:=log2(mu)-log2(mu_est)]
    dt[,logError2:=logError^2]
    dt[,logThetaError:=log2(trueTheta)-log2(disp_est)]
    dt[,logThetaError2:=logThetaError^2]
    dt[,logDevK:=log2(count+1)-log2(mu_est)]
    dt[,logDevK2:=logDevK^2]
    dt[,bins:=cut(log10(mu), c(min(log10(mu))-1, quantile(log10(mu), p=1:10/10)))]
    dt[,thetaBins:=cut(trueTheta, c(-0, quantile(trueTheta, p=1:10/10), Inf))]
    dt[,thetaEstBins:=cut(disp_est, c(-0, quantile(disp_est, p=1:10/10), Inf))]
    
    return(dt)
}


muErrWOS2 <- function(x, mu){
    rmerr <- rowMeans(normalizationFactors(x) / mu)
    (normalizationFactors(x) - mu*rmerr)^2
}

plotMuErrorWOSysErr <- function(mu, odsls){
    errls <- lapply(names(odsls), ols=odsls, mu=mu, function(n, ols, mu){
        data.table(error=abs(as.vector(muErrWOS2(ols[[n]], mu))), 
                   mu=as.vector(mu), Method=n)
    })
    dterr <- rbindlist(errls)
    
    ggplot(dterr, aes(mu, error)) + geom_hex() + 
        grids() + facet_wrap('Method') + 
        scale_x_log10() + 
        geom_abline(intercept = 0, slope = 0, col='red')
}


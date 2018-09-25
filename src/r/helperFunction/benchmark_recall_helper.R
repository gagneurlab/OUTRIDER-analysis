
pvalue_score_sign <- function(pvalue, zscore, inj){
    score <- pvalue
    score[inj == 'high' & zscore < 0] <- 1
    score[inj == 'low' & zscore > 0] <- 1
    score
}

readZscoreCutoff <- function(f, cutoffs=c(2,3)){
    dt <- fread(f)
    dt[,score:=score_sign(zScore, inj)]
    ans <- lapply(cutoffs, function(x){
        ans <- dt[score > x]
        if(nrow(ans) == 0){
            return(dt[,.(recall=0, precision=0, Nsamples, inj, inj_value, 
                    correction, pmethod, cutoff=x)][1])
        }
        ans[, .(
            recall=sum(trueOutlier != 0)/unique(numInjected),
            precision=sum(trueOutlier != 0)/.N,
            Nsamples = unique(Nsamples),
            inj = unique(inj),
            inj_value = unique(inj_value),
            correction = unique(correction),
            pmethod = unique(pmethod),
            cutoff=x)]
    })
    rbindlist(ans)
}

getBootstrapScore <- function(dt, pval, maxRows){
    if(isTRUE(pval)){
        dt[, score:= -pvalue_score_sign(pvalues, zScore, inj)]
    } else {
        dt[, score:=score_sign(zScore, inj)]
    }
    dt <- dt[order(score, decreasing=TRUE)][label == TRUE | 1:.N <= maxRows]
    dt
}

#
# gather data
#
readBootstrapData <- function(f, kzcut=Inf, pvalue=FALSE, zscoreForAll=FALSE, 
                    byFile=FALSE, maxRows=1e5, methodOnly=FALSE){
    if(is.data.table(f)){
         dt <- f
         f <- paste0('.../', unique(dt$correction), '/...')
    } else {
        dt <- fread(f)
    }
    dt[, label:=trueOutlier != 0]
    if(isTRUE(byFile)){
        pvalue <- !grepl('/correction=(pca|peer)/', f)
    }
    if(grepl('/correction=(baseCooks|basePearsonRes)/', f)){
        pvalue <- FALSE
    }
    
    goodEvents <- dt[,abs(kzScore) < kzcut]
    print(table(goodEvents))
    dt <- dt[goodEvents]
    
    dt1 <- getBootstrapScore(dt, !zscoreForAll & pvalue, maxRows)
    dt1 <- dt1[,                    bootstrapPR(score, label, total=max(numInjected), n=200, mc.cores=5)]
    dt2 <- getBootstrapScore(dt, pvalue, maxRows)
    dt2 <- dt2[padjust < FDR_LIMIT, bootstrapPR(score, label, total=max(numInjected), n=200, mc.cores=5)]
    
    dt1 <- rbind(dt1, data.table(score=-Inf, label=FALSE, rank=Inf,
            TP=max(dt1$TP), precision=0, recall=1, lower=0, upper=0))
    dt1 <- rbind(dt1, dt1[min(rank)==rank][1][,.(
            score, label, rank=0, TP, precision, recall=0, lower, upper)])
    if(nrow(dt2[!is.na(score)]) > 0){
        dt2 <- rbind(dt2, dt2[min(rank)==rank][1][,.(
                score, label, rank=0, TP, precision, recall=0, lower, upper)])
    }
    
    if(isTRUE(pvalue)){
        if(isTRUE(zscoreForAll)){
            dt1[, type:='Z score rank']
        } else {
            dt1[, type:='P-value rank']
        }
        dt2[, type:=paste('P-value rank \n& FDR <', FDR_LIMIT)]
    } else {
        dt1[, type:='Z score rank']
        dt2[, type:=paste('Z score rank \n& FDR <', FDR_LIMIT)]
    }
    if(isTRUE(methodOnly)){
        if(grepl('/correction=(pca|peer)/', f)){
            dt2 <- dt2[0]
        } else {
            dt1 <- dt1[0]
        }
    }
    
    dt2 <- dt2[!is.na(score)]
    ans <- rbind(
        dt1[,.(type, rank, TP, precision, recall, lower, upper)],
        dt2[,.(type, rank, TP, precision, recall, lower, upper)])
    ans$inj_value  <- unique(dt$inj_value)
    ans$inj        <- unique(dt$inj)
    ans$correction <- unique(dt$correction)
    ans
}

plotRibbonBenchmark <- function(data, title, linetype=c(1,2), maxRows=1e5,
                    wrap_function=function() facet_grid(inj_value ~ inj), 
                    zscoreData=NULL){
    #' Reduce plotting points
    prProbs <- c(1000/maxRows, 1-1000/maxRows)
    data <- data[rank < 1000 | sample(x=c(TRUE, FALSE), size=.N, replace=TRUE, prob=prProbs)]
    #data <- data[order(type, correction, inj_value, inj, precision, recall)]
    #dup <- duplicated(data[,.(type, correction, inj_value, inj, precision, recall)])
    #data <- data[!dup]
    print(head(data))
    
    #' Plot it
    gg <- ggplot(data, aes(recall, precision, color=correction, linetype=type)) +
        geom_line() +
        geom_ribbon(alpha=0.2, col=NA, aes(x=recall, fill=correction, 
                                           linetype=type,
                                           ymin = pmax(0, pmin(1, lower)), 
                                           ymax = pmax(0, pmin(1, upper)))) +
        scale_x_continuous(breaks = c(0,0.5,1)) +
        scale_y_continuous(breaks = c(0,0.5,1)) + 
        scale_linetype_manual(values=linetype) + 
        wrap_function() +
        grids(linetype="dashed") + 
        theme(strip.text.y = element_text(margin = margin(0,0.2,0,0.2, "cm")),
              strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"))) +
        theme(legend.key.height=unit(0.8, "cm")) + 
        labs(x='Recall', y='Precision', color="Correction", 
             fill='Correction', linetype="Method")
    if(!missing(title)){
        gg <- gg + ggtitle(title)
    }
    if(!is.null(zscoreData)){
        gg <- gg + 
            geom_point(data=zscoreData, 
                    aes(recall, precision, shape=factor(cutoff)), size=3) + 
            scale_shape_manual(values=c(1,2)) + 
            labs(shape='Z score cutoff')
    }
    gg
}


getBenchMarkDT <- function(ddir, regexs, name){
    # get files
    files <- list.files(path=file.path(ddir, "RecallAnalysis"), 
            pattern = '*_plot.tsv', recursive=TRUE, full.names=TRUE)
    
    # filter
    files <- files[rowAlls(sapply(regexs, grepl, x=files))]
    
    # get data
    dt <- bplapply(files, function(x) {
        a <- fread(x, stringsAsFactors=TRUE, nrows=1E5);
        a[,file:=x];
        a[,name:=gsub("_plot.tsv", "", name)]
        a
    }) %>% rbindlist(fill=TRUE)
    dt
    
    return(list(dt=dt, files=files))
}

getInjVarNum <- function(dt){
    files <- unique(dt[,file])
    dataset <- gsub("_plot.tsv", "", unique(basename(files)))
    
    injectFiles <- list.files(unique(dirname(dirname(dirname(dirname(files))))),
                              pattern=paste0(dataset, "_injectedcounts.RDS"), full.names=TRUE)
    numInjected <- sapply(injectFiles, function(x) {
        sum(assay(readRDS(x), "inj_mask") != 0)})
    
    dt <- data.table(
        numInjected=numInjected,
        inj=gsub("inj=", "", basename(dirname(dirname(injectFiles)))),
        inj_value=gsub("zscore=", "", basename(dirname(injectFiles)))
    )
    return(dt)
}

score_sign <- function(zScore, inj){
    score <- zScore
    score[which(inj=='low')] <- -score[which(inj=='low')]
    score[which(inj=='both')] <- abs(score[which(inj=='both')])
    return(score)
}

plotRecall <- function(dt, correctionMethods, kzCutoff=3, 
                    inj_values=c('zscore=2', 'zscore=3', 'zscore=4', 'zscore=6'),
                    filePrefix="plotBenchmark", numInjectedDt=NULL){
    dt1 <- dt[correction %in% correctionMethods]
    dt1 <- dt1[abs(kzScore) < kzCutoff]
    dt1[,Correction:=paste(name, Nsamples, q, correction, sep="_")]
    dt1[,score:= score_sign(zScore, inj)]
    dt1[,sig_score:= ifelse(padjust<FDR_LIMIT & score>0, score, -Inf)]
    
    dt1[,label:=as.logical(abs(trueOutlier))]
    
    PR <- dt1[,c(computePR(label, score), Method='Z-score'), 
              by=c('Correction', 'inj_value', 'inj')]
    PR2 <- dt1[,c(computePR(label, sig_score), Method=paste0('P-adj<', FDR_LIMIT)), 
               by=c('Correction', 'inj_value', 'inj')]
    dt2plot <- rbind(PR, PR2)
    
    # subset using diff==True for less accuarte but faster plotting
    dt2plot <- dt2plot[inj_value %in% inj_values]
    dt2plot[!duplicated(recall, fromLast=TRUE), diff:=TRUE, by='inj_value,Method']
    dt2plot[rank == 1, diff:=TRUE]
    
    # change wording
    dt2plot[!is.na(Cutoff), Cutoff:=gsub('Zscore>', '|Z| > ', Cutoff)]
    dt2plot[, inj_value:=as.factor(gsub('zscore=', 'Simulated |Z| = ', as.character(inj_value)))]
    dt2plot[,Method:=gsub('^(Z-score)?', '|Z| ranked', Method)]
    dt2plot[,Method:=gsub('P-adj<', '\n\t& FDR < ', Method)]
    dt2plot[,Correction:=gsub('autoCorrect', 'Autoencoder', Correction)]
    levels(dt2plot$inj) <- toCamelCase(levels(dt2plot$inj), capitalize=TRUE)
    
    if(is.null(numInjectedDt)){
        numInjectedDt <- getInjVarNum(dt1)
    }
    
    ggpr <- ggplot(dt2plot[diff==TRUE], aes(recall, precision, color=Correction, linetype=Method)) +
        geom_line() +
        geom_point(data=dt2plot[!is.na(Cutoff) & grepl('ranked$', Method)], aes(shape=Cutoff), size=2) +
        scale_x_continuous(breaks = c(0,0.5,1)) +
        scale_y_continuous(breaks = c(0,0.5,1)) + 
        scale_linetype_manual(values=c(2,1)) + 
        facet_grid(inj_value ~ inj) +
        grids(linetype="dashed") +
        theme(strip.text.y = element_text(margin = margin(0,0.2,0,0.2, "cm")),
              strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"))) +
        theme(legend.key.height=unit(0.8, "cm")) + 
        labs(title=paste("Precision-Recall:", paste(unique(dt1[,name]), collapse=" ")), x='Recall', y='Precision', 
             color="Correction", linetype="Method", shape="Cutoff")
    ggpr
    
    ggrr <- ggplot(dt2plot[rank < 1000], aes(rank, recall, color=Correction, linetype=Method)) +
        geom_line() +
        scale_linetype_manual(values=c(2,1)) + 
        geom_abline(data=numInjectedDt, intercept=0, slope=1/numInjectedDt$numInjected, linetype="dotted") + 
        facet_grid(inj_value ~ inj) +
        grids(linetype="dashed") +
        theme(strip.text.y = element_text(margin = margin(0,0.2,0,0.2, "cm")),
              strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"))) +
        theme(legend.key.height=unit(0.8, "cm")) +
        labs(title=paste("Recall-Rank:", paste(unique(dt1[,name]), collapse=" ")), x='Rank', y='Recall')
    ggrr
    
    ggmerged <- plot_grid(ggpr, ggrr)
    
    ggsave(filename=paste(filePrefix, '-merged.pdf'), plot=ggmerged,
           width=40, height=15, scale=0.7, units='in', device='pdf')
    
    #ggsave(filename=file.path(FIGDIR, paste0(FIGURE_NAME, '.png')),
    #       plot=ggrr, width=12, height=8, scale=0.7, units='in', device='png')
    
    return(list(ggmerged=ggmerged, ggpr=ggpr, ggrr=ggrr, dt2plot=dt2plot))
}


createEvalTable <- function(ods, q, correction, anMask='inj_mask', 
                    anTrueCts='trueCounts', Nsamples='All', inj='both',
                    inj_value='zscore=3', pmethod='BY', dataset='dataset',
                    nMeanBins=1){
    if(correction == "baseCooks"){
        Zscores <- assay(ods, 'cooks') * (2 * (counts(ods) - rowMeans(counts(ods)) > 0) - 1)
    } else if(correction == "basePearsonRes"){
        Zscores <- assay(ods, 'pearsonRes')
    } else {
        Zscores <- assay(ods, 'zScore')
    }
    
    Zscore_rank <- rank(-abs(Zscores))  # neg (-) for descending order
    trueOutlier <- assay(ods, anMask)
    totalTrueOut <- sum(abs(trueOutlier))
    padjust <- assays(ods)[['padjust']]
    
    # raw count zscore
    lk <- log2(assay(ods, anTrueCts) + 1)
    kzmat <- (lk - rowMeans(lk)) / rowSds(lk)
    
    # raw gene expression bin
    nMeanBins <- max(1, nMeanBins)
    rMeans <- rowMeans(counts(ods))
    quantilesBreaks <- quantile(rMeans, p=(1:nMeanBins/nMeanBins))
    meanBins <- cut(rMeans, breaks=c(0, quantilesBreaks), right=TRUE)
    
    ans <- data.table(
        numInjected = sum(trueOutlier != 0),
        trueOutlier = c(trueOutlier),
        pvalues = c(pValue(ods)),
        padjust = c(padjust),
        zScore = c(Zscores),
        zScore_rank = c(Zscore_rank),
        kzScore = c(kzmat),
        Nsamples=Nsamples,
        inj=inj, 
        inj_value=inj_value,
        correction=correction,
        pmethod=pmethod,
        dataset=dataset,
        q=q,
        geneMean=rMeans,
        geneMeanBin=meanBins)
    
    if('Ncov' %in% names(mcols(ods))){
        ans[, Ncovs := mcols(ods)[['Ncov']]]
    }
    
    setorder(ans, zScore_rank)
    ans[,recall := cumsum(abs(trueOutlier))/totalTrueOut]
    
    return(ans)
}

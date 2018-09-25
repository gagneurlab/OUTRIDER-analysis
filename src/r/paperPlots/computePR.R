computePR <- function(label, score){
    dt <- data.table(label, score)
    dt <- dt[order(-score)]
    dt[,rank := 1:length(score)]
    dt[,recall := cumsum(label)/sum(label)]
    dt[,precision := cumsum(label)/rank]
    
    dt <- dt[is.finite(score)]
    dt[,diff := c(diff(precision),1)>1E-10]
    #dt[1,diff:=TRUE] #always plot first point.
    
    # add common Z-score (2,3) cut-offs to table
    for(z in 2:3){
        dt[max(which(dt[,score]>z)), Cutoff:=paste0('Zscore>',z)]
    }
    
    return(dt)
}

getPRPlotData <- function(ods, trueOutlier, nBins=1){
    if(!is.data.table(ods)){
        # Get bins
        nBins <- max(1, nBins)
        rMeans <- rowMeans(counts(ods))
        mcols(ods)[['MeanBin']] <- cut(rMeans, 
                breaks=c(0,quantile(rMeans, p=(1:nBins/nBins))), right=TRUE)
        
        # Create data
        dt <- data.table(log_10_pvalue=-log10(c(assay(ods, 'pValue'))), 
                padjust=c(assay(ods, 'padjust')), 
                zScore=c(assay(ods, 'zScore')), 
                trueOutlier = unlist(trueOutlier), 
                bin=mcols(ods)[['MeanBin']])
    } else {
        dt <- ods
        dt[,bin:=1]
    }
    
    #' Compute Recall curves
    dt[,neg_zscore:= -zScore]
    dt[,sig_neg_zscore:= ifelse(padjust<FDR_LIMIT & neg_zscore>0, neg_zscore, -Inf)]
    
    #' Compute recall curves for individaul bins.
    dt2plot <- rbind(
            dt[,c(computePR(as.logical(abs(trueOutlier)), neg_zscore),     Method='Z-score'),                   by='bin'],
            dt[,c(computePR(as.logical(abs(trueOutlier)), sig_neg_zscore), Method=paste0('P-adj<', FDR_LIMIT)), by='bin']
    )
    
    # Naming
    dt2plot[!is.na(Cutoff), Cutoff:=gsub('Zscore>', '|Z| > ', Cutoff)]
    dt2plot[,Method:=gsub('^(Z-score)?', '|Z| ranked', Method)]
    dt2plot[,Method:=gsub('P-adj<', '\n\t& FDR < ', Method)]
    #levels(dt2plot$inj) <- toCamelCase(levels(dt2plot$inj), capitalize=TRUE)
    
    return(dt2plot)
}



plotSequencingCoverage <- function(ods, name){
    if(is.data.table(ods)){
        dt <- ods
    } else {
        cov <- colSums(counts(ods, normalized=FALSE))
        CovOds <- sort(cov)
        dt <- data.table(CovOds=CovOds, rank=seq_along(CovOds))
    }
    
    gg <- ggplot(dt, aes(rank,CovOds)) + geom_point() +
        xlab('Sample rank') + ylab('Sequencing depth') +
        ggtitle(paste0('Total sequencing depth (', name, ')')) +
        geom_hline(aes(yintercept=mean(CovOds) - 3*sd(CovOds), color = 'firebrick')) + 
        annotate('text', x=length(CovOds)*0.93, y=mean(CovOds) - 3*sd(CovOds), 
                label='Z-score = -3', vjust=-1, size=5) + 
        theme(legend.position = "none") + 
        grids(linetype="dashed") +
        scale_y_log10(
                breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x)))
    
    return(gg)
}


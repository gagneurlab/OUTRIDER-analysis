
plotPEERAlpha <- function(ods, plotLines=FALSE){
    require(peer)
    model <- metadata(ods)[['PEER_model']]
    data <- data.table(alpha=(model$alpha)[,1])
    data <- data[order(alpha, decreasing=TRUE)][,.(alpha, rank=seq_len(.N))]
    explainedVar <- data[,cut(cumsum(alpha), 
            c(0, sum(alpha)*0.8, sum(alpha)*0.9, sum(alpha)+1), 1:3)]
    explainedVar <- data.table(
            nFactor=sapply(1:2, function(x) max(which(x==explainedVar))))[,
            .(nFactor, label=paste0(c("80% explained", "90% explained"), "(n=", nFactor, ")"))]
    
    # plot it
    gg <- ggplot(data=data, aes(rank, alpha)) +
            geom_point() + geom_smooth(method="loess") + grids() +
            ggtitle("PEER: Inverse variance of factor weights") + 
            labs(x="Ranked factors", y="Inverse variance\nof factor weights")
    
    if(isTRUE(plotLines)){
        gg <- gg + 
            geom_vline(xintercept=explainedVar[,nFactor], show.legend = TRUE) +
            geom_text(data=explainedVar, aes(nFactor-0.5), 
                    label=explainedVar$label, y=max(data$alpha)*0.9, angle=90)
    }
    gg
}


plotZscoreInjSearch <- function(dt, digits=3){
    if(is(dt, "OutriderDataSet")){
        dt <- metadata(dt)[['optimalZscoreEncDim']]
    }
    dt <- dt[,.(enc=encodingDimension, z=zScore, loss=round(evaluationLoss, digits))]
    dt <- dt[order(enc, z)]
    dt[,bestEnc:=enc==getBestQDT(data.table(encodingDimension=enc, evaluationLoss=loss), digits=digits), by=z]
    
    # gg1 <- ggplot(dt, aes(factor(enc), factor(z), fill=loss)) + 
    #     geom_tile() + 
    #     scale_fill_gradient(low="black", high="white", breaks=200)
    
    gg <- ggplot(dt, aes(enc, loss, col=factor(z))) + 
        geom_line() + 
        geom_point(aes(enc, loss, size=bestEnc), pch=16) +
        scale_x_log10() + 
        grids(linetype='dotted') +
        scale_size_manual(values=c(1, 4)) + 
        labs(y='AUC: precision recall', x='Encoding dimension', 
                col='Z-score', size='Best q')
        
    #gg <- plot_grid(gg1, gg2, nrow = 2)
    
    return(gg)
}

#'---
#' title: Dataset summary
#' author: Christian Mertes
#' wb:
#'  input: 
#'   - odsList: '`sm expand(config["DATADIR"] + "/fitOutrider/{method}/{{dataset}}_ODS.RDS", method=list(set([config["AE_IMPLEMENTATION"]] + config["correction"])))`'
#'  output:
#'   - ggplots: '`sm config["htmlOutputPath"] + "/datasetSummary/{dataset}_ggplots.RDS"`'
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/datasetSummary/{dataset}.html"`'
#'  type: noindex
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---

if(FALSE){
    snakemake <- readRDS('tmp.snakemake.RDS'); snakemake 
}

source("./src/r/config.R")
library(plotly)

#' 
#' # Read datasets
#' 
ggplotFile <- snakemake@output$ggplots
odsFiles <- snakemake@input$odsList
dataset <- snakemake@wildcards$dataset
mc.cores <- min(10, length(odsFiles))
odsFiles
dataset

odsls <- mclapply(odsFiles, readRDS, mc.cores=mc.cores)
names(odsls) <- basename(dirname(odsFiles))

#' 
#' # Global QQ plot
#' 
dt <- rbindlist(mclapply(names(odsls), mc.cores=mc.cores, function(x){
    ans <- data.table(Method=x, pValue=c(pValue(odsls[[x]])))
    ans[is.na(pValue), pValue:=1]
    ans[pValue == 0, pValue:=min(pValue)/100]
    ans <- ans[order(pValue)][,.(Method, pValue, 
            nlo=-log10(pValue), nle=-log10(ppoints(.N)))]
}))

#' * sample for faster plotting
plotIDs <- dt[Method == unique(Method)[1],
        nle > 3 |
        nle > 2 & nle <= 3 & sample(c(TRUE, FALSE), .N, prob=c(0.01,  0.99), replace=TRUE) |
        nle > 0.5 & nle <= 2 & sample(c(TRUE, FALSE), .N, prob=c(0.0001, 0.9999), replace=TRUE) |
        nle <= 0.5 & sample(c(TRUE, FALSE), .N, prob=c(0.00001, 0.99999), replace=TRUE)]
dtPlot <- dt[rep(plotIDs, length(unique(Method)))]
gg <- ggplot(data=dtPlot, aes(nle, nlo, col=Method)) + geom_point() + 
    grids() +
    labs(title = paste('Global Q-Q plot for', dataset),  
         x=expression(paste(-log[10], " (expected ", italic(P), "-value)")),
         y=expression(paste(-log[10], " (obs. ", italic(P), "-value)"))) +
    geom_abline(intercept = 0, slope = 1)
    

#' * QQ plot + zoom
#+ fig.width=9, fig.height=9
gg 
gg + coord_cartesian(xlim=c(2,4), ylim=c(2,7))

#'
#' # Number of aberrant events per sample
#' 
getAberrantPlottingData <- function(ols, dataset='Kremer'){
    odls <- ols
    dt <- rbindlist(lapply(names(odls), function(x){
        dt <- data.table(Naberrant=aberrant(odls[[x]], by='sample'), Method=x)
        dt[, Sample_rank:=rank(Naberrant, ties.method = 'first')]
        dt[, medianNaberrant:=median(Naberrant)]
        dt[, nrow:=nrow(odls[[x]])]
        dt
    }))
    dt
}

plotAberrantFigure <- function(dt, dataset){
    cutoff <- max(dt$nrow) * 0.005
    ggp <- ggplot(dt, aes(Sample_rank, Naberrant, col=Method)) + 
        geom_line() +
        scale_y_log10() +
        geom_hline(linetype = "dotted", aes(yintercept=medianNaberrant, col=Method)) +
        geom_hline(linetype = "dashed", yintercept = cutoff, col='firebrick') + 
        labs(title=paste0('Aberrant genes per sample (', dataset, ')'),
             x='Sample rank',
             y='#Aberrantly\nexpressedgenes') + 
        geom_text(inherit.aes = FALSE, hjust='left', fontface='plain',
                  aes(x=2, y=cutoff*1.4, label='Aberrent samples')) + 
        geom_text(inherit.aes = FALSE, hjust='left', fontface='plain',
                  aes(x=2, y=max(dt$medianNaberrant)*1.4, label='Median'))
    ggp
}
plotDT <- getAberrantPlottingData(odsls, dataset)
ggAPS <- plotAberrantFigure(plotDT, dataset)
ggAPS

#' 
#' # Save ggplot object into RDS file
ggplotFile
saveRDS(list(gg, ggAPS), ggplotFile)

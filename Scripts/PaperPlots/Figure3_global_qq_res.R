#'---
#' title: Figure 2 Heatmap of GTEx data before/after correction
#' author: Christian Mertes
#' 
#' wb:
#'   input: 
#'     - odsg: '`sm expand(config["DATADIR"] + "/fitOutrider/{method}/Skin_Not_Sun_Exposed_Suprapubic_ODS.RDS", method=PAPER_METHODS)`'
#'     - odsk: '`sm expand(config["DATADIR"] + "/fitOutrider/{method}/Kremer_ODS.RDS", method=PAPER_METHODS)`'
#'   output:
#'     - pdf: '`sm config["FIGDIR"] + "/Figure3_global_result.pdf"`'
#'     - png: '`sm config["FIGDIR"] + "/Figure3_global_result.png"`'
#' output:
#'   html_document:
#'     code_folding: show
#'     code_download: TRUE
#'---

if(FALSE){
    snakemake <- readRDS('tmp.snakemake.RDS'); snakemake 
    getFiles <- function(dName){
        DIR <- '/data/ouga04b/ag_gagneur/project_local/scared/paper/revision/run0910-final/data/'
        paste0(DIR, "fitOutrider/", 
                c(CONFIG_YAML$FIGURE_IMPLEMENTATION, 'peer', 'pca'), "/", dName, "_ODS.RDS")
    }
    odsFiles <- list(
        GTEx =   getFiles("Skin_Not_Sun_Exposed_Suprapubic"),
        Kremer = getFiles("Kremer"))
    pdfOut <- "tmp.pdf"
}


#' 
#' # Read input and config
#' 
source('./src/r/config.R')

odsFiles <- list(
    GTEx = snakemake@input$odsg,
    Kremer = snakemake@input$odsk
)
pdfOut <- snakemake@output$pdf
pngOut <- snakemake@output$png

odsFiles
pdfOut
pngOut

#' 
#' # Read data
#' 
odsls <- lapply(odsFiles, function(x){
    ans <- lapply(x, readRDS)
    names(ans) <- basename(dirname(unlist(x)))
    ans
})

#' # Plotting functions
getQQPlottingData <- function(ols, dataset='Kremer', conf.alpha=0.05){
    odls <- ols[[dataset]]
    lk <- log2(1+counts(odls[[1]]))
    z <- (lk - rowMeans(lk))/rowSds(lk)
    qqPlotDT <- rbindlist(lapply(names(odls), function(x){
        o <- odls[[x]]
        dt <- data.table(
                z=c(z), 
                observedPvalue=c(assay(o, 'pValue')), 
                aberrant=c(aberrant(o)), Method=x)
        dt[order(observedPvalue)]
    }))
    
    # qqPlotDT <- qqPlotDT[abs(z)<2]
    # qqPlotDT <- qqPlotDT[aberrant!=TRUE]
    qqPlotDT[,expectedPvalue:= ppoints(observedPvalue), by=Method]
    
    # set confidence
    qqPlotDT[,nlupper:=-log10(qbeta(  conf.alpha/2, 1:.N, .N:1)), by=Method]
    qqPlotDT[,nllower:=-log10(qbeta(1-conf.alpha/2, 1:.N, .N:1)), by=Method]
    
    ## sample to avoid plotting problems.
    qqPlotDTSampled <- qqPlotDT[
            observedPvalue <  1E-3 |
            observedPvalue <  1E-2 & sample(c(TRUE, FALSE), nrow(qqPlotDT), prob = c(0.1,  0.9),  replace = TRUE)|
            observedPvalue >= 1E-2 & sample(c(TRUE, FALSE), nrow(qqPlotDT), prob = c(0.01, 0.99), replace = TRUE)]
    qqPlotDTSampled[,neglog10expectedPvalue := -log10(expectedPvalue)]
    qqPlotDTSampled[,neglog10observedPvalue := -log10(observedPvalue)]
    
    renameCorrectionMethods(qqPlotDTSampled, 'Method')
}

getAberrantPlottingData <- function(ols, dataset='Kremer'){
    odls <- ols[[dataset]]
    dt <- rbindlist(lapply(names(odls), function(x){
            dt <- data.table(Naberrant=aberrant(odls[[x]], by='sample'), Method=x)
            dt[, Sample_rank:=rank(Naberrant, ties.method = 'first')]
            dt[, medianNaberrant:=median(Naberrant)]
            dt[, nrow:=nrow(odls[[x]])]
            dt
    }))
    renameCorrectionMethods(dt, 'Method')
}


zoomtheme <- theme(legend.position="none", axis.title.x=element_blank(),
        axis.title.y=element_blank(), title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(color='white', fill="white"),
        plot.background = element_rect(color='white', fill = "white"),
        plot.margin = unit(c(0,0,-6,-6),"mm"))


plotFigureQQ <- function(dt, dataset, withInlet=TRUE, range=c(0.5, 3.5, 17, 37)){
    sdt <- dt[, .(nle=neglog10expectedPvalue, nlo=neglog10observedPvalue, Method=Method, nllow=nllower, nlup=nlupper)]
    ggp <- ggplot(sdt, aes(nle, nlo, col=Method)) + 
        geom_point(size=0.8) + 
        scale_color_brewer(palette='Dark2') + 
        geom_abline(intercept = 0, slope = 1) + 
        labs(title = paste(dataset),  
             x=expression(paste(-log[10], " (expected ", italic(P), "-value)")),
             y=expression(paste(-log[10], " (obs. ", italic(P), "-value)"))) + 
        geom_ribbon(data=sdt[Method=='PCA'], col=alpha('gray', 0.2), fill=alpha('gray', 0.5),
                aes(x=nle, ymin = nllow, ymax = nlup))
    
    if(isTRUE(withInlet)){
        ggZoom <- ggplotGrob(
            ggp + coord_cartesian(xlim=c(2, 3.5), ylim=c(2, 7)) + zoomtheme)
        ggp <- ggp + annotation_custom(grob=ggZoom, xmin=range[1], xmax=range[2],
                                       ymin=range[3], ymax=range[4])
    }
    ggp
}

plotAberrantFigure <- function(dt, dataset, textadjust=c(1.4, 1.4)){
    cutoff <- max(dt$nrow) * 0.005
    ggp <- ggplot(dt, aes(Sample_rank, Naberrant, col=Method)) + 
        geom_line() +
        scale_color_brewer(palette='Dark2') + 
        scale_y_log10() +
        geom_hline(linetype = "dotted", aes(yintercept=medianNaberrant, col=Method)) +
        geom_hline(linetype = "dashed", yintercept = cutoff, col='firebrick') + 
        labs(title=paste0( dataset),
            x='Sample rank',
            y='#Aberrant genes') + 
        annotate('text', x=2, y=cutoff*textadjust[1], label='Aberrant samples', 
                hjust='left', fontface='plain') + 
        annotate('text', x=2, y=max(dt$medianNaberrant)*textadjust[2], label='Median',
                hjust='left', fontface='plain')
    ggp
}


plotVolcanoFigure <- function(ods, id, name, ylim, xlim){
    dt <- data.table(
            zScore = c(zScore(ods[,id])), 
            aberrant=c(aberrant(ods[,id])),
            nlp=c(-log10(pValue(ods[,id]))))
    ggplot(dt, aes(zScore, nlp, col=aberrant)) + 
        geom_point() + 
        labs(title = paste(name), 
            y=expression(paste(-log[10], " (", italic(P), "-value)")),
            x='Z-score') + 
        scale_color_manual(values=c('gray', 'firebrick')) + 
        #grids(linetype='dotted') + 
        theme(legend.position="none")
}

#' 
#' # Plots 
#' 

#'
#' ## Q-Q plots
#'
plotdt <- getQQPlottingData(odsls, dataset='Kremer')
ggQQK <- plotFigureQQ(plotdt, 'Kremer', withInlet=FALSE)
ggQQK
plotdt <- getQQPlottingData(odsls, dataset='GTEx')
ggQQG <- plotFigureQQ(plotdt, 'GTEx', withInlet=FALSE, range=c(0.5, 3.7, 22, 49))
ggQQG

#' 
#' ## Aberrant per sample plot
#' 
plotdt <- getAberrantPlottingData(odsls, 'Kremer')
ggApsK <- plotAberrantFigure(plotdt, 'Kremer', textadjust = c(1.5, 1.5))
ggApsK
plotdt <- getAberrantPlottingData(odsls, 'GTEx')
ggApsG <- plotAberrantFigure(plotdt, 'GTEx')
ggApsG

#' 
#' ## Volcano examples
#' 
odsOI <- odsls$GTEx[['peer']]
abSample <- names(sort(aberrant(odsOI, by='s'))[ncol(odsOI)])
subjID <- gsub('(GTEX-[^-]+)-.*', '\\1', abSample)
ylimits <- c(0, 20)
xlimits <- c(-8, 8)
ggVA <- plotVolcanoFigure(odsOI, abSample, paste(subjID, '(PEER)')) +
    xlim(xlimits) + ylim(ylimits)
ggVA 
ggVN <- plotVolcanoFigure(odsls$GTEx[[3]], abSample, paste(subjID, '(OUTRIDER)')) + 
    xlim(xlimits) + ylim(ylimits)
ggVN


###
### merged plot
###
plotLegend <- get_legend(ggQQK + guides(
        col=guide_legend(override.aes=list(size=3), title='Method\n', order=1)))

ggfinal <- plot_grid(ncol=1, align='h',
    plot_grid(labels=c(LETTERS[1:2], ''), rel_widths=c(4,4,1.5), ncol=3,
        ggQQG + theme(legend.position="none"),
        ggQQK + theme(legend.position="none"),
        as_ggplot(plotLegend)),
    plot_grid(labels=LETTERS[3:4],
        ggApsG + theme(legend.position="none"),
        ggApsK + theme(legend.position="none")),
    plot_grid(labels=LETTERS[5:6],
        ggVN, 
        ggVA)
)
ggfinal


ggfinal <- plot_grid(ncol=2, labels=LETTERS[1:6], axis='l',align='v',
        ggQQG + theme(legend.position= c(.08, .63), 
                legend.direction = "vertical", legend.key.size= unit(0.5, "cm")),
        ggQQK + theme(legend.position="none"),
        ggApsG + theme(legend.position="none"),
        ggApsK + theme(legend.position="none"),
        ggVA,
        ggVN
)
ggfinal



#' 
#' # Save figures as PNG and PDF
#' 
ggsave(pdfOut, ggfinal, 'pdf', width=9, height=9)
ggsave(pngOut, ggfinal, 'png', width=9, height=9, dpi=900)


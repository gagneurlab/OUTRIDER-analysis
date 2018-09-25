#'---
#' title: Figure S4 mu simulation
#' author: Felix Brechtmann, Christian Mertes
#' wb:
#'   input: 
#'     - odsSimNB:       '`sm expand(config["DATADIR"] + "/RecallAnalysis/nsamples=All/inj={inj}/{inj_value}/correction={correction}/q=best/SimulationNBinom_fitN_Q10_OUTRIDERfit.RDS", inj=config["inj"], inj_value=config["inj_value"], correction=["pca", "peer", config["FIGURE_IMPLEMENTATION"]])`'
#'     - odsSimNBFit:    '`sm expand(config["DATADIR"] + "/RecallAnalysis/nsamples=All/inj={inj}/{inj_value}/correction={correction}/q=best/SimulationNBinom_fitY_Q10_OUTRIDERfit.RDS", inj=config["inj"], inj_value=config["inj_value"], correction=["pca", "peer", config["FIGURE_IMPLEMENTATION"]])`'
#'     - odsSimNorm:     '`sm expand(config["DATADIR"] + "/RecallAnalysis/nsamples=All/inj={inj}/{inj_value}/correction={correction}/q=best/SimulationNorm_fitN_Q10_OUTRIDERfit.RDS", inj=config["inj"], inj_value=config["inj_value"], correction=["pca", "peer", config["FIGURE_IMPLEMENTATION"]])`'
#'     - odsSimNormFit:  '`sm expand(config["DATADIR"] + "/RecallAnalysis/nsamples=All/inj={inj}/{inj_value}/correction={correction}/q=best/SimulationNorm_fitY_Q10_OUTRIDERfit.RDS", inj=config["inj"], inj_value=config["inj_value"], correction=["pca", "peer", config["FIGURE_IMPLEMENTATION"]])`'
#'   output:
#'     - pdfOut: '`sm config["FIGDIR"] + "/FigureS05_mu_plot.pdf"`'
#'     - pngOut: '`sm config["FIGDIR"] + "/FigureS05_mu_plot.png"`'
#' output:
#'   html_document:
#'     code_folding: show
#'     code_download: TRUE 
#'---
#'

if(FALSE){
    snakemake <- readRDS('tmp.snakemake.RDS'); snakemake 
    getFiles <- function(dName){
        DATADIR <- '/data/ouga04b/ag_gagneur/project_local/scared/paper/revision/run0910-final/data/'
        BENCHMARK_DIR <- 'RecallAnalysis'
        methodRegex <- paste0('/correction=(pca|peer|', CONFIG_YAML$FIGURE_IMPLEMENTATION, ')/')
        grep(methodRegex, value=TRUE,
                list.files(path=file.path(DATADIR, BENCHMARK_DIR), recursive=TRUE, 
                        full.names=TRUE, 
                        pattern=paste0('Simulation', dName, '_Q10_OUTRIDERfit.RDS')))
    }
    filesNB    <- getFiles('NBinom_fitN')
    filesNBFit <- getFiles('NBinom_fitY')
    filesLN   <- getFiles('Norm_fitN')
    filesLNFit <- getFiles('Norm_fitY')
    pdfOut <- './tmp.pdf'
}

source('./src/r/config.R')

pdfOut <- snakemake@output$pdfOut
pngOut <- snakemake@output$pngOut


computeMeanTable <- function(ods, injZscore, method){
    dt <- data.table('mu_true'=c(assay(ods, 'true_mean')),
                     'mu_hat'=c(assay(ods, 'normalizationFactors')),
                     'outlier'=c(assay(ods, 'inj_mask')),
                     'outlierGene' = rowAnys(assay(ods, 'inj_mask')!=0),             
                     'injZscore'=injZscore,
                     'Method'=method)
    return(dt)
}

plotDT <- function(l_files){
    methods <- unique(gsub('.*=', '', basename(dirname(dirname(l_files)))))
    l_files <- l_files[grep('/inj=both/', l_files)]
    l_files <- l_files[rowAnys(sapply(methods, 
            function(x) grepl(paste0('/correction=', x, '/'), l_files)))]
    
    odslist <- bplapply(l_files, readRDS, MulticoreParam(4, progressbar = TRUE))
    names(odslist) <- gsub('/.*$', '', gsub('^.*correction=', '', l_files, perl = TRUE), perl = TRUE)
    injZscores <- as.factor(gsub('/.*$', '', gsub('^.*zscore=', '', l_files, perl = TRUE), perl = TRUE))
    
    dt <- lapply(1:length(odslist), function(i) computeMeanTable(odslist[[i]], injZscores[i], names(odslist)[i])) %>% rbindlist()
    dt[,log2error:=(log2(mu_true) - log2(mu_hat))^2]
    dt[,mean_bins:=cut(mu_true, breaks=10^seq(-3, 7, length.out = 10))]
    
    dt <- renameCorrectionMethods(dt, 'Method')
    return(dt)
}

getMuPlot <- function(dt, title, colOrder=seq_along(unique(dt$Method))){
    
    dt <- simplifyBinRanges(dt, 'mean_bins')
    col <- brewer.pal(length(unique(dt$Method)), 'Dark2')
    
    gg <- ggplot(dt[injZscore==2], aes(mean_bins, log2error, col=Method)) + 
        geom_boxplot() + 
        scale_y_log10(labels = trans_format("log10", math_format(10^.x))) + 
        scale_color_manual(values=col[colOrder], 
                labels=c(expression(OUTRIDER[fitted]), 
                         expression(OUTRIDER[simulated]),
                         expression(PCA), 
                         expression(PEER))) + 
        theme(axis.text.x=element_text(angle=45, hjust=1)) +
        labs(y=expression((log[2](hat(mu)) - log[2](mu))^2), x='Gene mean',
                title=title)
    gg
}

#' 
#' # Read in NB simulation.
#' 
filesNB <- snakemake@input$odsSimNB
bestQNB <- getBestQ(readRDS(filesNB[[1]]))
dtnb <- plotDT(filesNB)
dtnb[, Method:=gsub('OUTRIDER', paste0('OUTRIDER (q==simualted)'), Method)]

filesNBFit <- snakemake@input$odsSimNBFit
bestQNBFit <- getBestQ(readRDS(filesNBFit[[1]]))
dtnbFit <- plotDT(filesNBFit)
dtnbFit[, Method:=gsub('OUTRIDER', paste0('OUTRIDER (q==fitted)'), Method)]

nbPlotBoth <- rbind(dtnb, dtnbFit[grepl('^OUTRIDER', Method)])
nbPlot <- getMuPlot(nbPlotBoth, paste0('Simulation (NB, q=', bestQNB, ')'), c(1,4,2,3))
nbPlot


#' 
#' # Read in Normal
#' 
filesLN <- snakemake@input$odsSimNorm
bestQLN <- getBestQ(readRDS(filesLN[[1]]))
dtnorm <- plotDT(filesLN)
dtnorm[, Method:=gsub('OUTRIDER', paste0('OUTRIDER (q[simulated])'), Method)]

filesLNFit <- snakemake@input$odsSimNormFit
bestQLNFit <- getBestQ(readRDS(filesLNFit[[1]]))
dtnormFit <- plotDT(filesLNFit)
dtnormFit[, Method:=gsub('OUTRIDER', paste0('OUTRIDER (q[fitted])'), Method)]

normalPlotBoth <- rbind(dtnorm, dtnormFit[grepl('^OUTRIDER', Method)])
normalPlot <- getMuPlot(normalPlotBoth, paste0('Simulation (LN, q=', bestQLN, ')'), c(1,4,2,3))
normalPlot


#' 
#' # Assemble figure
#' 
relHeight <- c(4,1)
ggAll <- plot_grid(rel_heights = relHeight, ncol=1,
    plot_grid(labels=LETTERS[1:2],
        nbPlot + theme(legend.position="none"), 
        normalPlot + theme(legend.position="none")),
    plot_grid(get_legend(nbPlot + theme(legend.direction="horizontal", 
            legend.justification="center", legend.box.just="bottom",
            legend.text.align=0, legend.position="bottom"))))
ggAll


#' 
#' # Facet plots 
#' 
ggplot(nbPlotBoth, aes(mean_bins, log2error, col=Method)) + 
    geom_boxplot() + scale_y_log10() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    facet_grid(injZscore ~ outlier) + 
    labs(title='Simulation\nNegative Binomial')

ggplot(normalPlotBoth, aes(mean_bins, log2error, col=Method)) + 
    geom_boxplot() + scale_y_log10() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    facet_grid(injZscore ~ outlier) + 
    labs(title='Simulation\nLog Normal')

#' 
#' # Save figures
#' 
ggsave(pdfOut, ggAll, 'pdf', width=8, height=8)
ggsave(pngOut, ggAll, 'png', width=8, height=8, dpi=900)



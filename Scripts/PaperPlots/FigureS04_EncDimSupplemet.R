#'---
#' title: Figure S3 zscore injection impact
#' author: Christian Mertes, Felix Brechtmann
#' wb:
#'   py: [ 'EncDimDataSets=["Skin_Not_Sun_Exposed_Suprapubic", "Kremer", "SimulationNBinom_fitN_Q10", "SimulationNorm_fitN_Q10"]' ]
#'   input: 
#'     - outrider: '`sm expand(config["DATADIR"] + "/zscore_with_encDim_{dataset}_imp_" + config["FIGURE_IMPLEMENTATION"] + ".RDS", dataset=EncDimDataSets)`'
#'     - peer:     '`sm expand(config["DATADIR"] + "/fitOutrider/peer/{dataset}_ODS.RDS", dataset=EncDimDataSets)`'
#'     - pca:      '`sm expand(config["DATADIR"] + "/zscore_with_encDim_{dataset}_imp_pca.RDS", dataset=EncDimDataSets)`'
#'   output:
#'     - pdfOut: '`sm config["FIGDIR"] + "/FigureS04_zscoreSearch.pdf"`'
#'     - pngOut: '`sm config["FIGDIR"] + "/FigureS04_zscoreSearch.png"`'
#' output:
#'   html_document:
#'     code_folding: show
#'     code_download: TRUE
#'---
#'

if(FALSE){
    
    snakemake <- readRDS('tmp.snakemake.RDS'); snakemake 
    
    odsFiles <- list(
        OUTRIDER = c(
            GTEx = 'Output/data/zscore_with_encDim_SimulationNBinom_imp_NLas_TCNo.RDS',
            Kremer = 'Output/data/zscore_with_encDim_SimulationNBinom_imp_NLas_TCNo.RDS',
            'Simulation NB (q=5)' = 'Output/data/zscore_with_encDim_SimulationNBinom_imp_NLas_TCNo.RDS'),
        PEER = c(
            GTEx = 'Output/data/fitOutrider/peer/Skin_Not_Sun_Exposed_Suprapubic_ODS.RDS',
            Kremer = 'Output/data/fitOutrider/peer/Kremer_ODS.RDS',
            'Simulation NB (q=5)' = 'Output/data/fitOutrider/peer/SimulationNBinom_ODS.RDS'))
    
    odsSim <- readRDS('/data/ouga04b/ag_gagneur/project_local/scared/paper/revision/run0822_v6p/data/bestQFitODS/SimulationNBinom_ODS.RDS')
    odsKremer <- readRDS('/data/ouga04b/ag_gagneur/project_local/scared/paper/revision/run0822_v6p/data/filteredODS/Kremer_ODS.RDS')
    
    
    odsList <- list(odsSim, odsKremer)
    encDimParams <- list('Sim'=c(2,4,5,6,8,9,10,11,12,13,15,20),
                         'Kremer'=c(seq(5,55,2), 60, 70, 90))
    
    register(MulticoreParam(threads, length(encDimParams[2])*6, recursive=TRUE, progressbar=TRUE))
    #register(SerialParam())
    
    odsList <- lapply(1:2, function(i) findInjectZscore(odsList[[i]], 
            evalAucPRLoss=TRUE, implementation='newED', 
            encDimParams = encDimParams[[i]], zScoreParams=c(1.5,2:5,'lnorm')))
    
    saveRDS(odsList, 'odsListEncDimJustinCase.RDS')
}

source('src/r/config.R')

#' 
#' # Read input parameters
#' 
pdfOut <- snakemake@output$pdfOut
pngOut <- snakemake@output$pngOut
odsFiles <- list(
    OUTRIDER = snakemake@input$outrider,
    PCA = snakemake@input$pca,
    PEER = snakemake@input$peer)

names(odsFiles[['OUTRIDER']]) <- gsub(
        'zscore_with_encDim_|_imp_.*.RDS', '', basename(odsFiles[['OUTRIDER']]))
names(odsFiles[['PCA']]) <- gsub(
    'zscore_with_encDim_|_imp_.*.RDS', '', basename(odsFiles[['PCA']]))
names(odsFiles[['PEER']]) <- gsub('_ODS.RDS', '', basename(odsFiles[['PEER']]))

#'
#' plot Functions
getPlotEncDim <- function(files){
    odsls <- lapply(files, readRDS)
    odslsDims <- lapply(names(files), function(x){
        fn <- file.path(CONFIG_YAML$DATADIR, 'bestQFitODS', paste0(x, '_ODS.RDS'))
        dim(readRDS(fn))
    })
    names(odslsDims) <- names(odsls)
    plotList <- lapply(names(odsls), function(x) {
        plotZscoreInjSearch(odsls[[x]], digits=3) + 
            theme(legend.position="none") +
            ggtitle(renameDatasetNames(x)) +
            scale_x_log10(limits=c(2, odslsDims[[x]][2]*0.45)) + 
            labs(y='AUC\nprecision-recall')
    })
    return(list(plots=plotList, pdsls=odsls))
}

#'
#' # Get data and plots
#' 
outriderPlots <- getPlotEncDim(odsFiles[['OUTRIDER']])
pcaPlots <- getPlotEncDim(odsFiles[['PCA']])
peerPlots <- lapply(names(odsFiles[['PEER']]), function(x){
    plotPEERAlpha(readRDS(odsFiles[['PEER']][[x]])) + 
        ggtitle(renameDatasetNames(x)) 
})

#' 
#' # Single method
plot_grid(plotlist = outriderPlots$plots)
plot_grid(plotlist = pcaPlots$plots)
plot_grid(plotlist = peerPlots)

plot_grid(ncol=2,
    plot_grid(plotlist = outriderPlots$plots),
    plot_grid(plotlist = pcaPlots$plots)
)
#' 
#' # Assamble figure
#' 
outriderPlots$plots[[2]] <- outriderPlots$plots[[2]] + theme(legend.position="right")
#gg <- plot_grid(rel_widths = c(1.3, 1), 
#    plot_grid(plotlist=outriderPlots$plots,  ncol=1, axis='r', align='v', labels=LETTERS[c(1,3,5,7)]),
#    plot_grid(plotlist=peerPlots, ncol=1, axis='r', align='v', labels=LETTERS[c(2,4,6,8)]))
gg <- plot_grid(
    plotlist=outriderPlots$plots,  ncol=1, axis='r', align='v', labels=LETTERS[1:4])
gg

#'
#' # Save figure
#' 
ggsave(pdfOut, gg, 'pdf', width=8, height=9)
ggsave(pngOut, gg, 'png', width=8, height=9, dpi=900)


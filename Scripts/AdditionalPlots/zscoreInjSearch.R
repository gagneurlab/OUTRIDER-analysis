#'---
#' title: Run Grid Search for Zscore and Encoding Dimension.
#' author: Christian Mertes
#' wb:
#'   threads: 40
#'   input: 
#'     - ods: '`sm config["DATADIR"] + "/bestQFitODS/{dataset}_ODS.RDS"`'
#'   output:
#'     - res: '`sm config["DATADIR"] + "/zscore_with_encDim_{dataset}_imp_{impl}.RDS"`'
#'     - wBhtml: 'Output/html/zscore_with_encDim_{dataset}_imp_{impl}.html'
#'   type: noindex
#' html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---

if(FALSE){
    odsFile <- "./Output/data/bestQFitODS/Kremer_ODS.RDS"
    odsFile <- "./Output/data/bestQFitODS/SimulationNBinom_ODS.RDS"    
    outFile <- "./Output/data/wbuild/zscore-with-encDim-GTEx_not_sun_exposed_Outrider.RDS"
    encDimParam <- seq(7, 13, 2)
    zScore <- c(2,3,4,'lnorm')
    impl <- "NLas_TCNo"
    threads <- 20
    
}

source("./src/r/config.R")


#' 
#' # Input params
threads <- snakemake@threads
odsFile <- snakemake@input$ods
outFile <- snakemake@output$res
impl    <- snakemake@wildcards$impl


#' 
#' load ODS
ods <- readRDS(odsFile)
q <- getBestQ(ods)
encDimParam <- unique(pmin(min(dim(ods)-1), pmax(2, sort(c(
    seq(-50, 50, 5) + q, 
    seq(-20, 20, 2) + q, 
    seq(-10, 10, 1) + q,
    60, 70, 90, 130)))))
zScore  <- c(1.5, 2:5, 'lnorm')

odsFile
impl
encDimParam
zScore
threads

#' # Load config
register(MulticoreParam(threads, length(encDimParam)*length(zScore), 
        recursive=TRUE, progressbar=TRUE))

#' # Find the optimal encoding dimension
ods <- findInjectZscore(ods, encDimParams=encDimParam, zScoreParams=zScore, 
        implementation=impl, BPPARAM=bpparam())

#' # Results
plotZscoreInjSearch(ods)
res <- metadata(ods)[['optimalZscoreEncDim']]

#' 
#' # Result table
#' 
DT::datatable(res)

#' save results
saveRDS(res, outFile)


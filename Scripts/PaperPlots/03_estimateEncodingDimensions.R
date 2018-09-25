#'---
#' title: Find best encoding dimension
#' author: Felix Brechtmann
#' wb:
#'  threads: 20
#'  params: 
#'    implementation: '`sm config["AE_IMPLEMENTATION"]`'
#'  input: 
#'    - ods: '`sm config["DATADIR"] + "/filteredODS/{dataset}_ODS.RDS"`'
#'  output:
#'   - ods:     '`sm config["DATADIR"] + "/bestQFitODS/{dataset}_ODS.RDS"`'
#'   - ggplots: '`sm config["DATADIR"] + "/bestQFitODS/{dataset}_ggplots.RDS"`'
#'   - wBhtml:  '`sm config["htmlOutputPath"] + "/bestQFitODS/{dataset}.html"`'
#'  type: noindex
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---

if(FALSE){
    snakemake <- readRDS("tmp.snakemake.RDS"); snakemake
}

#' # Load configuration
#' * Read input
odsFile <- snakemake@input$ods
rdsOut  <- snakemake@output$ods
ggplotsOut <- snakemake@output$ggplots
dataset <- snakemake@wildcards$dataset
threads <- snakemake@threads
implementation <- snakemake@params$implementation

#' source config
source("./src/r/config.R")
register(MulticoreParam(threads, 150))

#' 
#' # Read data
#' 
dataset
odsFile
ods <- readRDS(odsFile)
dim(ods)

#' 
#' # Set encoding search 
#' 
encDimParam <- c(
        -10:10 + estimateBestQ(ods),
        -10:10*2 + estimateBestQ(ods),
        seq(5, ceiling(ncol(ods)*0.33), 5),
        round(ncol(ods) * c(0.5, 0.75))
)

#' Make sure it stays in the bounds
topBound <- ceiling(min(dim(ods) - 1)*0.75)
bottomBound <- 2
encDimParam <- rev(unique(sort(pmin(topBound, pmax(bottomBound, encDimParam)))))

encDimParam

#' 
#' # Find best encoding dimension
#' 
ods <- findEncodingDim(ods, lnorm=TRUE, implementation=implementation, 
        params=encDimParam, evalAucPRLoss=TRUE, iterations=10)

#' # Result of encoding search
#' * Best encoding dimension
getBestQ(ods)
ggEncDim <- plotEncDimSearch(ods)
ggEncDim

DT::datatable(metadata(ods)[['encDimTable']])

#' 
#' # Save results
#' 
rdsOut
saveRDS(ods, rdsOut)
saveRDS(list(ggEncDim=ggEncDim), ggplotsOut)

#'---
#' title: Run OUTRIDER
#' author: Felix Brechtmann, Christian Mertes
#' wb:
#'  threads: 10
#'  input: 
#'   - ods:     '`sm config["DATADIR"] + "/bestQFitODS/{dataset}_ODS.RDS"`'
#'  output:
#'   - ods:    '`sm config["DATADIR"] + "/fitOutrider/{method}/{dataset}_ODS.RDS"`'
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/fitOutrider/{method}/{dataset}.html"`'
#'  type: noindex
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---
#'

if(FALSE){
    snakemake <- readRDS('tmp.snakemake.RDS'); snakemake
    method  <- 'pca'
    dataset <- 'Kremer'
    odsFile <- paste0('Output/data/bestQFitODS/', dataset, '_ODS.RDS')
    odsOut <- paste0('Output/data/fitOutrider/', method, '/', dataset, '_ODS.RDS')
    threads <- 10
}

#' # Load configuration
#' * Read input
odsFile <- snakemake@input$ods
odsOut  <- snakemake@output$ods
dataset <- snakemake@wildcards$dataset
method  <- snakemake@wildcards$method
threads <- snakemake@threads * 2

#' source config
source("./src/r/config.R")

#' * Read data
dataset
odsFile
ods <- readRDS(odsFile)

#' Read in best q
bestQ <- getBestQ(ods)
bestQ
plotEncDimSearch(ods)

#' # Run full OUTRIDER pipeline
register(MulticoreParam(threads))

message <- paste(date(), ":\nRun for: ", odsOut, 
        '\nwith impl:', method, '\nq: ', bestQ)
message(message)
message
start <- Sys.time()

ods <- OUTRIDER(ods, q=bestQ, autoCorrect=TRUE, implementation=method)

Sys.time() - start

#' * Save results
saveRDS(ods, odsOut)

#'
#' # Plot results
plotAberrantPerSample(ods, main=paste('Aberrant per samples:', method, dataset))
plotQQ(ods, global=TRUE, main=paste('Global QQ for:', method, dataset))

if(grepl('PEER', method, ignore.case=TRUE)){
    plotPEERAlpha(ods)
}


#'---
#' title: Recall All Method
#' author: Christian Mertes
#' wb:
#'   threads: 40
#'   input:
#'     - evalTable: '`sm expand(OUTDIR + "/nsamples={nsamples}/inj={inj}/{inj_value}/correction={correction}/q={q}/pmethod={pmethod}/{{dataset}}_plot.tsv", nsamples=config["N_samples"], inj=config["inj"], inj_value=config["inj_value"], q=config["Qs"], correction=config["correction"], pmethod=config["FDR_METHOD"])`'
#'   output: 
#'     - ggplots: '`sm config["htmlOutputPath"] + "/RecallBenchmark/da_{dataset}_all_ggplots.RDS"`'
#'     - wBhtml: '`sm config["htmlOutputPath"] + "/RecallBenchmark/da_{dataset}_all.html"`'
#'   type: noindex
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---


if(FALSE){
    snakemake <- readRDS('tmp.snakemake.RDS'); snakemake
}

#' load config 
source('./src/r/config.R')

register(MulticoreParam(20))

recallFiles <- snakemake@input$evalTable
ggplotFile <- snakemake@output$ggplots
dataset <- snakemake@wildcards$dataset
correction <- snakemake@wildcards$correction
maxRows <- 1e5

recallFiles
ggplotFile
dataset
correction


#' Extract data in parallel
prDatals <- mclapply(recallFiles, mc.allow.recursive = TRUE, mc.preschedule = FALSE,
        mc.cores=min(12, length(recallFiles)), FUN=readBootstrapData)
prData <- rbindlist(prDatals)
title <- paste("Precision-Recall:", unique(dataset), '(All)', collapse=" ")

ggAll <- plotRibbonBenchmark(prData, title)
ggAll


#' Extract data in parallel
prDatals <- mclapply(recallFiles, kzcut=3, mc.allow.recursive = TRUE, mc.preschedule = FALSE,
        mc.cores=min(12, length(recallFiles)), FUN=readBootstrapData)
prData <- rbindlist(prDatals)
title <- paste("Precision-Recall:", unique(dataset), '(kz<3)', collapse=" ")

ggKZ <- plotRibbonBenchmark(prData, title)
ggKZ


#' 
#' # Pvalue ranked
#' 
prDatals <- mclapply(recallFiles, mc.allow.recursive = TRUE, 
                     mc.cores=min(6, length(recallFiles)), FUN=readBootstrapData, pvalue=TRUE)
prData <- rbindlist(prDatals)
title <- paste("Precision-Recall:", unique(dataset), '(All & P-value)', collapse=" ")

ggPval <- plotRibbonBenchmark(prData, title)
ggPval

#' 
#' # Pvalue ranked cutoff
#' 
prDatals <- mclapply(recallFiles, kzcut=3, mc.allow.recursive = TRUE, 
                     mc.cores=min(6, length(recallFiles)), FUN=readBootstrapData, pvalue=TRUE)
prData <- rbindlist(prDatals)
title <- paste("Precision-Recall:", unique(dataset), '(kz<3 & P-value)', collapse=" ")

ggPvalKZ <- plotRibbonBenchmark(prData, title)
ggPvalKZ


#' 
#' # Pvalue/Zscore mix ranked
#' 
prDatals <- mclapply(recallFiles, mc.allow.recursive = TRUE, byFile=TRUE,
        mc.cores=min(6, length(recallFiles)), FUN=readBootstrapData, pvalue=TRUE)
prData <- rbindlist(prDatals)
title <- paste("Precision-Recall:", unique(dataset), '(All & Mix)', collapse=" ")

ggMix <- plotRibbonBenchmark(prData, title)
ggMix

#' 
#' # Pvalue/Zscore mix ranked cutoff
#' 
prDatals <- mclapply(recallFiles, kzcut=3, mc.allow.recursive = TRUE, byFile=TRUE,
        mc.cores=min(6, length(recallFiles)), FUN=readBootstrapData, pvalue=TRUE)
prData <- rbindlist(prDatals)
title <- paste("Precision-Recall:", unique(dataset), '(kz<3 & Mix)', collapse=" ")

ggMixKZ <- plotRibbonBenchmark(prData, title)
ggMixKZ



#' 
#' # Save data
#' 
saveRDS(list(ggAll, ggKZ, ggPval, ggPvalKZ, ggMix, ggMixKZ), file=ggplotFile)


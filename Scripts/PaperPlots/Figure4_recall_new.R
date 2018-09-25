#'---
#' title: Figure 4 new recall / bootstraped
#' author: Christian Mertes
#' wb:
#'   threads: 40
#'   input:
#'     - evalTable: '`sm expand(OUTDIR + "/nsamples=All/inj={inj}/{inj_value}/correction={correction}/q=best/pmethod={pmethod}/Skin_Not_Sun_Exposed_Suprapubic_plot.tsv", inj=config["inj"], inj_value=config["inj_value"], correction=PAPER_METHODS, pmethod=config["FDR_METHOD"])`'
#'   output: 
#'    - figurePdf: '`sm config["FIGDIR"] + "/Figure4_outrider_recall.pdf"`'
#'    - figurePng: '`sm config["FIGDIR"] + "/Figure4_outrider_recall.png"`'
#'    - ggplots:   '`sm config["FIGDIR"] + "/Figure4_outrider_recall_ggplots.RDS"`'
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---

if(FALSE){
    snakemake <- readRDS('tmp.snakemake.RDS'); snakemake
    max.mc <- 18
}

#' load config 
source('./src/r/config.R')

register(MulticoreParam(20))
max.mc <- 12 

recallFiles <- snakemake@input$evalTable
ggplotFile <- snakemake@output$ggplots
outPdf <- snakemake@output$figurePdf
outPng <- snakemake@output$figurePng

maxRows <- 1e5
dataset <- 'GTEx'

recallFiles
ggplotFile


#' Extract data in parallel
prDatals <- mclapply(recallFiles, mc.allow.recursive=TRUE, mc.preschedule=FALSE, methodOnly=TRUE,
        mc.cores=min(max.mc, length(recallFiles)), FUN=readBootstrapData, pvalue=TRUE, zscoreForAll=TRUE)
zsDatals <- mclapply(recallFiles, mc.cores=min(max.mc, length(recallFiles)), readZscoreCutoff)
prData <- rbindlist(prDatals)
zsData <- rbindlist(zsDatals)

zsData <- renameCorrectionMethods(zsData, 'correction')
zsData <- correctPrecisionRankPlotNames(zsData)
zsData <- zsData[correction %in% c('PCA', 'PEER')]
zsData
prData <- renameCorrectionMethods(prData, 'correction')
prData <- correctPrecisionRankPlotNames(prData)
prData

ggAll <- plotRibbonBenchmark(prData, zscoreData=zsData, linetype=c(2,1)) + 
    scale_color_brewer(palette='Dark2') + 
    scale_fill_brewer(palette='Dark2')
ggAll


#' 
#' # Save figures as PNG and PDF
#' 
ggsave(outPdf, ggAll, 'pdf', width=9, height=8)
ggsave(outPng, ggAll, 'png', width=9, height=8, dpi=900)
saveRDS(ggAll, file=ggplotFile)

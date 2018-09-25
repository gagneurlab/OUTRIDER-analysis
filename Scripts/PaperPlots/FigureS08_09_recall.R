#'---
#' title: Figure 4 new recall / bootstraped
#' author: Christian Mertes
#' wb:
#'   py: [ 'SUPL_CORRECTIONS=PAPER_METHODS + ["baseCooks", "basePearsonRes"]' ]
#'   threads: 40
#'   params: 
#'     - methods: '`sm SUPL_CORRECTIONS`'
#'   input:
#'     - evalTableGTEx:   '`sm expand(OUTDIR + "/nsamples=All/inj={inj}/{inj_value}/correction={correction}/q=best/pmethod={pmethod}/Skin_Not_Sun_Exposed_Suprapubic_plot.tsv", inj=config["inj"], inj_value=config["inj_value"], correction=SUPL_CORRECTIONS, pmethod=config["FDR_METHOD"])`'
#'     - evalTableKremer: '`sm expand(OUTDIR + "/nsamples=All/inj={inj}/{inj_value}/correction={correction}/q=best/pmethod={pmethod}/Kremer_plot.tsv", inj=config["inj"], inj_value=config["inj_value"], correction=SUPL_CORRECTIONS, pmethod=config["FDR_METHOD"])`'
#'   output: 
#'    - pdfGtex:   '`sm config["FIGDIR"] + "/FigureS08_recall_gtex.pdf"`'
#'    - pngGtex:   '`sm config["FIGDIR"] + "/FigureS08_recall_gtex.png"`'
#'    - pdfKremer: '`sm config["FIGDIR"] + "/FigureS09_recall_Kremer.pdf"`'
#'    - pngKremer: '`sm config["FIGDIR"] + "/FigureS09_recall_Kremer.png"`'
#'    - ggplots:   '`sm config["FIGDIR"] + "/FigureS08_09_recall_ggplots.RDS"`'
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---

if(FALSE){
    snakemake <- readRDS('tmp.snakemake.RDS'); snakemake
    mc.cores <- 20
}

#' load config 
source('./src/r/config.R')

ggplotFile <- snakemake@output$ggplots
maxRows <- 1e5
mc.cores <- 12

getBenchPlot <- function(recallFiles){
    
    prDatals <- mclapply(recallFiles, mc.allow.recursive=TRUE, 
            mc.preschedule=FALSE, mc.cores=mc.cores,
            FUN=readBootstrapData, pvalue=TRUE, maxRows=maxRows, zscoreForAll=TRUE)
    zsDatals <- mclapply(recallFiles, mc.cores=min(mc.cores, length(recallFiles)), readZscoreCutoff)
    prData <- rbindlist(prDatals)
    zsData <- rbindlist(zsDatals)
    
    zsData <- renameCorrectionMethods(zsData, 'correction')
    zsData <- correctPrecisionRankPlotNames(zsData)
    zsData <- zsData[correction %in% c('OUTRIDER', 'PCA', 'PEER')]
    
    prData <- renameCorrectionMethods(prData, 'correction')
    prData <- correctPrecisionRankPlotNames(prData)
    
    col <- brewer.pal(5, 'Dark2')[c(4,1,2,5,3)]
    breaks <- sort(unique(prData$correction))[c(2,3,5,1,4)]
    
    lt <- c(4,2,3,1)
    ltBreaks <- sort(unique(prData$type))[lt]
    
    gg <- plotRibbonBenchmark(prData, zscoreData=zsData, linetype=lt, maxRows=maxRows) + 
        scale_color_manual(values=col, breaks=breaks) + 
        scale_fill_manual(values=col, breaks=breaks) + 
        scale_linetype_manual(values=lt, breaks=ltBreaks)
    
    return(list(plot=gg, dt=prData))
}

#'
#' # GTEx benchmark
#' 
recallFiles <- snakemake@input$evalTableGTEx
recallFiles

#' Plot it
ggGTEx <- getBenchPlot(recallFiles)
ggGTEx$plot


#'
#' # Kremer benchmark
#' 
recallFiles <- snakemake@input$evalTableKremer
recallFiles

#' Plot it
ggKremer <- getBenchPlot(recallFiles)
ggKremer$plot


#' 
#' # Save figures as PNG and PDF
#' 
width  <- 9
height <- 8
ggsave(snakemake@output$pdfGtex, ggGTEx$plot,     'pdf', width=width, height=height)
ggsave(snakemake@output$pngGtex, ggGTEx$plot,     'png', width=width, height=height, dpi=900)
ggsave(snakemake@output$pdfKremer, ggKremer$plot, 'pdf', width=width, height=height)
ggsave(snakemake@output$pngKremer, ggKremer$plot, 'png', width=width, height=height, dpi=900)

saveRDS(list(ggGTEx$plot, ggKremer$plot), file=ggplotFile)

#'---
#' title: Figure 5 Recall
#' author: Felix Brechtmann
#' wb:
#'   threads: 40
#'   input: 
#'     - benchRes:  '`sm expand(config["DATADIR"] + "/RecallAnalysis/nsamples=All/inj=both/zscore=4/correction={method}/q=best/pmethod=BY/Skin_Not_Sun_Exposed_Suprapubic_OUTRIDERpVal.RDS", method=PAPER_METHODS)`'
#'   output:
#'     - figurePdf: '`sm config["FIGDIR"] + "/FigureS10_recall_single_bin.pdf"`'
#'     - figurePng: '`sm config["FIGDIR"] + "/FigureS10_recall_single_bin.png"`'
#' output:
#'   html_document:
#'     code_folding: show
#'     code_download: TRUE
#'---
#'

if(FALSE){
    snakemake <- readRDS('tmp.snakemake.RDS'); snakemake    
    odsFiles <- "Output/data/RecallAnalysis/nsamples=All/inj=both/zscore=6/correction=NLas_TCNo/q=best/pmethod=BY/Skin_Not_Sun_Exposed_Suprapubic_OUTRIDERpVal.RDS"
    FIGURE_NAME <- "FigureS4_outrider_recall_bin"
}

source('src/r/config.R')

maxRows <- 1e5
pdfOut <- snakemake@output$figurePdf
pngOut <- snakemake@output$figurePng

#'
#' Extraction and plotting functions
#' 
getPlotData <- function(file, bins=16, mc.cores=bins){
    message(date(), ': work on file: ', file)
    
    ods <- readRDS(file)
    q <- getBestQ(ods)
    params <- gsub('.*=', '', rev(unlist(strsplit(file[[1]], '/'))))[1:7]
    
    evalTable <- createEvalTable(ods=ods, q=q, correction=params[4], 
            Nsamples=params[7],
            inj=params[6], inj_value=paste0('zscore=', params[5]),
            pmethod=params[2], dataset=gsub('_OUTRIDER.*', '', params[1]), 
            nMeanBins=bins)
    
    bootstrapCall <- function(x) readBootstrapData(x, pvalue=TRUE, 
            zscoreForAll=TRUE, byFile=FALSE, maxRows=maxRows)
    
    dt1 <- bootstrapCall(evalTable)
    
    binlabels <- unique(evalTable$geneMeanBin)
    dt2ls <- mclapply(binlabels, mc.allow.recursive = TRUE, mc.cores=mc.cores, function(x){
        ansdt <- evalTable[geneMeanBin == x]
        ansdt <- ansdt[,numInjected:=sum(trueOutlier!=0)] 
        ans <- bootstrapCall(ansdt)
        ans[,geneMeanBin:=x]
        ans
    })
    
    return(list(
        plotDT1 = dt1,
        plotDT2 = rbindlist(dt2ls)
    ))
}

plotBinBenchmark <- function(datals, listName, wrapFun){
    data <- rbindlist(lapply(datals, '[[', listName))
    
    data <- renameCorrectionMethods(data, 'correction')
    data <- correctPrecisionRankPlotNames(data)
    if('geneMeanBin' %in% colnames(data)){
        data <- simplifyBinRanges(data, 'geneMeanBin')
    }
    
    gg <- plotRibbonBenchmark(data, linetype=c(2,1), wrap_function=wrapFun) +
        scale_color_brewer(palette='Dark2') + 
        scale_fill_brewer(palette='Dark2')
    gg
}


#' 
#' Script for plotting
#' 

#' get input
odsFiles <- snakemake@input$benchRes
names(odsFiles) <- gsub('.*=', '', basename(dirname(dirname(dirname(odsFiles)))))
odsFiles <- odsFiles[CONFIG_YAML$FIGURE_IMPLEMENTATION]
odsFiles

#' get data
allData <- mclapply(odsFiles, getPlotData, mc.cores=length(odsFiles),
        mc.allow.recursive=TRUE, bins=9)

#' plot panels
ggFullBin <- plotBinBenchmark(allData, 'plotDT1', function() facet_grid(inj_value ~ inj))
ggFullBin

ggBinWraps <- plotBinBenchmark(datals=allData, listName='plotDT2', 
        wrapFun=function() facet_wrap('geneMeanBin'))
ggBinWraps

#'
#' # Assembly figure
#' 
ggAll <- plot_grid(ggFullBin, ggBinWraps + theme(legend.position="none"), 
        labels=LETTERS[1:2], ncol=1, rel_heights=c(1.2,2))
ggAll

#'
#' # Save plots
#' 
ggsave(pdfOut, ggAll, width=10, height=11, device='pdf') 
ggsave(pngOut, ggAll, width=10, height=11, device='png', dpi=900)


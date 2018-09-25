#'---
#' title: Sample size analysis for Kremer (data generation)
#' author: Christian Mertes
#' 
#' wb:
#'   threads: 70
#'   params: 
#'     - methods: '`sm PAPER_METHODS`'
#'   input: 
#'     - odsk: '`sm config["DATADIR"] + "/bestQFitODS/Kremer_ODS.RDS"`'
#'   output:
#'     - rds: '`sm config["DATADIR"] + "/singleParameterAnaylsis/SampleSizeKremer.RDS"`'
#' output:
#'   html_document:
#'     code_folding: show
#'     code_download: TRUE
#'---

if(FALSE){
    snakemake <- readRDS('tmp.snakemake.RDS'); snakemake
    odsFile <- 'Output/data/bestQFitODS/Kremer_ODS.RDS'
    methods <- c('pca', 'peer', CONFIG_YAML$FIGURE_IMPLEMENTATION)
}


#' 
#' # Read input and config
#' 
source('./src/r/config.R')
register(MulticoreParam(50, 200, progressbar=TRUE))

odsFile <- snakemake@input$odsk
rdsOut <- snakemake@output$rds
methods <- snakemake@params$methods

odsFile
rdsOut

#' 
#' # Read data and set default values
#' 
ods <- readRDS(odsFile)
ods

set.seed(42)
numRuns <- 5
sampleSizes <- ceiling(ncol(ods) * c(1, 3/4, 2/3, 1/2, 1/3, 1/4, 1/10))

numRuns
sampleSizes

#' # Get samples of interest
benchmarkSet <- data.table(do.call(rbind, list(
    c(sampleID='35791', RNA_ID='MUC1344', geneID='TIMMDC1'),
    c('66744', 'MUC1365', 'TIMMDC1'),
    c('80256', 'MUC1404', 'ALDH18A1'),
    c('58955', 'MUC1350', 'CLPP'),
    c('73804', 'MUC1396', 'MGST1'),
    c('62346', 'MUC1361', 'MCOLN1')
)))

#' 
#' # Run the benchmark
#' 
runSampleSizeBenchmark <- function(i, nSamples, ods, method, BPPARAM=bpparam()){
    
    spick <- c(which(colnames(ods) %in% benchmarkSet$RNA_ID), sample(ncol(ods)))
    spick <- unique(spick)[seq_len(nSamples)]
    
    cods <- OutriderDataSet(countData=counts(ods[,spick]))
    cq <- ceiling(getBestQ(ods) * nSamples/ncol(ods))
    
    cods <- OUTRIDER(cods, cq, autoCorrect=TRUE, method, BPPARAM=BPPARAM)
    res <- results(cods, BPPARAM=BPPARAM)
    
    return(
        list(
            i = i,
            nSamples = nSamples,
            method = method,
            spick = spick,
            q = cq,
            ods = cods,
            res = res
    )) 
}

params <- expand.grid(run=seq_len(numRuns), size=sampleSizes, method=methods)
res <- bplapply(seq_len(nrow(params)), function(j){
    print(paste("running:", paste(params[j, ], collapse=', ')))
    ans <- runSampleSizeBenchmark(i=params[j,1], nSamples=params[j,2],
            ods=ods, method=params[j,3], BPPARAM=SerialParam())
    print(paste(date(), ': done'))
    ans
})

plotDtls <- lapply(res, function(x){
    dt <- rbindlist(apply(benchmarkSet, 1, function(i){
        data.table(
            sampleID = i['RNA_ID'],
            geneID = i['geneID'],
            pvalue = as.vector(pValue(x[['ods']][i['geneID'], i['RNA_ID']])),
            aberrant = aberrant(x[['ods']])[i['geneID'], i['RNA_ID']])
    }))
    dt[,Method:=x$method]
    dt[,q:=x$q]
    dt[,nSamples:=x$nSamples]
    dt[,i:=x$i]
})
plotDt <- rbindlist(plotDtls)
plotDt

#' 
#' # Plot it
#' 
plotDtRe <- renameCorrectionMethods(plotDt, 'Method')
gg <- ggplot(plotDtRe, aes(as.factor(nSamples), -log10(pvalue), col=Method)) + 
    geom_boxplot(position=position_dodge(0.75), outlier.shape = "") + 
    geom_point(aes(col=Method, shape=geneID, alpha=aberrant), position=position_dodge(0.75)) + 
    scale_color_brewer(palette='Dark2') + 
    labs(x='Sample size', y=expression(-log[10](italic(P)-value)),
            shape='HGNC symbol', alpha='Significant') + 
    theme(legend.text.align=0) + 
    grids(linetype='dotted') + 
    scale_alpha_discrete(range=c(0.3,1)) + 
    scale_shape_discrete(
        breaks=unique(plotDtRe$geneID),
        labels=sapply(unique(plotDtRe$geneID), function(xy) {
                parse(text=paste('italic(', xy, ')'))
        }))
gg


#'
#' # Save rds file
#' 
saveRDS(list(results = res, plotDt = plotDt), rdsOut)

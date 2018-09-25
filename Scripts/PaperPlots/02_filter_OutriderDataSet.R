#'---
#' title: Filter ODS file
#' author: Christian Mertes
#' wb:
#'  input: 
#'    - ods: '`sm config["DATADIR"] + "/rawODS/{dataset}_ODS.RDS"`'
#'  output:
#'    - ods:     '`sm config["DATADIR"] + "/filteredODS/{dataset}_ODS.RDS"`'
#'    - ggplots: '`sm config["DATADIR"] + "/filteredODS/{dataset}_ggplots.RDS"`'
#'    - wBhtml:  'Output/html/filteredODS/{dataset}.html'
#'  type: noindex
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---

#
# DEBUG values
# 
if(FALSE){
    snakemake <- readRDS("tmp.snakemake.RDS"); snakemake
    dataset <- "Skin_Not_Sun_Exposed_Suprapubic"
    odsFile <- "Output/data/rawODS/Skin_Not_Sun_Exposed_Suprapubic_ODS.RDS"
}

source('src/r/config.R')

#'
#' Read input
#' 
dataset    <- snakemake@wildcards$dataset
odsFile    <- snakemake@input$ods
ggplotsOut <- snakemake@output$ggplots
odsOut     <- snakemake@output$ods

odsFile
ods <- readRDS(odsFile)

#' 
#' # Filter samples
#' 

#' Remove low coverage samples from data
ggSeqCov <- plotSequencingCoverage(ods, 'GTEx')
ggSeqCov

cov <- colSums(counts(ods))
zsf <- cov > mean(cov)-3*sd(cov)
table(zsf)
ods <- ods[, zsf]

#' 
#' # Filter genes
#' 
ods <- filterExpression(ods, filterGenes=FALSE, savefpkm=TRUE)
ggfpkm <- plotFPKM(ods)
ggfpkm
ods <- ods[mcols(ods)$passedFilter == TRUE]
dim(ods)

#' Remove low coverage genes from data
plot(sort(rowSums(counts(ods) == 0)))
table(rowSums(counts(ods) == 0) < ncol(ods)*0.75)
ods <- ods[rowSums(counts(ods) == 0) < ncol(ods)*0.75]
dim(ods)


#' 
#' # Plot final count correlation
#' 
#+ plot correlation, fig.height=12, fig.width=12
if('AGE' %in% colnames(colData(ods))){
    plotCountCorHeatmap(ods, normalized=FALSE, 
            rowCoFactor="SEX", colCoFactor="AGE", 
            main=paste('Count correlation\nGTEx ', dataset))
} else {
    plotCountCorHeatmap(ods, normalized=FALSE,
            main=paste('Count correlation\nGTEx ', dataset))
}

#'
#' # Save results
#'
saveRDS(list(ggfpkm=ggfpkm, ggSeqCov=ggSeqCov), ggplotsOut)
saveRDS(ods, odsOut)



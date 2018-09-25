#'---
#' title: Supplement Figure 6 Kremer Benchmark
#' author: Christian Mertes
#' wb:
#'   input:
#'     - odsk: '`sm expand(config["DATADIR"] + "/fitOutrider/{method}/Kremer_ODS.RDS", method=PAPER_METHODS)`'
#'   output:
#'     - outPdf: '`sm config["FIGDIR"] + "/FigureS13_Kremer_benchmark.pdf"`'
#'     - outPng: '`sm config["FIGDIR"] + "/FigureS13_Kremer_benchmark.png"`'
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---
#'

if(FALSE){
    snakemake <- readRDS('tmp.snakemake.RDS'); snakemake 
    odsFiles <- c(
        'Output/data/fitOutrider/pca/Kremer_ODS.RDS',
        'Output/data/fitOutrider/peer/Kremer_ODS.RDS',
        'Output/data/fitOutrider/NLas_TCNo/Kremer_ODS.RDS')
        
}


PENAL_CEX <- 1.3
FIGURE_NAME <- basename(gsub('.pdf$', '', snakemake@output$outPdf))

source('src/r/config.R')
library(VennDiagram)
convertCMD <- 'convert -density 300 -trim -quality 100 -flatten'

KREMER_S7_FILE = 'https://media.nature.com/original/nature-assets/ncomms/2017/170612/ncomms15824/extref/ncomms15824-s7.txt'
KREMER_S3_FILE = 'https://media.nature.com/original/nature-assets/ncomms/2017/170612/ncomms15824/extref/ncomms15824-s3.txt'
KREMER_S8_FILE = 'https://media.nature.com/original/nature-assets/ncomms/2017/170612/ncomms15824/extref/ncomms15824-s8.txt'


#' # Get result data
#' 
#' ## Get OUTRIDER results
odsFiles <- snakemake@input$odsk
names(odsFiles) <- c('PCA', 'PEER', 'OUTRIDER')
odsls <- lapply(odsFiles, readRDS)
rawresls <- lapply(odsls, results, padj=FDR_LIMIT, zScore=Z_LIMIT)

DT::datatable(rawresls[['OUTRIDER']])

#'
#' ## Get Kremer et al results
mapping <- data.table(read.table(KREMER_S8_FILE, sep="\t", header=TRUE))
sampleData <- data.table(read.table(KREMER_S3_FILE, sep='\t', header=TRUE))
kremerRes <- data.table(read.table(KREMER_S7_FILE, sep='\t', header=TRUE))

# fix missannotation
sampleData[DIAGNOSIS_GENE == 'MGST1', DIAGNOSIS_GENE:=NA]

#' 
#' * Add fibroblast id
#' 
resls <- lapply(names(rawresls), function(n){
    x <- rawresls[[n]]
    x <- x[,.(RNA_ID=sampleID, geneID=geneID, zScore=zScore, n=TRUE)]
    x <- merge(mapping[,.(RNA_ID, FIBROBLAST_ID)], x, by='RNA_ID', all.y=TRUE)
    x <- x[,.(sampleID=as.character(FIBROBLAST_ID), geneID, RNA_ID, zScore, n)]
    x <- x[,.(zScore=mean(zScore)), by=c("sampleID", "geneID", "RNA_ID", "n")]
    setnames(x, 'n', n)
    setnames(x, 'zScore', paste0("zScore_", n))
    x
})
names(resls) <- names(odsFiles)
odsResMerge <- resls[['OUTRIDER']]

#'
#' ## Merge data
kremerResMerge <- unique(kremerRes[RNA_ABER_EXP_SIGNI == TRUE, .(
        sampleID=as.character(FIBROBLAST_ID), geneID=HGNCID, Kremer=TRUE)])
disease <- data.table(
    do.call(rbind, list(
        c(sampleID='35791', geneID='TIMMDC1'),
        c('66744', 'TIMMDC1'),
        c('80256', 'ALDH18A1'),
        c('58955', 'CLPP'),
        c('73804', 'MGST1'),
        c('62346', 'MCOLN1')
    )))[,.(sampleID, geneID, DISEASE_CAUSING=TRUE)]

resAll <- merge(merge(merge(merge(kremerResMerge, odsResMerge, all=TRUE, 
        by=c('sampleID','geneID')), disease, all=TRUE), 
        resls[['PCA']], by=c('sampleID','geneID', 'RNA_ID'), all=TRUE), 
        resls[['PEER']], by=c('sampleID','geneID', 'RNA_ID'), all=TRUE)

resAll <- merge(data.table(as.data.table(rowData(odsls[['OUTRIDER']])), 
        geneID=rownames(odsls[['OUTRIDER']]), 
        meanCorrected=rowMeans(counts(odsls[['OUTRIDER']], normalized=TRUE))), 
        resAll, by='geneID', all.y=TRUE)
# remove NHDF samples
resAll <- resAll[sampleID != 'NHDF']

#' add unsolved case information
unsolvedCases <- sampleData[is.na(DIAGNOSIS_GENE), FIBROBLAST_ID]
resAll[,unsolvedCase:=sampleID %in% unsolvedCases]

resAll
table(unsolvedCases %in% resAll$sampleID)

#' 
# create plotting data
data2plot <- list(
        Kremer   = resAll[unsolvedCase==TRUE & Kremer==TRUE,          paste(sampleID, geneID)],
        OUTRIDER = resAll[unsolvedCase==TRUE & OUTRIDER==TRUE,        paste(sampleID, geneID)], 
        KNOWN    = resAll[unsolvedCase==TRUE & DISEASE_CAUSING==TRUE, paste(sampleID, geneID)],
        PCA      = resAll[unsolvedCase==TRUE & PCA==TRUE,             paste(sampleID, geneID)],
        PEER     = resAll[unsolvedCase==TRUE & PEER==TRUE,            paste(sampleID, geneID)])
filterdData2plot <- list(
    Kremer   = resAll[unsolvedCase==TRUE & Kremer==TRUE, paste(sampleID, geneID)],
    OUTRIDER = resAll[unsolvedCase==TRUE & OUTRIDER==TRUE, paste(sampleID, geneID)], 
    KNOWN    = resAll[unsolvedCase==TRUE & DISEASE_CAUSING==TRUE, paste(sampleID, geneID)],
    PCA      = resAll[unsolvedCase==TRUE & PCA==TRUE, paste(sampleID, geneID)],
    PEER     = resAll[unsolvedCase==TRUE & PEER==TRUE, paste(sampleID, geneID)])

par(mfrow=c(1,2))
colECDF <- c('Only Kremer'='darkgreen', 'Only OUTRIDER'='firebrick', 'Both'='darkblue')
plot(ecdf(resAll[!is.na(theta) & is.na(OUTRIDER) & !is.na(Kremer), theta]), log='x', col=colECDF[1], do.points=FALSE, 
        main='ECDF of gene dispersion', xlab='Dispersion', xlim=range(resAll[,theta], na.rm=TRUE))
plot(ecdf(resAll[!is.na(theta) & !is.na(OUTRIDER) &  is.na(Kremer), theta]), col=colECDF[2], do.points=FALSE, add=TRUE)
plot(ecdf(resAll[!is.na(theta) & !is.na(OUTRIDER) & !is.na(Kremer), theta]), col=colECDF[3],  do.points=FALSE, add=TRUE)

plot(ecdf(resAll[!is.na(theta) &  is.na(OUTRIDER) & !is.na(Kremer), meanCorrected]), col=colECDF[1], do.points=FALSE, 
     main='ECDF of gene mean expression', xlab='Gene mean expression', log='x', xlim=range(resAll[,meanCorrected], na.rm=TRUE))
plot(ecdf(resAll[!is.na(theta) & !is.na(OUTRIDER) &  is.na(Kremer), meanCorrected]), col=colECDF[2], do.points=FALSE, add=TRUE)
plot(ecdf(resAll[!is.na(theta) & !is.na(OUTRIDER) & !is.na(Kremer), meanCorrected]), col=colECDF[3],  do.points=FALSE, add=TRUE)

ks.test(resAll[!is.na(theta) & !is.na(OUTRIDER) &  is.na(Kremer), theta], resAll[!is.na(theta) & is.na(OUTRIDER) & !is.na(Kremer), theta])
ks.test(resAll[!is.na(theta) & !is.na(OUTRIDER) &  is.na(Kremer), meanCorrected], resAll[!is.na(theta) & is.na(OUTRIDER) & !is.na(Kremer), meanCorrected])

legend('topleft', names(colECDF), col=colECDF, pch=20, lty=1)


#'
#' # Overlap with Kremer et al 
#' 
resAll[unsolvedCase==TRUE,.(
    NumOverlap=sum(OUTRIDER == TRUE & Kremer == TRUE, na.rm=TRUE),
    PercentOverlap=round(sum(OUTRIDER == TRUE & Kremer == TRUE, na.rm=TRUE)/sum(Kremer == TRUE, na.rm=TRUE)*100, 1),
    NumOnlyOUTRIDER=sum(OUTRIDER == TRUE & is.na(Kremer), na.rm=TRUE),
    NumOnlyOUTRIDERDownReg=sum(OUTRIDER == TRUE & is.na(Kremer) & zScore_OUTRIDER < 0, na.rm=TRUE),
    PercentOnlyOUTRIDER=round(sum(OUTRIDER == TRUE & is.na(Kremer), na.rm=TRUE)/sum(OUTRIDER == TRUE, na.rm=TRUE)*100, 1),
    NumOutrider=sum(OUTRIDER == TRUE, na.rm=TRUE),
    NumPCA=sum(PCA == TRUE, na.rm=TRUE),
    NumPEER=sum(PEER == TRUE, na.rm=TRUE),
    NumPCAMoreOutrider=round(sum(PCA == TRUE, na.rm=TRUE)/sum(OUTRIDER == TRUE, na.rm=TRUE), 1),
    NumPEERMoreOutrider=round(sum(PEER == TRUE, na.rm=TRUE)/sum(OUTRIDER == TRUE, na.rm=TRUE), 1))]

resAll[DISEASE_CAUSING == TRUE, .(geneID, sampleID, RNA_ID, Kremer, OUTRIDER, PEER, PCA)]

#upload library
height=1000
width=1000
main.pos=c(0.5,1.0)
getVenPlot <- function(d2p=data2plot, main='Vendiagram',
                    clusterNames=c("Expression outliers by\nKremer et al.",
                            "Expression outliers by\nOUTRIDER", 
                            "Pathogenic by\nKremer et al."),
                    fill = RColorBrewer::brewer.pal(3, 'Dark2')[3:1],
                    cat.pos = c(200, 160, 0), ...,
                    cat.dist = c(0.0250, 0.0250, -0.04)){
    venn.diagram(
        category.names = clusterNames,
        fill = fill,
        cat.pos = cat.pos,
        cat.dist = cat.dist,
        main = main,
        main.cex = 3,
        main.fontface = 'bold',
        main.fontfamily = 'sans',
        main.just = c(1,0),
        main.pos = main.pos,
        x = d2p,
        filename=NULL,
        output = TRUE,
        imagetype="png",
        height = height, 
        width = width, 
        resolution = 1200,
        compression = "lzw",
        lwd = 2,
        lty = 'blank',
        cex = 2.5,
        fontface = "bold",
        fontfamily = "sans",
        cat.cex = 2.5,
        cat.fontface = "bold",
        cat.fontfamily = "sans",
        cat.default.pos = "outer",
        ...
    )
}

figFuns <- c()
figFuns['A'] <- c(function(){
    plotAsImage(height=height, width=width, mar=c(0,0,0,0)+0.1, expr={
        grid.draw(getVenPlot(filterdData2plot[c('Kremer', 'OUTRIDER', 'KNOWN')], 
                main="Filtered Kremer data set\nOverlap between Kremer and OUTRIDER",
                rotation = 1))
    })
})
figFuns['B'] <- c(function(){
    plotAsImage(height=height, width=width, mar=c(0,0,0,0)+0.1, expr={
        grid.draw(getVenPlot(filterdData2plot[c('PCA', 'PEER', 'OUTRIDER', 'KNOWN')],
                main="Filtered Kremer data set\nOverlap between PCA, PEER and OUTRIDER",
                clusterNames = c('PCA', 'PEER', 'OUTRIDER', 'Pathogenic by\nKremer et al.'), 
                fill=RColorBrewer::brewer.pal(4, 'Dark2')[4:1],
                cat.pos = c(-20,20,320,1), cat.dist = c(0.2,0.2,0.13,0.1)
                ))
    })
})


#########################
# create layout 
LAYOUT_MATRIX_FOR_FIGURE <- matrix(1:2, ncol=1, byrow=TRUE)
PENAL_CEX <- rep(1,length(figFuns))
names(PENAL_CEX) <- names(figFuns)

# Main figure
main.pos=c(0.5,1.0)
make_paper_figure(figFuns, LAYOUT_MATRIX_FOR_FIGURE, 2, penal_cex = PENAL_CEX, 
        figure_name=FIGURE_NAME, figdir=FIGDIR, upperLetter=TRUE, window_ratio = 0.5)


system(paste(convertCMD, snakemake@output$outPdf, snakemake@output$outPng))

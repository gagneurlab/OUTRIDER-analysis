#'---
#' title: Figure Supl 3 Q-Q plots for GTEx and Kremer
#' author: Christian Mertes
#' wb:
#'   input: 
#'     - gtexods:   '`sm config["DATADIR"] + "/fitOutrider/" + config["FIGURE_IMPLEMENTATION"] + "/Skin_Not_Sun_Exposed_Suprapubic_ODS.RDS"`' 
#'     - kremerods: '`sm config["DATADIR"] + "/fitOutrider/" + config["FIGURE_IMPLEMENTATION"] + "/Kremer_ODS.RDS"`'
#'   output:
#'     - gtexPNG: '`sm config["FIGDIR"] + "/FigureS06_qqplots_gtex.png"`'
#'     - gtexPDF: '`sm config["FIGDIR"] + "/FigureS06_qqplots_gtex.pdf"`'
#'     - kremerPNG: '`sm config["FIGDIR"] + "/FigureS07_qqplots_kremer.png"`'
#'     - kremerPDF: '`sm config["FIGDIR"] + "/FigureS07_qqplots_kremer.pdf"`'
#' output:
#'   html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---

if(FALSE){
    snakemake <- readRDS('tmp.snakemake.RDS'); snakemake 
}

source("./src/r/config.R")

###############################
#load data (2 OUTRIDER objects)
###############################
gtex_data <- readRDS(snakemake@input$gtexods)
kremer_data <- readRDS(snakemake@input$kremerods)

nc <- 6
nr <- 4 
ods <- kremer_data
name <- "Kremer"
RAND_SEED <- 42

plotData <- function(ods, name){
    set.seed(RAND_SEED)
    layout(widths = rep(c(2, 4), c(1, nc)), heights = rep(c(2,4,2), c(1,nr,1)),
        rbind(rbind(nc*nr+2, cbind(nc*nr+1, matrix(1:(nc*nr), ncol=nc))), nc*nr+3))
    ids <- sample(1:nrow(ods), nc*nr)
    for(i in ids){
        par(cex=1, mar=c(1,1,1,1))
        plotQQ(ods, i, main='', legendPos=NA, zScore=Z_LIMIT, padj=FDR_LIMIT)
    }
    plot.new()
    plot.new()
    plot.new()
    
    xlab <- expression(paste(-log[10], " (expected ", italic(P), "-value)"))
    ylab <- expression(paste(-log[10], " (observed ", italic(P), "-value)"))
    
    mtext(xlab, side=1, line = -1,   at = 0.55, cex=1.5)
    mtext(ylab, side=2, line = -1.5, at = 13-nc/2,   cex=1.5)
    mtext(paste(nc*nr, 'random Q-Q plots for', name), at=0.55, side=3, 
            line=54-nc*3.35, cex=2)
}

# plot it
data2plot <- list(GTEx = gtex_data, Kremer = kremer_data)
for(idx in seq_along(data2plot)){
    name <- names(data2plot)[idx]
    outFile <- snakemake@output[[paste0(tolower(name), 'PDF')]]
    fig_name <- gsub('.pdf$', '', basename(outFile))
    data <- data2plot[[idx]]
    
    paper_png(file.path(FIGDIR, fig_name), height=8, width=12, units = 'in')
        plotData(data, name)
        put_letter_topleft(LETTERS[idx], y_shift = -34, cex = 2)
    dev.off()
    paper_pdf(file.path(FIGDIR, fig_name), height=8, width=12)
        plotData(data, name)
        put_letter_topleft(LETTERS[idx], y_shift = -34, cex = 2)
    dev.off()
}



#'---
#' title: Fig S1
#' author: Felix Brechtmann
#' wb:
#'   params: 
#'     - methods: '`sm PAPER_METHODS`' 
#'   input: 
#'     - odsKremer: '`sm expand(config["DATADIR"] + "/fitOutrider/{method}/Kremer_ODS.RDS", method=PAPER_METHODS)`'
#'     - odsGTEx:   '`sm expand(config["DATADIR"] + "/fitOutrider/{method}/Skin_Not_Sun_Exposed_Suprapubic_ODS.RDS", method=PAPER_METHODS)`'
#'   output:
#'     - pdfKremer: '`sm config["FIGDIR"] + "/FigureS02_heatmaps_Kremer.pdf"`'
#'     - pngKremer: '`sm config["FIGDIR"] + "/FigureS02_heatmaps_Kremer.png"`'
#'     - pdfGTEx:   '`sm config["FIGDIR"] + "/FigureS03_heatmaps_GTEx.pdf"`'
#'     - pngGTEx:   '`sm config["FIGDIR"] + "/FigureS03_heatmaps_GTEx.png"`'
#' output:
#'   html_document:
#'     code_folding: show
#'     code_download: TRUE
#'---
#'

if(FALSE){
    snakemake <- readRDS('tmp.snakemake.RDS'); snakemake 
    
    odsFiles <- list(
        GTEx = list(
            PCA = 'Output/data/fitOutrider/pca/Skin_Not_Sun_Exposed_Suprapubic_ODS.RDS',
            PEER = 'Output/data/fitOutrider/peer/Skin_Not_Sun_Exposed_Suprapubic_ODS.RDS',
            OUTRIDER = 'Output/data/fitOutrider/NLas_TCNo/Skin_Not_Sun_Exposed_Suprapubic_ODS.RDS'
        ), 
        Kremer = list(
            PCA = 'Output/data/fitOutrider/pca/Kremer_ODS.RDS',
            PEER = 'Output/data/fitOutrider/peer/Kremer_ODS.RDS',
            OUTRIDER = 'Output/data/fitOutrider/NLas_TCNo/Kremer_ODS.RDS'
        )
    )
    
    pdfOutGTEx <- './Output/'
}


source('./src/r/config.R')
library(ComplexHeatmap)
library(circlize)

#' 
#' # Get input
#' 
odsFiles <- list(
    GTEx = snakemake@input$odsGTEx,
    Kremer = snakemake@input$odsKremer)
names(odsFiles[[1]]) <- c('PCA', 'PEER', 'OUTRIDER')
names(odsFiles[[2]]) <- c('PCA', 'PEER', 'OUTRIDER')
pdfOutGTEx <- snakemake@output$pdfGTEx
pngOutGTEx <- snakemake@output$pngGTEx
pdfOutKremer <- snakemake@output$pdfKremer
pngOutKremer <- snakemake@output$pngKremer
convertCMD <- 'convert -density 300 -trim -quality 100 -flatten'
odsFiles

#' 
#' # Functions
#' 
getCor <- function(ods, normalized=TRUE){
    lc <- log2(counts(ods, normalized=normalized) + 1)
    clc <- lc - rowMeans(lc)
    corc <- cor(clc, method="spearman")
    corc
}

getHeatmap <- function(corc, title){
        Heatmap(corc, 
            col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
            show_row_names=FALSE,
            show_column_names=FALSE,
            column_title='Samples',
            column_title_side = 'bottom',
            row_title='Samples',
            name=title
        )
}

plotLabel <- function(label){
    grid.text(label, 
            gp=gpar(fontface="bold", fontsize=20), 
            x=unit(0.03, "npc"), 
            y=unit(0.97, "npc"))
}

plotGrid <- function(h1, h2, h3, h4, labels=LETTERS[1:4]){
    
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nr=2, nc=2)))
    pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
    draw(h1, newpage=FALSE)
    plotLabel(labels[1])
    popViewport()
    
    pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
    draw(h2, newpage=FALSE)
    plotLabel(labels[2])
    popViewport()
    
    pushViewport(viewport(layout.pos.row=2, layout.pos.col=1))
    draw(h3, newpage=FALSE)
    plotLabel(labels[3])
    popViewport()
    
    pushViewport(viewport(layout.pos.row=2, layout.pos.col=2))
    draw(h4, newpage=FALSE)
    plotLabel(labels[4])
    popViewport()
}

#' 
#' # Kremer data set
#' 
odslsk <- lapply(odsFiles$Kremer, readRDS)

h1k <- getHeatmap(getCor(odslsk[['OUTRIDER']], FALSE), 'Kremer\nraw counts')
h2k <- getHeatmap(getCor(odslsk[['OUTRIDER']], TRUE), 'Kremer\nOUTRIDER')
h3k <- getHeatmap(getCor(odslsk[['PCA']], TRUE), 'Kremer\nPCA')
h4k <- getHeatmap(getCor(odslsk[['PEER']], TRUE), 'Kremer\nPEER')

plotGrid(h1k, h2k, h3k, h4k, LETTERS[1:4])

pdf(pdfOutKremer, width=10, height=8)
    plotGrid(h1k, h2k, h3k, h4k, LETTERS[1:4])
dev.off()
system(paste(convertCMD, pdfOutKremer, pngOutKremer))


#' 
#' # GTEx data set
#' 
odslsg <- lapply(odsFiles$GTEx, readRDS)

h1g <- getHeatmap(getCor(odslsg[['OUTRIDER']], FALSE), 'GTEx\nraw counts')
h2g <- getHeatmap(getCor(odslsg[['OUTRIDER']], TRUE), 'GTEx\nOUTRIDER')
h3g <- getHeatmap(getCor(odslsg[['PCA']], TRUE), 'GTEx\nPCA')
h4g <- getHeatmap(getCor(odslsg[['PEER']], TRUE), 'GTEx\nPEER')

plotGrid(h1g, h2g, h3g, h4g, LETTERS[1:4])

pdf(pdfOutGTEx, width=10, height=8)
    plotGrid(h1g, h2g, h3g, h4g, LETTERS[5:8])
dev.off()

system(paste(convertCMD, pdfOutGTEx, pngOutGTEx))

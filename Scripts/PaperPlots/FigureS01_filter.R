#'---
#' title: Fig S1
#' author: Felix Brechtmann
#' wb:
#'   input: 
#'     - ggplotKremer: '`sm config["DATADIR"] + "/filteredODS/Kremer_ggplots.RDS"`'
#'     - ggplotGTEx: '`sm config["DATADIR"] + "/filteredODS/Skin_Not_Sun_Exposed_Suprapubic_ggplots.RDS"`'
#'   output:
#'     - pdfout: '`sm config["FIGDIR"] + "/FigureS01_filtering.pdf"`'
#'     - pngout: '`sm config["FIGDIR"] + "/FigureS01_filtering.png"`'
#' output:
#'   html_document:
#'     code_folding: show
#'     code_download: TRUE
#'---
#'

if(FALSE){
    snakemake <- readRDS('tmp.snakemake.RDS'); snakemake 
    inputFiles <- c(
        Kremer = 'Output/data/filteredODS/Kremer_ggplots.RDS',
        GTEx = 'Output/data/filteredODS/Skin_Not_Sun_Exposed_Suprapubic_ggplots.RDS')
}

#' 
#' # Load input and config
#' 
source('src/r/config.R')

pdfOut <- snakemake@output$pdfout
pngOut <- snakemake@output$pngout

inputFiles <- c(
    Kremer = snakemake@input$ggplotKremer,
    GTEx = snakemake@input$ggplotGTEx)
ggplots <- lapply(inputFiles, readRDS)

#'
#' # Load figures
#' 
genefilterKremer <- ggplots$Kremer$ggfpkm + theme(legend.position = c(0.02, 0.8)) + 
    scale_y_continuous(labels = fancy_scientific) + 
    labs(title='Kremer')
genefilterGtex <- ggplots$GTEx$ggfpkm + theme(legend.position = c(0.02, 0.8)) +
    scale_y_continuous(labels = fancy_scientific) +
    labs(title='GTEx')

#'
#' # Assemble figure
#' 
gg <- plot_grid(genefilterGtex, genefilterKremer, labels = c("A", "B"), ncol=2)
gg

#' 
#' # Save plots
#' 
ggsave(filename=pdfOut, plot=gg, width=12, height=6, device='pdf')
ggsave(filename=pngOut, plot=gg, width=12, height=6, device='png', dpi=900)

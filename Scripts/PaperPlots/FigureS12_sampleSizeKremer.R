#'---
#' title: Figure S9 Sample size analysis for Kremer
#' author: Christian Mertes
#' 
#' wb:
#'   threads: 70
#'   params: 
#'     - methods: '`sm PAPER_METHODS`'
#'   input: 
#'     - rdsIn: '`sm config["DATADIR"] + "/singleParameterAnaylsis/SampleSizeKremer.RDS"`'
#'   output:
#'     - pdf: '`sm config["FIGDIR"] + "/FigureS12_sample_size_analysis.pdf"`'
#'     - png: '`sm config["FIGDIR"] + "/FigureS12_sample_size_analysis.png"`'
#' output:
#'   html_document:
#'     code_folding: show
#'     code_download: TRUE
#'---

if(FALSE){
    snakemake <- readRDS('tmp.snakemake.RDS'); snakemake
    resFile <- 'Output/data/singleParameterAnaylsis/SampleSizeKremer.RDS'
    methods <- c('pca', 'peer', CONFIG_YAML$FIGURE_IMPLEMENTATION)
}


#' 
#' # Read input and config
#' 
source('./src/r/config.R')

resFile <- snakemake@input$rdsIn
pdfOut <- snakemake@output$pdf
pngOut <- snakemake@output$png

resFile
pdfOut
pngOut

#' 
#' # Read data
#' 
res <- readRDS(resFile)
plotDt <- res$plotDt
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
        breaks=unique(plotDt$geneID),
        labels=sapply(unique(plotDt$geneID), function(xy) {
                parse(text=paste('italic(', xy, ')'))
        })) + 
    ggtitle('Sample size analysis (Kremer)')
gg
    

#' 
#' # Save figures as PNG and PDF
#' 
ggsave(pdfOut, gg, 'pdf', width=8, height=8)
ggsave(pngOut, gg, 'png', width=8, height=8, dpi=900)


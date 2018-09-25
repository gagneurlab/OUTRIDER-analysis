#'---
#' title: Figure2
#' author: Felix Brechtmann, Christian Mertes
#' wb:
#'   input: 
#'     - gtexods: '`sm config["DATADIR"] + "/fitOutrider/" + config["FIGURE_IMPLEMENTATION"] + "/Skin_Not_Sun_Exposed_Suprapubic_ODS.RDS"`'
#'   output:
#'     - figurePdf: '`sm config["FIGDIR"] + "/Figure2_nb_examples.pdf"`'
#'     - figurePng: '`sm config["FIGDIR"] + "/Figure2_nb_examples.png"`'
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---
#'

if(FALSE){
    gtexodsFile <- 'Output/data/fitOutrider/ed_NRob_NCR_NTC_YLAS/Skin_Not_Sun_Exposed_Suprapubic_ODS.RDS'
}

# Paper Figure 3 examples of a gene with and without ourliers.
#
# Script to generate the full figure including:
#
# - QQplots for the two cases
# - normalized expression for the same genes.
source("./src/r/config.R")

FIGURE_NAME <- "Figure2_nb_examples"
gtexodsFile <- snakemake@input$gtexods
ods <- readRDS(gtexodsFile)
################################
#load data (2 OUTRIDER objects)
###############################

normal_gene <- c(TRIM33='ENSG00000197323.6')
aberrant_gene <- c(SLC39A4='ENSG00000147804.5')
#ENSG00000077713.14'#"ENSG00000237973.1"# 'ENSG00000162572.15' 

### fig functions ###
letterCex   <- 1.9
upperLetter <- TRUE
PENAL_CEX   <- rep(1.3, 4)
mar_bottom  <- 5
mar_left    <- 5.5
mar_top     <- 2.5
mar_right   <- 0.5
ylab_left_line= 3.5
xlab_line= 2.5
main_line=0.5


figFuns <- c()

figFuns['a'] <- c(function(){
    par(mar=c(mar_bottom, mar_left, mar_top, mar_right))
    plotExpressionRank(ods, normal_gene, basePlot=TRUE, 
            padjCut=FDR_LIMIT, zScoreCut=Z_LIMIT,
            main=bquote(paste("Expression plot: ", italic(.(names(normal_gene))))))
})

figFuns['b'] <- c(function(){
    par(mar=c(mar_bottom, mar_left, mar_top, mar_right))
    plotQQ(ods, normal_gene, legendPos="topleft", padj=FDR_LIMIT,
            zScore=Z_LIMIT, main=bquote(paste("Q-Q plot: ", italic(.(names(normal_gene))))))
})

figFuns['c'] <- c(function(){
    par(mar=c(mar_bottom, mar_left, mar_top, mar_right))
    plotExpressionRank(ods, aberrant_gene, basePlot=TRUE, 
            padjCut=FDR_LIMIT, zScoreCut=Z_LIMIT,
            main=bquote(paste("Expression plot: ", italic(.(names(aberrant_gene))))))
})

figFuns['d'] <- c(function(){
    par(mar=c(mar_bottom, mar_left, mar_top, mar_right))
    plotQQ(ods, aberrant_gene, legendPos="topleft", padj=FDR_LIMIT,
           zScore=Z_LIMIT, main=bquote(paste("Q-Q plot: ", italic(.(names(aberrant_gene))))))
})
################################
#Output:
plotIt <- TRUE


#########################
# create layout 

names(figFuns) <- LETTERS[1:length(figFuns)]
LAYOUT_MATRIX_FOR_FIGURE <- matrix(1:4, nrow=2, byrow = TRUE)
names(PENAL_CEX) <- names(figFuns)

make_paper_figure(figFuns, LAYOUT_MATRIX_FOR_FIGURE, 1, penal_cex = PENAL_CEX, 
        FIGURE_NAME, FIGDIR, upperLetter=upperLetter, cex=letterCex)

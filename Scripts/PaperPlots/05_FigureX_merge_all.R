#'---
#' title: Merge Main and Supplement figures
#' author: Christian Mertes
#' wb:
#'  input: 
#'   - Figure1:   '`sm config["FIGDIR"] + "/Figure1_scheme.pdf"`'
#'   - Figure2:   '`sm config["FIGDIR"] + "/Figure2_nb_examples.pdf"`'
#'   - Figure3:   '`sm config["FIGDIR"] + "/Figure3_global_result.pdf"`'
#'   - Figure4:   '`sm config["FIGDIR"] + "/Figure4_outrider_recall.pdf"`'
#'   - Figure5:   '`sm config["FIGDIR"] + "/Figure5_gtex_enrichment.pdf"`'
#'   - FigureS1:  '`sm config["FIGDIR"] + "/FigureS01_filtering.pdf"`'
#'   - FigureS2:  '`sm config["FIGDIR"] + "/FigureS02_heatmaps_Kremer.pdf"`'
#'   - FigureS3:  '`sm config["FIGDIR"] + "/FigureS03_heatmaps_GTEx.pdf"`'
#'   - FigureS4:  '`sm config["FIGDIR"] + "/FigureS04_zscoreSearch.pdf"`'
#'   - FigureS5:  '`sm config["FIGDIR"] + "/FigureS05_mu_plot.pdf"`'
#'   - FigureS6:  '`sm config["FIGDIR"] + "/FigureS06_qqplots_gtex.pdf"`'
#'   - FigureS7:  '`sm config["FIGDIR"] + "/FigureS07_qqplots_kremer.pdf"`'
#'   - FigureS8:  '`sm config["FIGDIR"] + "/FigureS08_recall_gtex.pdf"`'
#'   - FigureS9:  '`sm config["FIGDIR"] + "/FigureS09_recall_Kremer.pdf"`'
#'   - FigureS10: '`sm config["FIGDIR"] + "/FigureS10_recall_single_bin.pdf"`'
#'   - FigureS11: '`sm config["FIGDIR"] + "/FigureS11_gtex_enrichment_all.pdf"`'
#'   - FigureS12: '`sm config["FIGDIR"] + "/FigureS12_sample_size_analysis.pdf"`'
#'   - FigureS13: '`sm config["FIGDIR"] + "/FigureS13_Kremer_benchmark.pdf"`'
#'  output:
#'  - figure:  '`sm config["FIGDIR"] + "/Figure_all_main.pdf"`'
#'  - figureS: '`sm config["FIGDIR"] + "/Figure_all_supl.pdf"`'
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---
source("./src/r/config.R")


# source config and load package
mainFigures <- unlist(unique(snakemake@input[grepl("Figure[0-9]",  snakemake@input)]))
suplFigures <- unlist(unique(snakemake@input[grepl("FigureS[0-9]", snakemake@input)]))


message(date(), ': Start with main figures ... ')
system(paste0(
    'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dNumRenderingThreads=10 -dNOGC ',
    ' -dBandBufferSpace=5000000000 -dBufferSpace=10000000000 -sBandListStorage=memory ',
    '-sOutputFile=', snakemake@output$figure, ' ',
    paste(mainFigures, collapse = ' ')))
system(paste0(
    'gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/prepress -dNOPAUSE ',
    '-dQUIET -dBATCH -sOutputFile=', 
    gsub(".pdf$", "_reduced.pdf", snakemake@output$figure), ' ', snakemake@output$figure))

message(date(), ': Start with supplement figures ... ')
system(paste0(
    'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dNumRenderingThreads=10 -dNOGC ',
    ' -dBandBufferSpace=5000000000 -dBufferSpace=10000000000 -sBandListStorage=memory ',
    '-sOutputFile=', snakemake@output$figureS, ' ',
    paste(suplFigures, collapse = ' ')))
system(paste0(
    'gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/prepress -dNOPAUSE ',
    '-dQUIET -dBATCH -sOutputFile=', 
    gsub(".pdf$", "_reduced.pdf", snakemake@output$figureS), ' ', snakemake@output$figureS))




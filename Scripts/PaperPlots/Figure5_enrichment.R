#'---
#' title: Figure 5 GTEX enrichment
#' author: Christian Mertes
#' wb:
#'   params:
#'     - methods: '`sm PAPER_METHODS`'
#'   input: 
#'     - rdsFiles: '`sm expand(config["DATADIR"] + "/GTEx_variant_enrichment/{dataset}_featureset.RDS", dataset=config["GTEx_tissues"])`'
#'     - odsFiles: '`sm expand(config["DATADIR"] + "/fitOutrider/{method}/{dataset}_ODS.RDS", method=PAPER_METHODS, dataset=config["GTEx_tissues"])`'
#'   output: 
#'     - pdf: '`sm config["FIGDIR"] + "/Figure5_gtex_enrichment.pdf"`'
#'     - png: '`sm config["FIGDIR"] + "/Figure5_gtex_enrichment.png"`'
#' output:
#'   html_document:
#'     code_folding: show
#'     code_download: TRUE
#'---

# Debug data
if(FALSE){
    snakemake <- readRDS('tmp.snakemake.RDS'); snakemake
    
    rdsFiles <- list.files('./Output/data/GTEx_variant_enrichment/', full.names = TRUE, pattern = '_featureset.RDS')
    mainMethod <- CONFIG_YAML$FIGURE_IMPLEMENTATION
    methods4Plotting <- c('pca', 'peer', mainMethod)
    outPdf <- 'tmp.pdf'
}

# load config
source('src/r/config.R')
register(MulticoreParam(60, progressbar = TRUE))


#' # Extract data
#' 
#' ## Input
#' 
rdsFiles <- snakemake@input$rdsFiles
mainMethod <- snakemake@config$FIGURE_IMPLEMENTATION
methods4Plotting <- snakemake@params$methods
outPdf <- snakemake@output$pdf
outPng <- snakemake@output$png

datasets <- gsub('_featureset.RDS', '', basename(rdsFiles))
names(rdsFiles) <- datasets
rdsFiles
length(rdsFiles)

#' 
#' ## Read data
#' 
cutoffList <- list(c1=c(0.005, 2), c2=c(0.0005, 3), c3=c(0.00005, 5))
res <- rbindlist(lapply(cutoffList, function(curCut){
    x <- names(rdsFiles)
    rbindlist(bplapply(x, getEnrichmentForTissues, rdsFiles=rdsFiles, cutoffs=curCut))
}))
res

#' Set correct boundaries 
res <- res[!grepl('RAND', Method)]
res[,enrichNAs:=sum(is.na(enrichment)), by=Tissue]
res[enrichNAs > 1, unique(Tissue)]


newres <- merge(res[Method != mainMethod], res[Method == mainMethod],
                by=c('Type', 'Cutoff', 'Tissue'))

setnames(newres, 'enrichment.x', 'Enrichment')
setnames(newres, 'enrichment.y', 'AE_Enrichment')
setnames(newres, 'Method.x', 'Method')
newres[,Type:=gsub('([><])', ' \\1 ', gsub('(-value| score) \\(|\\)', '', Type))]

newres <- newres[Method %in% methods4Plotting]
newres <- renameCorrectionMethods(newres, 'Method')

#'
#' # Scatter plots of enrichments per method and cutoff
#' 
#' ## Pvalue plot
#' 
#+ fig.width=14, fig.height=14
gg <- ggplot(newres[grepl('^P', Type)], aes(Enrichment, AE_Enrichment, col=Method)) + geom_point() + 
    geom_abline(intercept = 0,slope = 1, col='gray') + 
    grids(linetype = 'dotted') +
    scale_y_log10() + 
    scale_x_log10() + 
    scale_color_manual(values=brewer.pal(length(methods4Plotting) + 1, 'Dark2')[1:length(methods4Plotting)+1]) + 
    facet_grid(Type ~ Method) + 
    geom_smooth(method='lm') + 
    theme(strip.text.y = element_text(margin = margin(0,0.2,0,0.2, "cm")),
          strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"))) + 
    labs(y=paste('OUTRIDER enrichment')) + 
    theme(legend.position='none')
gg


#' 
#' # Save figures as PNG and PDF
#' 
ggsave(outPdf, gg, 'pdf', width=6, height=9)
ggsave(outPng, gg, 'png', width=6, height=9, dpi=900)

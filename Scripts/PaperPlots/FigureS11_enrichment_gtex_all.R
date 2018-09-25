#'---
#' title: Figure 5 GTEX enrichment
#' author: Christian Mertes
#' wb:
#'   input: 
#'     - rdsFiles: '`sm expand(config["DATADIR"] + "/GTEx_variant_enrichment/{dataset}_featureset.RDS", dataset=config["GTEx_tissues"])`'
#'   output: 
#'     - pdf: '`sm config["FIGDIR"] + "/FigureS11_gtex_enrichment_all.pdf"`'
#'     - png: '`sm config["FIGDIR"] + "/FigureS11_gtex_enrichment_all.png"`'
#' output:
#'   html_document:
#'     code_folding: show
#'     code_download: TRUE
#'---

# Debug data
if(FALSE){
    snakemake <- readRDS('tmp.snakemake.RDS'); snakemake
    rdsFiles <- list.files('./Output/data/GTEx_variant_enrichment/', full.names = TRUE, pattern = '_featureset.RDS')
    mainMethod <- 'ed_NRob_NCR_NTC_YLAS'
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
newres[,Type:=gsub('Z', '|Z|', Type)]

newres <- newres[Method %in% c('pca', 'peer', 'LiZscore')]
newres[,Method:=gsub('pca', 'PCA', Method)]
newres[,Method:=gsub('peer', 'PEER', Method)]
newres[,Method:=gsub('LiZscore', 'Li et al.', Method)]

#'
#' # Scatter plots of enrichments per method and cutoff
#' 
#' ## Pvalue plot
#' 
#+ fig.width=14, fig.height=14
gg <- ggplot(newres, aes(Enrichment, AE_Enrichment, col=Method)) + 
    geom_point() + 
    geom_abline(intercept = 0,slope = 1, col='gray') + 
    grids(linetype = 'dotted') +
    scale_y_log10() + 
    scale_x_log10() + 
    scale_color_manual(values=RColorBrewer::brewer.pal(4, 'Dark2')[c(4,2:3)]) + 
    facet_wrap(facets = c('Method', 'Type'), scales='free', ncol=3) + 
    geom_smooth(method='lm') + 
    theme(strip.text.y = element_text(margin = margin(0,0.2,0,0.2, "cm")),
          strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"))) + 
    labs(title='SNP enrichment in GTEx', y=paste('OUTRIDER enrichment'))
gg


#' 
#' # Save figures as PNG and PDF
#' 
ggsave(snakemake@output$pdf, gg, 'pdf', width=9, height=11)
ggsave(snakemake@output$png, gg, 'png', width=9, height=11, dpi=900)

#'---
#' title: Kremer Data Set
#' author: Felix Brechtmann
#' wb:
#'  input: 
#'    - kremer: '`sm config["KREMER_COUNTS"]`'
#'    - gtfFile: '`sm config["KREMER_GENE_ANNOTATION_FILE"]`'
#'    - anno: '`sm config["SAMPLE_ANNOTATION_FILE"]`'
#'  output:
#'    - ods: '`sm config["DATADIR"] + "/rawODS/Kremer_ODS.RDS"`'
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---

if(FALSE){
    snakemake <- readRDS('tmp.snakemake.RDS'); snakemake
    countFile <- CONFIG_YAML$KREMER_COUNTS
    gtfFile <- CONFIG_YAML$KREMER_COUNTS
    annoFile <- CONFIG_YAML$SAMPLE_ANNOTATION_FILE
}

source('src/r/config.R')

countFile <- snakemake@input$kremer
gtfFile   <- snakemake@input$gtfFile
annoFile  <- snakemake@input$anno
odsOut    <- snakemake@output$ods

#' Load read counts
countFile
countdata <- read.table(countFile, check.names=FALSE, sep='\t')
dim(countdata)

#' Read Sample Annotation
annoFile
sampleAnno <- fread(annoFile)
sa <- as.data.frame(sampleAnno, check.names=FALSE)
dim(sa)

common_samples <- intersect(sa$RNA_ID, colnames(countdata))
sa <- sa[sa$RNA_ID %in% common_samples,]

#' Create Outrider data set.
ods <- OutriderDataSet(countData=countdata)
colData(ods) <- cbind(colData(ods), sa)

#' Import GTF and filter by fpkm.
library(RMySQL)
library(AnnotationDbi)
txdb <- loadDb(gtfFile)
con <- dbConnect(MySQL(), host='genome-mysql.cse.ucsc.edu', dbname="hg19", user='genome')
map <- dbGetQuery(con, 'select kgId AS TXNAME, geneSymbol from kgXref')

# Filter, storing the fpkm values and not subsetting the ods object.
ods <- computeGeneLength(ods, gtfFile=txdb, mapping=map)

#' Summary of filtered data and plot.
print(paste0('Number of samples: ', dim(ods)[2]))
print(paste0('Number of genes: ', dim(ods)[1]))
print(paste0('Number of genes after filtering: ', sum(mcols(ods)[['passedFilter']])))


#save Outrider data set to file.
saveRDS(ods, odsOut)


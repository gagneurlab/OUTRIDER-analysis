#'---
#' title: Create GTEx tissue object
#' author: Felix Brechtmann, Christian Mertes
#' wb:
#'  input:
#'    - counts:     '`sm config["GTEX_COUNTS"]`'
#'    - annoSample: '`sm config["GTEX_SAMPLE_ANNOTATION_FILE"]`'
#'    - annoPheno:  '`sm config["GTEX_PHENO_ANNOTATION_FILE"]`'
#'    - gff:        '`sm config["GTEX_GENE_ANNOTATION_FILE"]`'
#'  output:
#'    - ods: '`sm config["DATADIR"] + "/rawODS/{tissue}_ODS.RDS"`'
#'    - wBhtml: 'Output/html/GTEx_tissues/rawODS/{tissue}.html'
#'  type: noindex
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---

#
# DEBUG values
# 
tissue    <- "Skin_Not_Sun_Exposed_Suprapubic"
countFile <- "Data/rawdata/V6P/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz"
annoFile  <- "Data/rawdata/V6P/GTEx_Data_V6_Annotations_SampleAttributesDS.txt"
phenoFile <- "Data/rawdata/V6P/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt"
gffFile   <- "Data/rawdata/V6P/gencode.v19.genes.v6p_model.patched_contigs.gtf.gz"
countFile <- "Data/rawdata/V7/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.gz"
annoFile  <- "Data/rawdata/V7/GTEx_v7_Annotations_SampleAttributesDS.txt"
phenoFile <- "Data/rawdata/V7/GTEx_v7_Annotations_SubjectPhenotypesDS.txt"
gffFile   <- "Data/rawdata/V7/gencode.v19.genes.v7.patched_contigs.gtf"

source('src/r/config.R')

#'
#' Read input
#' 
tissue    <- snakemake@wildcards$tissue
countFile <- snakemake@input$counts
annoFile  <- snakemake@input$annoSample
phenoFile <- snakemake@input$annoPheno
gffFile   <- snakemake@input$gff
outFile   <- snakemake@output$ods

# Multicore 
register(MulticoreParam(5, progressbar=TRUE))

# load gtex
gtexdt <- fread(annoFile)
gtexdt <- gtexdt[!grepl('^K-562', SAMPID)]
gtexdt[,SMTSD:=gsub("_$", "", gsub("[\\s-()]+", "_", SMTSD, perl=TRUE))]
gtexdt[,subjectID:=gsub("^([^-]+-[^-]+)-.*", "\\1", SAMPID, perl=TRUE)]
dim(gtexdt)

#' 
#' # Load countdata from gtex
#' 
countFile
countreadCMD <- paste('zcat ', countFile)
if(Sys.info()['sysname'] == 'Darwin'){
    countreadCMD <- paste('zcat < ', countFile)
}
countdata <- fread(countreadCMD)
gtexdt[, ctsExists:=SAMPID %in% colnames(countdata)]
dim(countdata)
countdata[1:5,1:5]

#' 
#' ## Select tissue `{r snakemake@wildcards$tissue}`
#' 
colData <- gtexdt[SMTSD == tissue & ctsExists == TRUE]
colData <- colData[!duplicated(subjectID)]
dim(colData)

#' 
#' ## Clean sample data
#' 
hist(colData$SMRIN, breaks=100)
abline(v=5.7, col='red')
colData <- colData[SMRIN >= 5.7]
dim(colData)

table(colData$SMAFRZE)
colData <- colData[SMAFRZE %in% c('USE ME', 'RNASEQ')]
dim(colData)

table(colData$SMATSSCR)
colData <- colData[is.na(SMATSSCR) | SMATSSCR < 3]
dim(colData)

#' 
#' ## Add Phenotype data
#' 
phenoAnno <- fread(phenoFile)
colData <- merge(colData, phenoAnno, by.x="subjectID", by.y="SUBJID", all.x=TRUE)
if('GENDER' %in% colnames(colData)){
   setnames(colData, 'GENDER', 'SEX') 
}

#'
#' ## Final colData object
#' 
tissue
dim(colData)
DT::datatable(colData)

#'
#' ## Create OUTRIDER object
#' 
sampleIdx <- colnames(countdata) %in% colData[, SAMPID]
cts <- countdata[, sampleIdx, with=FALSE]

ods <- OutriderDataSet(countData=cts, colData=colData)
rownames(ods) <- countdata[,Name]
mcols(ods)$gene_symbol <- countdata[,Description]
ods


#' 
#' ## Annotate genes (basepair length)
#' 
ods <- computeGeneLength(ods, gffFile)


# 
#save Outrider data set to file.
outFile
saveRDS(ods, outFile)


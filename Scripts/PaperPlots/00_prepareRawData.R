#'---
#' title: Run OUTRIDER
#' author: Felix Brechtmann, Christian Mertes
#' wb:
#'  output:
#'    - KREMER_COUNTS:               '`sm config["KREMER_COUNTS"]`'
#'    - KREMER_GENE_ANNOTATION_FILE: '`sm config["KREMER_GENE_ANNOTATION_FILE"]`'
#'    - SAMPLE_ANNOTATION_FILE:      '`sm config["SAMPLE_ANNOTATION_FILE"]`'
#'    - GTEX_COUNTS:                 '`sm config["GTEX_COUNTS"]`'
#'    - GTEX_GENE_ANNOTATION_FILE:   '`sm config["GTEX_GENE_ANNOTATION_FILE"]`'
#'    - GTEX_SAMPLE_ANNOTATION_FILE: '`sm config["GTEX_SAMPLE_ANNOTATION_FILE"]`'
#'    - GTEX_PHENO_ANNOTATION_FILE:  '`sm config["GTEX_PHENO_ANNOTATION_FILE"]`'
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---
#'

source('src/r/config.R')

GTEX_VERSION <- basename(dirname(snakemake@output$GTEX_COUNTS))

#' 
#' Download raw data files for the full analysis
#' 

GTEX_URLS <- list(
    V7 = list(
        cts   = "https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.gz",
        gtf   = "https://storage.googleapis.com/gtex_analysis_v7/reference/gencode.v19.genes.v7.patched_contigs.gtf",
        anno  = "https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SampleAttributesDS.txt",
        pheno = "https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SubjectPhenotypesDS.txt"
    ),
    V6P = list(
        cts   = "https://storage.googleapis.com/gtex_analysis_v6p/rna_seq_data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz",
        gtf   = "https://storage.googleapis.com/gtex_analysis_v6p/reference/gencode.v19.genes.v6p_model.patched_contigs.gtf.gz",
        anno  = "https://storage.googleapis.com/gtex_analysis_v6p/annotations/GTEx_Data_V6_Annotations_SampleAttributesDS.txt",
        pheno = "https://storage.googleapis.com/gtex_analysis_v6p/annotations/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt"
    ),
    V6 = list(
        cts = "https://storage.googleapis.com/gtex_analysis_v6/rna_seq_data/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz",
        gtf  = "https://storage.googleapis.com/gtex_analysis_v6/reference/gencode.v19.genes.patched_contigs.gtf.gz",
        anno = "https://storage.googleapis.com/gtex_analysis_v6/annotations/GTEx_Data_V6_Annotations_SampleAttributesDS.txt",
        pheno = "https://storage.googleapis.com/gtex_analysis_v6/annotations/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt"
    ),
    inhouse = list(
        cts = "/s/project/sra-download/files/dbGaP-11206/phe000006.v2.GTEx_RNAseq.expression-data-matrixfmt.c1/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_reads.gct.gz",
        gtf  = "",
        anno = "/data/nasif12/home_if12/mertes/projects/OUTRIDER-analysis-final/Data/GTEx_SRA_anno_file.tsv",
        pheno = "https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SubjectPhenotypesDS.txt"
    )
)

URLS <- c(
    GTEX_COUNTS                 = GTEX_URLS[[GTEX_VERSION]][['cts']], 
    GTEX_GENE_ANNOTATION_FILE   = GTEX_URLS[[GTEX_VERSION]][['gtf']],
    GTEX_SAMPLE_ANNOTATION_FILE = GTEX_URLS[[GTEX_VERSION]][['anno']],
    GTEX_PHENO_ANNOTATION_FILE  = GTEX_URLS[[GTEX_VERSION]][['pheno']],
    KREMER_COUNTS               = 'https://media.nature.com/original/nature-assets/ncomms/2017/170612/ncomms15824/extref/ncomms15824-s1.txt',
    SAMPLE_ANNOTATION_FILE      = 'https://media.nature.com/original/nature-assets/ncomms/2017/170612/ncomms15824/extref/ncomms15824-s8.txt',
    KREMER_GENE_ANNOTATION_FILE = 'https://i12g-gagneurweb.in.tum.de/public/paper/mitoMultiOmics/ucsc.knownGenes.db')

FILE_NAMES <- snakemake@output[names(URLS)]

for(i in names(URLS)){
    file2save <- FILE_NAMES[[i]]
    message('Downloading file:', i)
    message('From: ', URLS[i])
    message('To: ', file2save)
    download.file(URLS[i], destfile=file2save, method='wget', extra='-c')
}

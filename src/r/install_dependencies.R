#'
#' Install requirements for the OUTRIDER analysis pipeline
#'
threads <- 4
reqirements <- c("png", "AnnotationDbi", "BBmisc", "BiocParallel",
        "BiocStyle", "circlize", "ComplexHeatmap", "DESeq2", "ensemblVEP", 
	"GenomicFeatures", "ggbeeswarm", "knitrBootstrap", "matrixStats",
	"plotly", "microbenchmark", "R.utils", "RColorBrewer", "reshape2",
        "RMySQL", "SRAdb", "VennDiagram", "cowplot", "data.table", "devtools",
        "edgeR", "ggplot2", "ggpubr", "ggthemes", "gplots", "gsubfn", "knitr",
        "magrittr", "plyr", "readr", "reticulate", 
        "rmarkdown", "scales", "stringr", "sva", "tidyr", "VariantAnnotation",
	"yaml", "zebrafishRNASeq"
)
githubReq <- c(patchwork="thomasp85/patchwork", OUTRIDER='gagneurlab/OUTRIDER')
sourcPkg <- c(peer='https://github.com/downloads/PMBio/peer/R_peer_source_1.3.tgz')

#' 
#' Setup bioconductor
#' 
instPkg <- installed.packages()[,'Package']
if(!'BiocInstaller' %in% instPkg){
    source("https://bioconductor.org/biocLite.R")
    biocLite("BiocInstaller")
}
library(BiocInstaller)

#' install missing packages
pkg2beInstalled <- reqirements[!reqirements %in% instPkg]
if(length(pkg2beInstalled) > 0){
    biocLite(pkg2beInstalled, Ncpu=threads)
}

#' install remaining packages from github
library(devtools)
pkg2beInstalled <- githubReq[!names(githubReq) %in% instPkg]
if(length(pkg2beInstalled) > 0){
    install_github(pkg2beInstalled)
}

#' install remaining packages from source 
pkg2beInstalled <- sourcPkg[!names(sourcPkg) %in% instPkg]
if(length(pkg2beInstalled) > 0){
    try(install_url(pkg2beInstalled))
}


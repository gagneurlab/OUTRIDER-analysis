##--------------------------------------------
## Rscript
## 

##--------------------------------------------
message("Loading required packages ...", appendLF=FALSE)
suppressPackageStartupMessages({
    ## Knitr config
    library(knitr)
    library(rmarkdown)
    library(GenomicFeatures)
    library(BiocStyle)
    opts_knit$set(root.dir=getwd())
    opts_chunk$set(cache=F)
    
    ##
    library(R.utils)
    library(yaml)
    
    library(plyr)
    
    ## Libraries for Analysis
    library(BiocParallel)
    library(cowplot)
    library(data.table)
    library(devtools)
    library(gplots)
    library(ggplot2)
    library(ggpubr)
    library(ggthemes)
    library(magrittr)
    library(matrixStats)
    library(RColorBrewer)
    library(reticulate)
    library(readr)
    library(tidyr)
    library(reshape2)
})
message(" DONE")

if(file.exists("../OUTRIDER/DESCRIPTION")){
    suppressPackageStartupMessages(devtools::load_all("../OUTRIDER/"))
} else if(!require(OUTRIDER)){
    devtools::install_github("gagneurlab/OUTRIDER")
    message("Loading OUTRIDER")
    suppressPackageStartupMessages(library(OUTRIDER))
}

CONFIG_YAML <- yaml.load_file("./wbuild.yaml")

## SMALL FUNCTIONS
## 
sourceDirectory('src/r/paperPlots')
sourceDirectory('src/r/helperFunction')
if(file.exists(".wBuild/wBuildParser.R")){
    source(".wBuild/wBuildParser.R")
}

# create pseudo snakemake object
if(!exists('snakemake')){
    setClass("snakemake",
            slots = list(config="list"),
            prototype = list(config=CONFIG_YAML)
    )
    snakemake <- new("snakemake")
} else {
    # use snakemake parameters if present
    saveRDS(snakemake, 'tmp.snakemake.RDS')
}

CONFIG_YAML <- snakemake@config
if("CONFIG_YAML" %in% search()){
    detach("CONFIG_YAML")
}
attach(CONFIG_YAML)


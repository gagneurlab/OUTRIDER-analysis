if(FALSE){
    snakemake <- readRDS('tmp.snakemake.RDS'); snakemake
    in_inj_counts <- "Output/data/RecallAnalysis/nsamples=All/inj=both/zscore=3/Simulation_injectedcounts.RDS"
    bestQFile     <- "Output/data/GTEx_not_sun_exposed_OutriderDONE.RDS"
    out_ods_obj   <- "Output/data/RecallAnalysis/nsamples=All/inj=both/zscore=3/correction=baseCooks/q=best/GTEx_not_sun_exposed_OUTRIDERfit.RDS"
    correction    <- "baseCooks"
    q4fit         <- "best"
    threads       <- 10
}

source("./src/r/config.R")

# inputs:
in_inj_counts <- snakemake@input[['inj_counts']]
bestQFile <-  snakemake@input[['bestQ']]
kremerCoVars <- fread(snakemake@input$KremerCovars)
gtexCoVars <- fread(snakemake@input$GTExCovars)

# outputs:
out_ods_obj <- snakemake@output[['outrider_fit']]
# ------

# Params:
correction <- snakemake@wildcards[['correction']]
q4fit <-  snakemake@wildcards[['q']]
dataset <- snakemake@wildcards$dataset

# Run params:
threads <- snakemake@threads

# Inferred params
# --------------------------------------------
autoCorrect <- TRUE
if(correction=='None'){
    autoCorrect <- FALSE
    correction <- NULL
}

autoCorrect
correction

register(MulticoreParam(threads))


#'
#' # Read data
#' Read count table.
ods <- readRDS(in_inj_counts)

#' Get best q
bestQ <- getBestQ(readRDS(bestQFile))
if(q4fit != "best"){
    bestQ <- as.integer(q4fit)
}
bestQ <- min(bestQ, dim(ods)-1)
bestQ

if(is.null(rownames(ods))){
    rownames(ods) <- paste("gene", 1:nrow(ods))
}

getDdsFit <- function(ods, dataset){
    library(DESeq2)
    dds <- DESeqDataSetFromMatrix(counts(ods), colData = colData(ods), ~1)
    if(dataset == 'Kremer'){
        colData(dds)$SEX <- factor(colData(dds)$SEX)
        colData(dds)$RNA_HOX_GROUP <- factor(colData(dds)$RNA_HOX_GROUP)
        design(dds) <- ~ 1 + SEX + RNA_HOX_GROUP
    } else if(dataset %in% CONFIG_YAML$GTEx_tissues){
        colData(dds)$SEX <- factor(colData(dds)$SEX)
        colData(dds)$AGE <- factor(gsub('-', '_', colData(dds)$AGE))
        colData(dds)$DTHHRDY <- factor(colData(dds)$DTHHRDY)
        design(dds) <- ~ 1 + SEX + AGE + DTHHRDY
    }
    dds <- estimateSizeFactors(dds)
    dds <- estimateDispersions(dds)
    dds    
}

#'
#'  # Run OUTRIDER 
if(correction == "baseCooks"){
    message(date(), ": Use baseCooks ranking")
    dds <- getDdsFit(ods, dataset)
    dds <- nbinomWaldTest(dds)
    
    cooksD <- assay(dds, 'cooks')
    assay(ods, 'cooks') <- cooksD
    
} else if(correction == "basePearsonRes"){
    message(date(), ": Use basePearson residual ranking")
    dds <- getDdsFit(ods, dataset)
    
    # DESeq2:::calculateCooksDistance
    modelMatrix <- DESeq2:::getModelMatrix(dds)
    dispersions <- DESeq2:::robustMethodOfMomentsDisp(dds, modelMatrix)
    V <- assay(dds, "mu") + dispersions * assay(dds, "mu")^2
    
    # https://support.bioconductor.org/p/57617
    PearsonRes <- (counts(dds) - assay(dds, "mu"))/sqrt(V)
    assay(ods, 'pearsonRes') <- PearsonRes
    
} else {
    message(date(), ": Use standard OUTRIDER pipeline")
    ods <- OUTRIDER(ods, q=bestQ, autoCorrect=autoCorrect, 
            implementation=correction)
}

saveRDS(ods, out_ods_obj)

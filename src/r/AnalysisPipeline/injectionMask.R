if(FALSE){
    snakemake <- readRDS('tmp.snakemake.RDS'); snakemake
}

source("./src/r/config.R")

# --------------------------------------------
# Snakemake
in_ods <- snakemake@input[['counts']]
out_true_outliers <- snakemake@output[['true_outliers']]


### parameters
inj <- snakemake@wildcards[['inj']]
Nsamples <- snakemake@wildcards[['nsamples']]
p = snakemake@config$P_BENCH
FDR_LIMIT = snakemake@config$FDR_LIMIT
dataset <- snakemake@wildcards[['dataset']]

# --------------------------------------------
# load samples with the least aberrant genes
ods <- readRDS(in_ods)

#'
#' Remove aberrant samples
numAberrant <- aberrant(ods, padj=FDR_LIMIT, zScore=Z_LIMIT, by='sample')
least_aberrant <- numAberrant < length(ods) * snakemake@config$OUTLIER_RATIO
if(grepl('old250', dataset)){
    least_aberrant <- !colnames(ods) %in% paste0('GTEX-', c(
            "11DZ1-2426-SM-5GZZX", "11WQK-1026-SM-5EQLX", "12584-0726-SM-5FQTK",
            "12WS9-0326-SM-59HJV", "12WSN-1426-SM-5GCO6", "1313W-0626-SM-5EQ4H", 
            "13OVI-2726-SM-5KM56", "14753-0726-SM-5QGQO", "V1D1-2026-SM-3GAF4", 
            "VJYA-1526-SM-3GIJV",  "Y5V6-2326-SM-4VDSA",  "YF7O-1526-SM-5IFI6", 
            "ZAB5-2626-SM-5KM3Y",  "ZDYS-2426-SM-4WKGI",  "ZQG8-1926-SM-4YCEO", 
            "ZV7C-0926-SM-59HL1",  "ZY6K-1426-SM-5GZX2",  "ZYFD-0626-SM-5E44E", 
            "ZYVF-0926-SM-5E44J"))
}
ods <- ods[,least_aberrant]

#' 
#' Down sample if requested
## Option to downsample
if(!Nsamples == 'All'){
    maxSamples <- min(ncol(ods), as.integer(Nsamples))
    if(maxSamples < ncol(ods)){
        ods <- ods[,sample(1:ncol(ods), maxSamples)]
    }
}
dim(ods)

## Filter genes with aberrant counts.
#ods <- ods[aberrant(ods, padj=0.01, by='gene')>0,]

rawtable <- counts(ods)

#inject outliers
size <- dim(rawtable)[1] * (dim(rawtable)[2])

index_mat <- matrix(nrow=nrow(rawtable), 
        data=sample(c(0,1,-1), size, prob=c(1 - p, p/2, p/2), replace=TRUE))

# create outrider object
mask_ods <- OutriderDataSet(countData=counts(ods), colData=colData(ods))
sizeFactors(mask_ods) <- NULL
thetaCorrection(mask_ods) <- NULL
metadata(mask_ods)$optimalEncDim <- metadata(ods)$optimalEncDim
metadata(mask_ods)$encDimTable <- metadata(ods)$encDimTable
metadata(mask_ods)$dim <- dim(mask_ods)

# use simulated values if present
if(all(c("true_mean", "true_sd") %in% assayNames(ods))){
    colData(mask_ods)[['true_sizeFactor']] <- colData(ods)[['true_sizeFactor']]
    assay(mask_ods, 'true_sd') <- assay(ods, "true_sd")
    assay(mask_ods, 'true_mean') <- assay(ods, "true_mean")
    if('Ncov' %in% colnames(mcols(ods))){
        mcols(mask_ods)[['Ncov']] <- mcols(ods)[['Ncov']]
    }
}

counts(mask_ods) <- rawtable
assay(mask_ods, "inj_mask") <- index_mat

# save it
saveRDS(mask_ods, out_true_outliers)

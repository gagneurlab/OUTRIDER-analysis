if(FALSE){
    snakemake <- readRDS('tmp.snakemake.RDS'); snakemake
}

source("./src/r/config.R")

# --------------------------------------------
# Snakemake
in_true_outliers <- snakemake@input[['true_outliers']]
out_inj_counts <- snakemake@output[['inj_counts']]

### parameters
inj_value <- snakemake@wildcards[['inj_value']]
inj <- snakemake@wildcards[['inj']]
Nsamples <- snakemake@wildcards[['nsamples']]
p = snakemake@config$P_BENCH
dataset <- snakemake@wildcards[['dataset']]


# --------------------------------------------
# Read data

inj_value <- unlist(strsplit(inj_value, '\\='))
zScore <- as.numeric(inj_value[2])
ods <- readRDS(in_true_outliers)
k <- counts(ods)
index_mat <- assay(ods, "inj_mask")
size <- nrow(k) * ncol(k)

zScore
size
ods

# --------------------------------------------
# check injection mode

if(inj=='low'){
    index_mat <- -abs(index_mat)
}
if(inj=='high'){
    index_mat <- abs(index_mat)
}


# --------------------------------------------
# inject outliers
 
# estimate from data if not simulated
sf <- DESeq2::estimateSizeFactorsForMatrix(k)
normtable <- t(t(k)/sf)
datasd <- matrix(rowSds(log2(normtable + 1)), nrow=nrow(ods), ncol=ncol(ods), byrow=FALSE)
lmu <- matrix(rowMeans(log2(normtable+1)), nrow=nrow(ods), ncol=ncol(ods), byrow=FALSE)

# use simulated values if present
if(all(c("true_mean", "true_sd") %in% assayNames(ods))){
    sf <- colData(ods)[['true_sizeFactor']]
    datasd <- assay(ods, "true_sd")
    lmu <- log2(assay(ods, "true_mean"))
}

# inject outliers
max_out <- 1E2 * min(max(k), .Machine$integer.max/1E3)
n_rejected <- 0
list_index <- which(index_mat != 0, arr.ind = TRUE)
for(i in seq_len(nrow(list_index))){
    row <- list_index[i,'row']
    col <- list_index[i,'col']
    fc <- zScore * datasd[row,col]
    clcount <- index_mat[row,col] * fc + lmu[row,col]
    
    #multiply size factor again
    art_out <- round(sf[col]*2^clcount)
    if(art_out < max_out){
        k[row,col] <- art_out
    }else{
        #remove super large outliers
        index_mat[row,col] <- 0 
        n_rejected <- n_rejected + 1
    }
}
mode(k) <- "integer"


# create outrider object
ods_out <- copy(ods)
assay(ods_out, 'trueCounts') <- counts(ods)
counts(ods_out) <- k
assay(ods, "inj_mask") <- index_mat

#save tables
#save rawtable which is now the injected one.
saveRDS(ods_out, out_inj_counts)


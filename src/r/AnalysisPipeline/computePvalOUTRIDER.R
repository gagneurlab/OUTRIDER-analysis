in_ods      <- "Output/data/RecallAnalysis/nsamples=All/inj=both/zscore=3/correction=baseCooks/q=best/GTEx_not_sun_exposed_OUTRIDERfit.RDS"
out_ods     <- "Output/data/RecallAnalysis/nsamples=All/inj=both/zscore=3/correction=baseCooks/q=best/alternative=two.sided_pmethod=BY/GTEx_not_sun_exposed_OUTRIDERpVal.RDS"
correction  <- "baseCooks"
alternative <- "two.sided"
pmethod     <- "BY"
threads     <- 10

# inputs:
in_ods <- snakemake@input[['outrider_fit']]

# outputs:
out_ods <- snakemake@output[['outrider_pval']]

# Params:
alternative <- snakemake@wildcards[['alternative']]
pmethod <- snakemake@wildcards[['pmethod']]
correction <- snakemake@wildcards[['correction']]

# Run params:
threads <- snakemake@threads

# Inferred params
# --------------------------------------------
source("./src/r/config.R")
register(MulticoreParam(threads))

ods <- readRDS(in_ods)
message(paste0(date(), ": P-value calculation ..."))
if(correction %in% c("baseCooks", "basePearsonRes")){
    assay(ods, 'pValue') <- matrix(1, ncol=ncol(ods), nrow=nrow(ods))
    assay(ods, 'padjust') <- assay(ods, 'pValue')
} else {
    ods <- computePvalues(ods, method=pmethod, alternative=alternative)
}

saveRDS(ods, out_ods)

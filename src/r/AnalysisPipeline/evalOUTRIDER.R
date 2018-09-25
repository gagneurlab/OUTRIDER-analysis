# debug files
if(FALSE){
    snakemake <- readRDS('tmp.snakemake.RDS'); snakemake 
    in_ods       <- "Output/data/RecallAnalysis/nsamples=All/inj=both/zscore=3/correction=baseCooks/q=best/alternative=two.sided_pmethod=BY/GTEx_not_sun_exposed_OUTRIDERpVal.RDS"
    out_plot_tsv <- "Output/data/RecallAnalysis/nsamples=All/inj=both/zscore=3/correction=baseCooks/q=best/alternative=two.sided_pmethod=BY/GTEx_not_sun_exposed_plot.tsv"
    
    Nsamples    <- "All"
    inj         <- "both"
    inj_value   <- "3"
    correction  <- "beseCooks"
    alternative <- "two.sided"
    pmethod     <- "BY"
    q           <- "best"
    dataset     <- "GTEx_not_sun_exposed"
}

source("./src/r/config.R")

# inputs:
in_ods <- snakemake@input[['outrider_pval']]

# outputs:
out_plot_tsv <- snakemake@output[['plot_table']]

# Params:
threads <- snakemake@threads

# 
# Wildcards
Nsamples    <- snakemake@wildcards[['nsamples']]
inj         <- snakemake@wildcards[['inj']]
inj_value   <- snakemake@wildcards[['inj_value']]
correction  <- snakemake@wildcards[['correction']]
pmethod     <- snakemake@wildcards[['pmethod']]
q           <- snakemake@wildcards[['q']]
dataset     <- snakemake@wildcards[['dataset']]

#------------------------

ods <- readRDS(in_ods)

dt_eval <- createEvalTable(ods, q=q, correction=correction,
        anMask='inj_mask', anTrueCts='trueCounts', Nsamples=Nsamples, 
        inj=inj, inj_value=inj_value, pmethod=pmethod, dataset=dataset)

##save plot data table
write.table(dt_eval, file = out_plot_tsv, row.names = FALSE, sep = '\t', quote = FALSE)

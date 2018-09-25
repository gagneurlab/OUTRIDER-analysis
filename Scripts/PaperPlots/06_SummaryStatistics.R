#'---
#' title: Summary Statistics for OUTRIDER paper
#' author: Christian Mertes
#' wb:
#'  input: 
#'  - ODS_GTEx:       '`sm config["DATADIR"] + "/fitOutrider/" + config["FIGURE_IMPLEMENTATION"] + "/Skin_Not_Sun_Exposed_Suprapubic_ODS.RDS"`'
#'  - ODS_GTExpca:    '`sm config["DATADIR"] + "/fitOutrider/pca/Skin_Not_Sun_Exposed_Suprapubic_ODS.RDS"`'
#'  - ODS_GTExpeer:   '`sm config["DATADIR"] + "/fitOutrider/peer/Skin_Not_Sun_Exposed_Suprapubic_ODS.RDS"`'
#'  - ODS_Kremer:     '`sm config["DATADIR"] + "/fitOutrider/" + config["FIGURE_IMPLEMENTATION"] + "/Kremer_ODS.RDS"`'
#'  - ODS_Kremerpca:  '`sm config["DATADIR"] + "/fitOutrider/pca/Kremer_ODS.RDS"`'
#'  - ODS_Kremerpeer: '`sm config["DATADIR"] + "/fitOutrider/peer/Kremer_ODS.RDS"`'
#'  - ODS_GTEx_mask:  '`sm config["DATADIR"] + "/RecallAnalysis/nsamples=All/Skin_Not_Sun_Exposed_Suprapubic_outlier_mask.RDS"`'
#'  - ODS_GTEx_bench: '`sm config["DATADIR"] + "/RecallAnalysis/nsamples=All/inj=low/zscore=6/correction=" + config["FIGURE_IMPLEMENTATION"] + "/q=best/pmethod=BY/Skin_Not_Sun_Exposed_Suprapubic_OUTRIDERpVal.RDS"`'
#'  - ODS_RAW_GTEx:   '`sm config["DATADIR"] + "/rawODS/Skin_Not_Sun_Exposed_Suprapubic_ODS.RDS"`'
#'  - ODS_RAW_Kremer: '`sm config["DATADIR"] + "/rawODS/Kremer_ODS.RDS"`'
#'  - ODS_TISSUE_GTEx_outrider: '`sm expand(config["DATADIR"] + "/fitOutrider/" + config["FIGURE_IMPLEMENTATION"] + "/{tissue}_ODS.RDS", tissue=config["GTEx_tissues"])`'
#'  - ODS_TISSUE_GTEx_pca:      '`sm expand(config["DATADIR"] + "/fitOutrider/pca/{tissue}_ODS.RDS", tissue=config["GTEx_tissues"])`'
#'  - ODS_TISSUE_GTEx_peer:     '`sm expand(config["DATADIR"] + "/fitOutrider/peer/{tissue}_ODS.RDS", tissue=config["GTEx_tissues"])`'
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---

if(FALSE){
    snakemake <- readRDS('tmp.snakemake.RDS'); snakemake 
}

source("src/r/config.R")

#
# functions
# 
getDimensionsODS <- function(ods){
    c('Number of genes'=nrow(ods),
            'Number of samples'=ncol(ods))
}

getNumberOfAberrantSamples <- function(ods){
    hitPerSample <- aberrant(ods, padj=FDR_LIMIT, zScore=Z_LIMIT, by="s")
    list('#Samples with at least one hit'=sum(hitPerSample > 0),
            'Quantiles'=quantile(hitPerSample, p=0.9),
            'Median'=median(hitPerSample))
}

getNumberOfOutlierSamples <- function(ods){
    geneLimit <- length(ods) * OUTLIER_RATIO
    hitPerSample <- aberrant(ods, padj=FDR_LIMIT, zScore=Z_LIMIT, by="s")
    goodSamples <- !(hitPerSample > geneLimit)
    list(data=goodSamples, values=round(c(
            'Number of good samples'=sum(goodSamples), 
            '% of good samples'=sum(goodSamples)/length(goodSamples)*100), 1))
}

#' 
#' # Cutoff used
c(FDR_LIMIT=FDR_LIMIT, OUTLIER_RATIO=OUTLIER_RATIO, ZSCORE_LIMIT=Z_LIMIT)

#' # Data before filtering
#' GTEx data size after filtering
raw_ods_gtex <- readRDS(snakemake@input$ODS_RAW_GTEx)
getDimensionsODS(raw_ods_gtex)

#' Kremer data size before filtering
raw_ods_kremer <- readRDS(snakemake@input$ODS_RAW_Kremer)
getDimensionsODS(raw_ods_kremer)

#' # Best encoding dimension for each data set
#' GTEx skin best encoding
odsls_gtex <- list(
    OUTRIDER=readRDS(snakemake@input$ODS_GTEx),
    PCA=readRDS(snakemake@input$ODS_GTExpca),
    PEER=readRDS(snakemake@input$ODS_GTExpeer))
ods_gtex <- odsls_gtex[[1]]
getDimensionsODS(ods_gtex)
getBestQ(ods_gtex)
plotEncDimSearch(ods_gtex)

#' Kremer best encoding
odsls_kremer <- list(
    OUTRIDER=readRDS(snakemake@input$ODS_Kremer),
    PCA=readRDS(snakemake@input$ODS_Kremerpca),
    PEER=readRDS(snakemake@input$ODS_Kremerpeer))
ods_kremer <- odsls_kremer[[1]]
getDimensionsODS(ods_kremer)
getBestQ(ods_kremer)
plotEncDimSearch(ods_kremer)


#' # Detection of aberrant events
#' GTEx number of samples with at least one hit
sapply(odsls_gtex, getNumberOfAberrantSamples)

#' Kremer number of samples with at least one hit
sapply(odsls_kremer, getNumberOfAberrantSamples)

#' # Detection of outlier samples
#' GTEx number of outlier samples
sapply(odsls_gtex, function(x) getNumberOfOutlierSamples(x)$values)

#' Kremer number of outlier samples
sapply(odsls_kremer, function(x) getNumberOfOutlierSamples(x)$values)


#' # Benchmark with GTEx
#' Data used for benchmark
gtex_good_samples <- getNumberOfOutlierSamples(ods_gtex)$data
c("Number of outlier samples"=sum(!gtex_good_samples), 
    "Number of good samples"=sum(gtex_good_samples),
    "Number of genes"=length(ods_gtex))

#' Number of outliers injected
outlierMask <- readRDS(snakemake@input$ODS_GTEx_mask)
sum(assay(outlierMask, "inj_mask") != 0)

#' Get recall for mean ~ 80
ods <- readRDS(snakemake@input$ODS_GTEx_bench)
tmpdt1 <- createEvalTable(ods, correction="OUTRIDER", q=getBestQ(ods))
tmpdt2 <- unique(rbind(
    tmpdt1[zScore < -2, .(recall=sum(trueOutlier != 0)/max(numInjected), precision=sum(trueOutlier != 0)/.N, Correction=correction, ztype="2")],
    tmpdt1[zScore < -3, .(recall=sum(trueOutlier != 0)/max(numInjected), precision=sum(trueOutlier != 0)/.N, Correction=correction, ztype="3")]))
tmpdt2
tmpdt <- readBootstrapData(tmpdt1)
plotRibbonBenchmark(tmpdt, linetype = c(2,1)) +
    geom_point(data=tmpdt2, inherit.aes = FALSE, 
            aes(x=recall, precision, color=Correction, shape=ztype), size=3)

#trueOutlier <- read.table(file.path(DATADIR, BENCHMARK_DIR, ODS_BENCH_MASK_FILE), sep = '\t')
#dt2plot_bin <- getPRPlotData(ods, trueOutlier, 12)
#dt2plot_bin[,LowerBond:=as.numeric(gsub('\\(|,.*', '', bin, perl=TRUE))]
#binOfInt <- dt2plot_bin[LowerBond < 100 & LowerBond > 50 & grepl('FDR', Method)]
#binOfInt[recall == max(recall)]
#plot(binOfInt$recall, binOfInt$precision, type='l', xlim=c(0,1), ylim=c(0,1), 
#        xlab='Recall', ylab='Precision')
#grid()

#' ## Benchmark against Kremer (ALDH18A1)
gene='ALDH18A1'
sample="MUC1404"
padj_val <- assay(ods_kremer[gene], 'padjust')[,sample]
change_val <- counts(ods_kremer[gene], normalize=TRUE)[,sample] / 
        mean(counts(ods_kremer[gene], normalize=TRUE))
signif(c("FDR_value"=padj_val, "change"=change_val), 2)


plotVolcano(ods_kremer, sample, padjCut=FDR_LIMIT, zScoreCut=Z_LIMIT)
plotExpressionRank(ods_kremer, geneID = gene, padjCut=FDR_LIMIT, zScoreCut=Z_LIMIT)

#'
#' # Aberrant sample statistics
#' 
#' ## All GTEx tissues
#' 
getNumberOfAberrantSamples <- function(file){
    o <- readRDS(file)
    ans <- sum(aberrant(o, by='s') > nrow(o) *  OUTLIER_RATIO)
    return(list(nas=ans, dim=dim(o)))
}

files <- snakemake@input$ODS_TISSUE_GTEx_outrider
nasout <- mclapply(files, getNumberOfAberrantSamples, mc.cores=15)
files <- snakemake@input$ODS_TISSUE_GTEx_pca
naspca <- mclapply(files, getNumberOfAberrantSamples, mc.cores=15)
files <- snakemake@input$ODS_TISSUE_GTEx_peer
naspeer <- mclapply(files, getNumberOfAberrantSamples, mc.cores=15)

hist(sapply(nasout, '[[', 'nas'), xlab='Aberrant per tissue')
sum(sapply(nasout, '[[', 'nas'))
sum(sapply(nasout, function(x) x[['dim']][2]))
hist(sapply(naspca, '[[', 'nas'), xlab='Aberrant per tissue')
sum(sapply(naspca, '[[', 'nas'))
sum(sapply(naspca, function(x) x[['dim']][2]))
hist(sapply(naspeer, '[[', 'nas'), xlab='Aberrant per tissue')
sum(sapply(naspeer, '[[', 'nas'))
sum(sapply(naspeer, function(x) x[['dim']][2]))

getNumberOfAberrantSamples(snakemake@input$ODS_Kremer)
getNumberOfAberrantSamples(snakemake@input$ODS_Kremerpca)
getNumberOfAberrantSamples(snakemake@input$ODS_Kremerpeer)


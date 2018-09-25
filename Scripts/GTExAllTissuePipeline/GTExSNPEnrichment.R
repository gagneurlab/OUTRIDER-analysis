#'---
#' title: GTEx Rare variant enrichtment
#' author: Felix Brechtmann, Christian Mertes
#' wb:
#'   threads: 5
#'   input: 
#'     - variantTable:  '`sm config["DATADIR"] + "/GTEx_variant_enrichment/rareModerateNHighVariantsTable.tsv"`'
#'     - odsFiles: '`sm expand(config["DATADIR"] + "/fitOutrider/{method}/{{dataset}}_ODS.RDS", method=list(set(PAPER_METHODS + config["correction"])))`'
#'   output:
#'     - rds: '`sm config["DATADIR"] + "/GTEx_variant_enrichment/{dataset}_featureset.RDS"`'
#'     - ggplots: '`sm config["DATADIR"] + "/GTEx_variant_enrichment/{dataset}_featureset_ggplots.RDS"`'
#'     - wBhtml: '`sm config["htmlOutputPath"] + "/GTEx_variant_enrichment/{dataset}_run.html"`'
#'   type: noindex
#' output:
#'   html_document:
#'     code_folding: show
#'     code_download: TRUE
#'---

# Debug data
if(FALSE){
    snakemake <- readRDS('tmp.snakemake.RDS'); snakemake 
    
    tissue <- "Skin_Not_Sun_Exposed_Suprapubic"
    datadir <- "/s/project/scared/paper/revision/run2107/data"
    odsFiles <- list(
        PCA      = file.path(datadir, '/fitOutrider/pca', paste0(tissue, "_ODS.RDS")),
        PEER     = file.path(datadir, '/fitOutrider/peer', paste0(tissue, "_ODS.RDS")),
        OUTRIDER = file.path(datadir, '/fitOutrider/ed_NRob_NCR_NTC_YLAS', paste0(tissue, "_ODS.RDS")))
    vcfFile <- "Output/data/GTEx_variant_enrichment/rareModerateNHighVariantsTable.tsv"
}

# source config
source("./src/r/config.R")

#' # Input variables
LiDir <- '/s/project/scared/GTExV6PRareVariationData/GTExV6PRareVariationData'
odsFiles <- snakemake@input$odsFiles
sinnames <- basename(dirname(odsFiles))
names(odsFiles) <- sinnames
tissue <- snakemake@wildcards$dataset
vcfFile <- snakemake@input$variantTable
outRDS <- snakemake@output$rds
ggplotFile <- snakemake@output$ggplots

MAF_LIMIT <- 0.05
register(MulticoreParam(length(odsFiles)))
ggplots <- list()
curAEVersion <- snakemake@config$FIGURE_IMPLEMENTATION

odsFiles
tissue
vcfFile

###########################################
#' 
#' # VCF parsing
#' 
###########################################

#' Read in vcf file
#+ read vcf
variants <- fread(vcfFile)

#' Process the MAF column
#' * MAF for standard annotated variants
#' * 1/(2*#samples) for not annotated variants
#' * max(freq(allels)) for multi allele locations
variants[,MAF:=pmax(AF, maf2number(GMAF), na.rm=TRUE)]

hist(variants[,.(MAF=max(MAF)), by=c('subjectID', 'variantID')][,log10(MAF)],
     breaks=30, main='MAF distribution over all variants of interest', 
     xlab='log10(MAF)')
abline(v=log10(MAF_LIMIT), col='red')

#' MAF filter
variants <- variants[MAF<MAF_LIMIT]


#' Number of variants per sample
hist(variants[,.N,by=c('variantID', 'subjectID')][,.N,by=subjectID][,N],
        main='Number of variants per sample', xlab='Number of variants')

#' 
#' * Simplify Consequence annotation
#+ simplify consequences
variants <- simplifyConsequences(variants) 

#' Keep only the gene level and most severe variant annotation 
#' (remove transcript level)
variants <- variants[order(subjectID, Gene, rank)]
dupVars <- duplicated(variants[,.(subjectID, Gene)])
table(dupVars)
variants <- variants[!dupVars]
sort(table(variants[,.N,by=c('variantID', 'simple_conseq')][,simple_conseq], useNA='always'))

#' 
#' # Final var overview
#' 

#' * Final number of variants per sample
hist(variants[,.N,by=c('variantID', 'subjectID')][,.N,by=subjectID][,N])

#' Number of variants sample combinations
nrow(variants[,.N,by=c('variantID', 'subjectID')])

#' * Number of affected samples:
length(unique(variants$subjectID))

#' * Number of affected genes:
length(unique(variants$Gene))

###########################################
#' 
#' # Read in outlier calls
#' 
###########################################

#' 
#' ## Read the OUTRIDER data
#' 
#' Read in ods file
#+ read ods results
odsls <- readODSResults4Enrich(odsFiles, threads=length(odsFiles))
odsls

#' 
#' ## Load data from Li Rare variation paper
#+ read li et al results
LiData <- loadLiValues(LiDir, tissue)
dim(LiData$LiZscore)


#'
#' ## Filter as in Li for samples with less than 50 
#' 
sampleIDs <- LiData$LiMultTissOutliers[N<50, sampleID]
geneIDs <- gsub('\\..*$','', LiData$LiGeneType[type %in% c('protein_coding', 'lincRNA'), geneID])
features <- CJ('subjectID'=sampleIDs, 'geneID'=geneIDs)

#' 
#' ## Remove samples which are sequenced more than once
#' 
tmpdt <- odsls[[1]][,.(subjectID, geneID)]
dupSubjects <- unique(tmpdt[duplicated(tmpdt),subjectID])
rm(tmpdt)
if(length(dupSubjects) > 0){
    features <- features[!subjectID %in% dupSubjects]
}


#'
#' ## Merge data sets into one table
#+ merge data
var2merge <- unique(variants[, .(subjectID, geneID=Gene, simple_conseq, back_simple_conseq=simple_conseq, MAF, IMPACT)])
features2merge <- c(list(features=features, vars=var2merge, Li=LiData$LiZscore), odsls)
featuresSubset <- Reduce(x=features2merge, 
        function(x, y) merge(x, y, by=c('subjectID', 'geneID'), all.x=TRUE))
methodZ <- grep('_z$', names(featuresSubset), value=TRUE)
rm(var2merge)
rm(features2merge)
gc()

#' * Remove genes not tested by Li et al and by OUTRIDER
table(featuresSubset[,!is.na(simple_conseq)])
dim(featuresSubset)
featuresSubset <- featuresSubset[!is.na(LiZscore_z)]
dim(featuresSubset)
featuresSubset <- featuresSubset[!rowAlls(is.na(featuresSubset[,8:ncol(featuresSubset)]))]
dim(featuresSubset)

#' 
#' ## Number of variants we test for enrichment
#' 
ggplots[['NumVariants']] <- ggplot(featuresSubset[,.(Consequence=simple_conseq)], aes(Consequence)) + geom_bar() + 
    scale_y_log10() + theme(axis.text.x=element_text(angle=45, hjust=1))
ggplots[['NumVariants']]
table(featuresSubset[, !is.na(simple_conseq)])


#' 
#' # Overview of enrichment
#+ Zscore overview per variant type
#+ zscore distribution, fig.width=14
dt <- melt(featuresSubset, 
        measure.vars=grep('_z$', colnames(featuresSubset), value=TRUE),
        variable.name='Method', value.name='Zscore')
ggplots[['ZscorePerVarType']] <- ggplot(dt[abs(Zscore) > 2], aes(simple_conseq, Zscore, col=Method)) + 
    stat_summary(fun.data=function(x, ...){c(y=-15, label=length(x))}, 
                 position=position_dodge(width=.85), geom="text") +
    geom_boxplot(position=position_dodge(width=.85)) + 
    geom_violin( position=position_dodge(width=.85)) + 
    grids() + 
    theme(axis.text.x = element_text(angle=45, hjust=1))
ggplots[['ZscorePerVarType']]

#' 
#' #Add random samples
featuresSubset[,'RAND_z':=rnorm(.N, 0, 5)]

#'
#' # Zscore enrichment (all)
enrichdt <- rbindlist(lapply(methodZ, function(x){
    dt <- calculateEnrichment(featuresSubset, x, 2)
    dt$dt[,c('Method', 'enrichment'):=list(x, dt$enrichment)]
    dt$dt
}))
ggplots[['enrichAll']] <- ggplot(enrichdt[cutoff==TRUE], aes(Method, enrichment)) +
    geom_point() + 
    coord_flip() + 
    labs(title = 'Enrichment (All)')
ggplots[['enrichAll']]

#'
#' # Zscore enrichment (HIGH)
enrichdt <- rbindlist(lapply(methodZ, function(x){
    highsub <- copy(featuresSubset)
    highsub[IMPACT != 'HIGH', simple_conseq:=NA]
    dt <- calculateEnrichment(highsub, x, 2)
    dt$dt[,c('Method', 'enrichment'):=list(x, dt$enrichment)]
    dt$dt
}))
ggplots[['enrichHigh']] <- ggplot(enrichdt[cutoff==TRUE], aes(Method, enrichment)) +
    geom_point() + 
    coord_flip() +
    grids() + 
    labs(title = 'Enrichment (High)')
ggplots[['enrichHigh']]

#'
#' # Zscore enrichment (rare: MAF < 0.001)
enrichdt <- rbindlist(lapply(methodZ, function(x){
    highsub <- copy(featuresSubset)
    highsub[MAF < 0.005, simple_conseq:=NA]
    dt <- calculateEnrichment(highsub, x, 2)
    dt$dt[,c('Method', 'enrichment'):=list(x, dt$enrichment)]
    dt$dt
}))
ggplots[['enrichHigh']] <- ggplot(enrichdt[cutoff==TRUE], aes(Method, enrichment)) +
    geom_point() + 
    coord_flip() +
    grids() + 
    labs(title = 'Enrichment (MAF < 0.001)')
ggplots[['enrichHigh']]

#'
#' # P-value enrichment
calculateEnrichment(featuresSubset, paste0(curAEVersion, '_p'), 0.00005, FALSE)
calculateEnrichment(featuresSubset, paste0(curAEVersion, '_p'), 0.0005, FALSE)
calculateEnrichment(featuresSubset, paste0(curAEVersion, '_p'), 0.005, FALSE)
calculateEnrichment(featuresSubset, paste0(curAEVersion, '_p'), 0.05, FALSE)

calculateEnrichment(featuresSubset, 'peer_p', 0.00005, FALSE)
calculateEnrichment(featuresSubset, 'peer_p', 0.0005, FALSE)
calculateEnrichment(featuresSubset, 'peer_p', 0.005, FALSE)
calculateEnrichment(featuresSubset, 'peer_p', 0.05, FALSE)

#' 
#' * Z-score scatter plot for AE versus Li et al
#+ scatter z-score
ggplots[['scatterAE_Li']] <- ggplot(featuresSubset, aes(x=LiZscore_z, y=get(paste0(curAEVersion, '_z')))) + geom_hex(bins=50) + 
    geom_abline(intercept = 0, slope = 1, col='red') +
    scale_fill_gradient(trans = "log") + grids() + 
    labs(y=curAEVersion, title='Scatter OUTRIDER versus Li Zscore')
ggplots[['scatterPeer_Li']] <- ggplot(featuresSubset, aes(LiZscore_z, peer_z)) + geom_hex(bins=50) + 
    geom_abline(intercept = 0, slope = 1, col='red') +
    scale_fill_gradient(trans = "log") + grids() + 
    labs(y='PEER', title='Scatter PEER versus Li Zscore')
ggplots[['scatterAE_Li']]
ggplots[['scatterPeer_Li']]

#' 
#' # Recall Rank plots all (HIGH/MODERATE)
#' 
#+ create recall data 1
#' add random 
methods2plot <- colnames(featuresSubset)[7:ncol(featuresSubset)]
rrdt <- rbindlist(bplapply(methods2plot, dt=featuresSubset,
    function(x, dt){
        dt <- calculateRecallRank(dt, x, grepl('_p$', x))  
        dt <- data.table(Method=x, dt[,.(
            recall=get(paste0(x, "_recall")), 
            rank=get(paste0(x, '_rank')))])
        dt[,Type:=ifelse(grepl('_p$', Method), 'P-value', 'Z-score')]
        dt[,Method:=gsub('_[pz]$', '', Method)]
   }
))
rrdt

#+ recall plots 1
ggplots[['recallAll']] <- plotRecallRankForEnrichment(rrdt, maxRank=1e5, maxPoints=1e4, logx=F, logy=F) + labs(title=paste(tissue, '\nMODERATE + HIGH IMPACT rank < 100k'))
ggplots[['recallAll']]
plotRecallRankForEnrichment(rrdt, maxRank=1e6, maxPoints=1e4, logx=F, logy=F) + labs(title=paste(tissue, '\nMODERATE + HIGH IMPACT rank < 1mio'))
plotRecallRankForEnrichment(rrdt, maxRank=1e8, maxPoints=1e4, logx=T, logy=T) + labs(title=paste(tissue, '\nMODERATE + HIGH IMPACT all'))

#' 
#' # Recall Rank plots all (HIGH)
#' 
#+ create recall data 2
#' add random 
methods2plot <- colnames(featuresSubset)[7:ncol(featuresSubset)]
featuresSubset[IMPACT != 'HIGH', simple_conseq:=NA]
rrdt <- rbindlist(bplapply(methods2plot, dt=featuresSubset,
    function(x, dt){
        dt <- calculateRecallRank(dt, x, grepl('_p$', x))  
        dt <- data.table(Method=x, dt[,.(
            recall=get(paste0(x, "_recall")), 
            rank=get(paste0(x, '_rank')))])
        dt[,Type:=ifelse(grepl('_p$', Method), 'P-value', 'Z-score')]
        dt[,Method:=gsub('_[pz]$', '', Method)]
    }
))
rrdt

#+ recall plots 2
ggplots[['recallHigh']] <- plotRecallRankForEnrichment(rrdt, maxRank=1e5, maxPoints=1e4, logx=F, logy=F) + labs(title=paste(tissue, '\nHIGH IMPACT rank < 100k'))
ggplots[['recallHigh']]
plotRecallRankForEnrichment(rrdt, maxRank=1e6, maxPoints=1e4, logx=F, logy=F) + labs(title=paste(tissue, '\nHIGH IMPACT rank < 1mio'))
plotRecallRankForEnrichment(rrdt, maxRank=1e8, maxPoints=1e4, logx=F, logy=F) + labs(title=paste(tissue, '\nHIGH IMPACT all'))


#' 
#' # Recall Rank plots all (MAF < 0.01)
#' 
#+ create recall data 3
#' add random
methods2plot <- colnames(featuresSubset)[7:ncol(featuresSubset)]
featuresSubset[, simple_conseq:=back_simple_conseq]
featuresSubset[MAF > 0.01, simple_conseq:=NA]
rrdt <- rbindlist(bplapply(methods2plot, dt=featuresSubset,
    function(x, dt){
        dt <- calculateRecallRank(dt, x, grepl('_p$', x))  
        dt <- data.table(Method=x, dt[,.(
            recall=get(paste0(x, "_recall")), 
            rank=get(paste0(x, '_rank')))])
        dt[,Type:=ifelse(grepl('_p$', Method), 'P-value', 'Z-score')]
        dt[,Method:=gsub('_[pz]$', '', Method)]
    }
))
rrdt

#+ recall plots 3
ggplots[['recallRare']] <- plotRecallRankForEnrichment(rrdt, maxRank=1e5, maxPoints=1e4, logx=F, logy=F) + labs(title=paste(tissue, '\nMAF < 0.01 rank < 100k'))
ggplots[['recallRare']]
plotRecallRankForEnrichment(rrdt, maxRank=1e6, maxPoints=1e4, logx=F, logy=F) + labs(title=paste(tissue, '\nMAF < 0.01 rank < 1mio'))
plotRecallRankForEnrichment(rrdt, maxRank=1e8, maxPoints=1e4, logx=T, logy=T) + labs(title=paste(tissue, '\nMAF < 0.01 all'))


#+ Save results
saveRDS(featuresSubset, outRDS)
saveRDS(ggplots, ggplotFile)


#'---
#' title: Simulations for normal distributed counts
#' author: Felix Brechtmann, Christian Mertes
#' wb:
#'   output:
#'     - odsBestQ: '`sm config["DATADIR"] + "/bestQFitODS/SimulationNorm_fitN_Q{q}_ODS.RDS"`'
#'     - odsFit:   '`sm config["DATADIR"] + "/filteredODS/SimulationNorm_fitY_Q{q}_ODS.RDS"`'
#'     - wBhtml: 'Output/html/Simulations/sim_norm_q{q,\d+}.html'
#'   type: noindex
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#' threads: 1
#'---


q <- 5

source('src/r/config.R')
outFileNFit <- snakemake@output$odsBestQ
outFileYFit <- snakemake@output$odsFit
q <- as.numeric(snakemake@wildcards$q)
set.seed(500)

#' 
#' # Simulate data
#' 
n <- 80                                   # samples
p <- 5000                                 # genes
q <- q                                    # latent space dimension
sdSim <- pmax(0.01, rnorm(p, 0.7, 0.3))   # dispersion
logMean <- 5                              # log offset for mean expression level.
s <- rnorm(n,mean=1, sd = 0.1)            # size factors.
sdVec <- rep(0.5, p)                      # sd for H matrix

#'
#' # Simulate covariates.
#'
H_true  <- matrix(rnorm(n*q), nrow=n, ncol=q)
D_true <- matrix(rnorm(p*q, sd = sdVec), nrow=p, ncol=q)
plot(rowSds(D_true), log='y')
y_true  <- D_true %*% t(cbind(H_true))
mu      <- t(t(exp(rnorm(p, logMean) + y_true))*s)

#'
#' # Simulate count Matrix with specified means.
#'
k <- matrix(round(rlnorm(n*p, meanlog = log(mu), sdlog = sdSim)), nrow=p, ncol=n)
mode(k) <- 'integer'
hist(log10(k)+1, breaks=100)
hist(log10(1/sdSim), breaks=100)

#'
#' Create Outrider data set
#'
ods <- OutriderDataSet(countData=k)
assay(ods, "true_mean")           <- mu
assay(ods, "true_sd")             <- matrix(sdSim/log(2), nrow=p, ncol=n)
mcols(ods)[,"true_theta"]         <- 1/sdSim
colData(ods)[['true_sizeFactor']] <- s

#' # Data summary
print(paste0('Number of samples: ', dim(ods)[2]))
print(paste0('Number of genes: ', dim(ods)[1]))

metadata(ods)[['optimalEncDim']] <- q
metadata(ods)[['encDimTable']] <- data.table(
        encodingDimension=q, evaluationLoss=1, evalMethod='simulation')


#
# write out files
# 
saveRDS(ods, outFileNFit)
saveRDS(ods, outFileYFit)

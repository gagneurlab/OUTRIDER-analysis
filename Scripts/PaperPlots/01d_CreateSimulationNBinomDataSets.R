#'---
#' title: Simulation for NBinom data
#' author: Felix Brechtmann, Christian Mertes
#' wb:
#'   output:
#'     - odsBestQ: '`sm config["DATADIR"] + "/bestQFitODS/SimulationNBinom_fitN_Q{q}_ODS.RDS"`'
#'     - odsFit:   '`sm config["DATADIR"] + "/filteredODS/SimulationNBinom_fitY_Q{q}_ODS.RDS"`'
#'     - wBhtml: 'Output/html/Simulations/sim_nbinom_q{q,\d+}.html'
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
n <- 80                                        # samples
p <- 5000                                      # genes
q <- q                                         # latent space dimension
theta <- rlnorm(p, meanlog=log(180), sdlog=2)  # dispersion
logMean <- 5                                   # log offset for mean expression level.
s <- rnorm(n,mean=1, sd = 0.1)                 # size factors.
sdVec <- rep(0.5, p)                           # sd for H matrix

#'
#' # Simulate covariates.
#'
H_true  <- matrix(rnorm(n*q), nrow=n, ncol=q)
D_true <- matrix(rnorm(p*q, sd = sdVec), nrow=p, ncol=q)
plot(rowSds(D_true), log='y')
y_true  <- D_true %*% t(cbind(H_true))
mu      <- t(t(exp(rnorm(p, logMean) + y_true))*s)

sdLogScale <- function(mu, theta){
    sqrt(mu*(1+mu/theta))/(mu + 1)/log(2)
}
# scale it up to overcome the estimation error
true_sd <- sdLogScale(rowMeans(mu), theta) * 2.5

#'
#' # Simulate count Matrix with specified means.
#'
k <- matrix(rnbinom(n*p, mu=mu, size=theta), nrow=p, ncol=n)
mode(k) <- 'integer'
hist(log10(k)+1, breaks=100)
hist(log10(theta), breaks=100)

#'
#' Create Outrider data set
#'
ods <- OutriderDataSet(countData=k)
assay(ods, "true_mean")           <- mu
assay(ods, "true_sd")             <- matrix(true_sd, nrow=p, ncol=n)
mcols(ods)[,"true_theta"]         <- theta
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

getGTExDispEstimation <- function(){
    ods <- readRDS('/s/project/scared/paper/revision/run2107/data/GTEx_not_sun_exposed_maskDONE.RDS')
    hist(log10(mcols(ods)$disp), breaks=60, main='GTEx sun exposed disp')
    lines(log10(1:5000),  dnbinom(1:5000, mu=180, size=2)*350000, lwd=3, col='red')
}

createSimulationData <- function(samples=80, genes=5000, encDim=5, 
                    theta=c(mu=180, size=2, meanlog=25, sdlog=0.8), 
                    lognorm=TRUE, outlierFreq=1E-3, outlierZscore=5,
                    maxCount=1E8, ctsNBinom=TRUE){
    #' ## Simulation of data
    #' # samples
    n = samples 
    # genes
    p = genes
    # latent space dimension
    q = encDim 
    s = rnorm(n,mean=1, sd = 0.1)
    
    if(isScalarNumeric(theta)){
        theta <- rep(theta,p)
    } else {
        if(isTRUE(lognorm)){
            if(all(c('meanlog', 'sdlog') %in% names(theta))){
                theta <- theta[c('meanlog', 'sdlog')]
            }
            theta <- rlnorm(p, meanlog =log(theta[1]), sdlog = tetha[2])
        } else {
            theta <- pmax(0.5, rnbinom(p, mu=theta[1], size=theta[2]))
        }
    }
    
    freq = outlierFreq
    #zScore=10
    zScore= outlierZscore
    
    #' Aproximates standard deviation of counts in log2 space.
    #'
    #'@noRd
    sdLogScale <- function(mu, disp){
        sqrt(mu*(1+mu/disp))/(mu + 1)/log(2)
    }
    
    
    h_true <- matrix(rnorm(n*q), nrow=n, ncol=q)
    
    sdVec <- rep(1, p)
    #sdVec <- c(rep(0.01, 15E2), rep(0.5, 15E2), rep(0.1,20E2))
    Wd_true <- matrix(rnorm(p*q, sd = sdVec), nrow=p, ncol=q)
    plot(rowSds(Wd_true), log='y')
    y_true <- Wd_true%*%t(cbind(h_true))
    mu <- rnorm(p,5)
    mu <- t(t(exp(mu + y_true))*s)
    #theta <- 1/(0.1 + 0.1/rowMeans(mu))
    plot(rowMeans(mu), 1/theta, log = 'xy')
    
    k <- matrix(rnbinom(n*p, mu = mu, size = theta), nrow=p, ncol=n)
    k[k > maxCount] <- maxCount
    mode(k) <- 'integer'
    
    #mu <- mu[rowSums(k)>0 & !rowAnyNAs(k) & 
    #                       !rowAnys(k>0.1*.Machine$integer.max),]
    
    #k <- k[rowSums(k)>0 & !rowAnyNAs(k) & 
    #           !rowAnys(k>0.1*.Machine$integer.max),]
    
    
    dim(k)
    hist(log10(k), breaks=20)
    
    
    ## generate in-silico outliers.
    # generate index of injected counts
    index <- matrix(sample(c(0,1,-1), p*n, prob = c(1 - freq, freq/2, freq/2), 
                           replace = TRUE), nrow = p)
    # inject on low, high or both sides
    # if(inj=='low'){
    #     index <- -abs(index)
    # }
    # if(inj=='high'){
    #     index <- abs(index)
    # }
    list_index <- which(index != 0, arr.ind = TRUE)
    
    max_out <- 1E2 * min(max(k), .Machine$integer.max/1E3)
    n_rejected <- 0
    for(i in seq_len(nrow(list_index))){
        row <- list_index[i,'row']
        col <- list_index[i,'col']
        fc <- zScore * sdLogScale(mu[row,col], theta[row])
        clcount <- index[row,col]*fc + 
            log2(1 + mu[row,col])
        #multiply size factor again
        art_out <- round(s[col]*2^clcount)
        if(art_out < max_out){
            k[row,col] <- art_out
        }else{
            #remove super large outliers
            index[row,col] <- 0 
            n_rejected <- n_rejected + 1
        }
        
    }
    mode(k) <- "integer"
    
    return(list(k=k, mu=mu, s=s, theta=theta, index=index, sdVec=sdVec, 
            w=Wd_true, q=q))
}

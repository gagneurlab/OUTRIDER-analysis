renameDatasetNames <- function(x, column='x'){
    isDT <- is.data.table(x)
    if(isFALSE(isDT)){
        x <- data.table(col=x)
        setnames(x, 'col', column)
    }
    naming <- c(
        GTEx                    = 'Skin_Not_Sun_Exposed_Suprapubic', 
        Kremer                  = 'Kremer',
        'Simulation (NB, q=10)' = 'SimulationNBinom_fitN_Q10',
        'Simulation (LN, q=10)' = 'SimulationNorm_fitN_Q10')
    
    x <- x[,as.character(get(column))]
    for(i in seq_along(naming)){
        x <- gsub(paste0('^', naming[i], '$'), names(naming)[i], x)
    }
    
    if(isTRUE(isDT)){
        x[,c(column):=list(x)]
    }
    
    return(x)
}


renameCorrectionMethods <- function(dt, column){
    naming <- c(
        PCA               = 'pca', 
        PEER              = 'peer',
        "Cook's"          = 'baseCooks',
        Pearson           = 'basePearsonRes',
        OUTRIDER          = CONFIG_YAML$FIGURE_IMPLEMENTATION,
        'robust OUTRIDER' = CONFIG_YAML$FIGURE_ROB_IMPLEMENTATION)
    
    x <- dt[,as.character(get(column))]
    for(i in seq_along(naming)){
        x <- gsub(paste0('^', naming[i], '$'), names(naming)[i], x)
    }
    dt[,c(column):=list(x)]
    
    return(dt)
}

correctPrecisionRankPlotNames <- function(dt){
    # Injection name
    dt[inj == 'both', inj:='Both']
    dt[inj == 'high', inj:='High']
    dt[inj == 'low',  inj:='Low']
    
    # Injection value
    dt[,inj_value:=gsub('zscore=', 'Simulated |Z| = ', inj_value)]
    
    # Rank name
    dt[correction == 'Pearson', type:='|Pearson| ranked']
    dt[correction == "Cook's",   type:="Cook's ranked"]
    
    dt[type == 'Z score rank', type:='|Z| ranked']
    dt[type == 'Z score rank', type:='|Z| ranked']
    dt[type == 'FDR filtered', type:=paste('Ranked & FDR <', FDR_LIMIT)]
    
    return(dt)
}

simplifyBinRanges <- function(dt, column, signifVal=2){
    fac2simp <- dt[,get(column)]
    if(!is.factor(fac2simp)){
        fac2simp <- factor(as.character(fac2simp))
    }
    
    binBorders <- strsplit(gsub('\\(|\\]', '', levels(fac2simp)), ',')
    levels(fac2simp) <- unlist(lapply(binBorders, function(x){
            paste0('(', paste(signif(as.numeric(x), signifVal), collapse=','), ']') }))
    
    dt[,c(column):=list(fac2simp)]
    return(dt)
}

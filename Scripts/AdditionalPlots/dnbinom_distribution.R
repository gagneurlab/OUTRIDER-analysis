#'---
#' title: dnbinom distribution for different theta and mu
#' author: Christian Mertes
#' wb:
#'  input: 
#'  output:
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---

library(data.table)

#+ fig.width=12, fig.height=6
ks <- 1:1e3
par(mfrow=c(1,2), mar=c(4,4,3,0) + 0.1)
plot(NA, log='yx',ylim=c(1e-6, 0.5), xlim=range(ks), xlab='k', ylab='dnbinom')
dtconf <- data.table(
    size = rep(c(0.25, 0.75, 1, 1.25, 2), 3),
    mu   = rep(c(1, 10, 1000), each=5),
    col  = rep(c('darkorange', 'darkgreen', 'firebrick'), each=5),
    lty  = rep(1:5, 3)
)
dtconf
sapply(1:nrow(dtconf), function(i){
    lines(ks, dnbinom(ks, size=dtconf[i,size], mu=dtconf[i,mu]), 
            col=dtconf[i,col], lty=dtconf[i,lty])
})
grid()
par(mar=c(1,1,1,1) + 1.1)
plot.new()
legend('center', title = 'Legend', bty = 'n',
       dtconf[,sprintf('S: %2.2f; M: %4i', round(size, 3), mu)],
       col=dtconf[,col], lty=dtconf[,lty], ncol=2, xpd = TRUE, seg.len = 1,
       text.width = 0.3)
mtext('Theta distribution', at = -0.2)

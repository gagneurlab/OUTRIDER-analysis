#'---
#' title: Figure 1a Outlier detection scheme
#' author: Christian Mertes
#' wb:
#'  output:
#'  - pdf: '`sm config["FIGDIR"] + "/Figure1a_outlier_scheme.png"`'
#'  - rds: '`sm config["DATADIR"] + "/Figure1a_outlier_scheme.rds"`'
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---


if(FALSE){
    snakemake <- readRDS('tmp.snakemake.RDS'); snakemake
    fig1aFilename <- 'tmp.pdf'
    fig1aRdsFilename <- 'tmp.rds'
}

fig1aFilename <- snakemake@output$pdf
fig1aRdsFilename <- snakemake@output$rds

source('src/r/config.R')
library(data.table)
library(ggplot2)
library(cowplot)
library(ggbeeswarm)

colors <- c('#e41a1c', '#377eb8', '#4daf4a')

set.seed(42)
dt1 <- data.table(
    Expression=rep(c(60, 200), each=5) + rnorm(10, 10, 10),
    Condition=rep(c('A', 'B'), each=5),
    col=rep(colors[1:2], each=5)
)
dt2 <- data.table(
    Expression=c(42, rep(200, 99) + rnorm(99, 0, 15)),
    Population='A',
    col=colors[3]
)

pt <- ggplot() + 
    annotate("text", x=0, y=0, size=8, fontface='bold', label='Typical experiment setup') +
    theme_nothing()
pt

p1 <- ggplot(data=dt1, aes(x=Condition, y=Expression)) + 
    geom_quasirandom(size=3, aes(col=col)) + 
    scale_color_manual(values=brewer.pal(3, 'Dark2')[1:2]) + 
    theme(legend.position="none",axis.ticks.y=element_blank(), axis.text.y=element_blank()) + 
    labs(title='Differential\nExpression Analysis\n(eg: DESeq2/edgeR)')
p1

p2 <- ggplot(dt1) + 
    geom_quasirandom(data=dt2, size=3, aes(x=Population, y=Expression, col=col)) + 
    theme(legend.position="none", axis.ticks.x=element_blank(), 
            axis.ticks.y=element_blank(), axis.text.y=element_blank()) + 
    scale_color_manual(values=brewer.pal(3, 'Dark2')[2]) + 
    scale_x_discrete(labels=' ') + 
    labs(title='Outlier\nDetection\n(OUTRIDER)')
p2

#' 
#' # Panel A Figure 1
#' 
#' fig.width=5, fig.height=6
mainp <- plot_grid(p1, p2)
fig1a <- plot_grid(pt, mainp, ncol=1, rel_heights=c(1,5))
fig1a

#' 
#' Save data
#' 
ggsave(fig1a, filename = fig1aFilename, width=4.5, height=5.3, device=cairo_pdf)
saveRDS(fig1a, fig1aRdsFilename)

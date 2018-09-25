# R function
# 
# Author: mertes, baderd
###############################################################################


plot_png_as_background = function(png_file, ...){
    require(png)
    #Get the plot information so the image will fill the plot box, and draw it
    limits_plot_dev <- par()$usr
    ima = readPNG(png_file)
    rasterImage(ima, limits_plot_dev[1], limits_plot_dev[3], limits_plot_dev[2], limits_plot_dev[4], ...)
}

plotAsImage <- function(expr, width=1200, height=1200, mar=c(1,1.5,1,1), ...){
    tmpF <- R.utils::tmpfile()
    png(tmpF, width=width, height=height)
        eval(expr)
    dev.off()
    par(mar=mar, ...)
    plot.new()
    plot_png_as_background(tmpF)
}


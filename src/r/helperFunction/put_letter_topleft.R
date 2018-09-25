#!/usr/bin/env Rscript
# TODO: Add comment
# 
# Author: baderda, mertes
###############################################################################

put_letter_topleft <- function(the_letter, y_shift=1.5, letter_width_factor=1.0,
                    cex=1.5, font=2, x_adjust=NULL, y_adjust=NULL, ...){
    x_pos = par('usr')[1] - diff(par('usr')[1:2]) / diff(par('plt')[1:2]) * 
            par('plt')[1] + strwidth(the_letter)*1.5
    x_pos = x_pos  + strwidth(the_letter)*letter_width_factor
    y_line= floor(par('mar')[3]) - y_shift
    if(!is.null(x_adjust)){
        x_pos <- x_pos - x_adjust
    }
    if(!is.null(y_adjust)){
        y_line <- y_line - y_adjust
    }
    mtext(text=the_letter, side=3, line=y_line, font=font, cex=cex, 
            at=x_pos, ...)
}

#plot(1)
#put_letter_topleft('A')

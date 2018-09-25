#
# for figures
#
suppressPackageStartupMessages(library(png))
DINA4_WIDTH  <-  8.27 * 0.95
DINA4_LENGTH <- 11.70 * 0.95

paper_pdf= function(filename, width=DINA4_WIDTH, ...){
    filename <- paste0(filename,'.pdf')
    message('Plotting PDF to ', filename)
    cairo_pdf(filename=filename, width=width, ..., family='Arial')
}

paper_png= function(filename, ...){
    filename <- paste0(filename,'.png')
    message('Plotting PNG to ', filename)
    png(filename=filename, ..., family='Arial', res=300)
}

#'
#' make complete paper figure given as pdf and png
#' also create each panel in a subfolder 
#' 
make_paper_figure <- function(plot_fcts, layout_matrix, window_ratio, 
                    penal_cex=1, figure_name, figdir=FIGDIR,
                    upperLetter=FALSE, ...){
	
	# order the functions by name so they are plotted the right way
	plot_fcts <- plot_fcts[sort(names(plot_fcts))]
	
	# get defined names
	figure_prefix = file.path(figdir, figure_name)
	panel_prefix = file.path(
	        paste0(figure_prefix, '_panel'), 
	        paste0(figure_name, '_panel_')
	)
	
	# check if the folder exists
	if(!file.exists(dirname(panel_prefix)))
		dir.create(dirname(panel_prefix), recursive = T)
	
	# get correct ratio and inches
	if(window_ratio < 1){
		height <- 12 
		width  <- height * window_ratio
	} else {
		width  <- 12 
		height <- width * 1 / (window_ratio)
	}
	
	# png part
	message(paste("Plot main figure png to dir:", figure_prefix))
	paper_png(figure_prefix, units = 'in', height=height, width=width)
	do_the_plotting(plot_fcts, layout_matrix, penal_cex, 
            upperLetter=upperLetter, ...)
	
	# pdf part
	message(paste("Plot main figure pdf to dir:", figure_prefix))
	paper_pdf(figure_prefix, height=height, width=width)
	do_the_plotting(plot_fcts, layout_matrix, penal_cex, 
            upperLetter=upperLetter, ...)
	
	
	# plot single panels
	message("Start plotting each panel ...")
	for(i in 1:length(plot_fcts)){
		fct = plot_fcts[[i]]
		letter = names(plot_fcts)[i]
		panel_name = paste0(panel_prefix, letter)
		message(paste("Plot panel ", letter, "with prefix", panel_prefix))
		
		# png
		paper_png(panel_name, units = 'in', height=height / 2, width=width / 2)
		fct()
		graphics.off()
		
		# pdf
		paper_pdf(panel_name, height=height / 2, width= width / 2)
		fct()
		graphics.off()
	}
}


#'
#' plot each panel given by the function
#' 
do_the_plotting <- function(plot_fcts, layout_matrix, penal_cex, 
            upperLetter=FALSE, ...){
	
	layout(layout_matrix)
	for(letter in names(plot_fcts)){
		message("Plot panel ", letter)
		
		# plot the panel
		par(cex = penal_cex[letter])
		plot_fcts[[letter]]()
		
		# print label
		letter_to_put = tolower(letter) 
		if(isTRUE(upperLetter)){
			letter_to_put = toupper(letter)
		}
		put_letter_topleft(letter_to_put, ...)
		
	}
	graphics.off()
}


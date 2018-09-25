
#' 
#' Continues scale with the form of 
#' a x 10^b aka 6 x 10^4
#' 
fancy_scientific <- function(l) {
    l <- format(l, scientific = TRUE)
    l <- gsub("^(.*)e", "'\\1'e", l)
    l <- gsub("e[+]?", "%*%10^", l)
    l <- gsub("'0'%*%10^00", "0", l, fixed=TRUE)
    parse(text=l)
}

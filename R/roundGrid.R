#' Round Grid
#' 
#' Round values to nearest centered value, so as to put on a grid
#' 
#' @param x numeric vector
#' @param frac numeric indicating fraction for rounding;
#' @details For example, 0.5 means to put things on a 0.5 degree grid, or to round everythign to nearest 0.5, such that all values end in .25 or .75 (the center of intervals of 0.5)
#' @export
roundGrid <- function(x, frac=1){
	# if frac is 1, then place in a 1º grid
	# if frac is 0.5, then place in the 0.5º grid
	floor(round(x/frac, 6))*frac+frac/2
}
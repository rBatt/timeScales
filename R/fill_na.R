#' Fill-in NA's	
#' 
#' Fill NA valus with linearly interpolated values
#' 
#' @param x vector with NA's that is to be interpolated
#' 
#' @details
#' If starting or ended values are NA, repeats nearest non-NA value. If \code{x} is a time series (object of class \code{ts}), the output is also a \code{ts}.
#' 
#' @return a numeric vector or a \code{ts} object
#' @seealso \code{\link{approx}} \code{\link{ts}}
#' @export
fill_na <- function(x){
	xvec <- seq_along(x)
	navec <- is.na(x)
	filled_x <- stats::approx(xvec[!navec], y=x[!navec], xout=xvec, rule=2)$y
	if(stats::is.ts(x)){
		filled_x <- stats::ts(filled_x, freq=stats::frequency(x))
	}
	return(filled_x)
}
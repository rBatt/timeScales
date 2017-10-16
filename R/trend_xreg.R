#' Trend X-Reg
#' 
#' A matrix of (polynomial) trend regressors (X) over the length of a response variable (Y)
#' 
#' @param exp.order integer, the order of the polynomial. 1 is a linear trend, 2 a quadratic, and so on
#' @param y a vector that is presumably the response variable, only used to determine length/ size/ number of rows in regressor matrix
#' 
#' @details no intercept column is included, and the output matrix is scaled. The standardization (scaling) occurs *after* the squaring, such that the "trend1" (linear) column squared does *not* equal the "trend2" column (this should be obvious because the "trend2" column will have negative values due to the standardization post-squaring.)
#' 
#' @return a matrix with column names indicating the order of the polynomial
#' @seealso \code{\link{interaction_xreg}} \code{forecast::fourier}
#' @export
trend_xreg <- function(exp.order, y){
	x <- seq_along(y)
	hnames <- paste0("trend",1:exp.order)
	
	raise_order <- function(o){
		matrix(x^o, ncol=1)
	}
	omat <- sapply(1:exp.order, raise_order)
	dimnames(omat) <- list(NULL, hnames)
	return(scale(omat))
}
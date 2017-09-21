#' Embed and Subsample a Time Series
#' 
#' Embed a time series like in \code{\link{embed}}, and then subsample the rows of that resulting matrix.
#' 
#' @param x numeric vector
#' @param width integer size of window/ embedding dimension
#' @param n integer, sample every nth element of \code{x}
#' @param phase integer, start subsampling sequency on the nth element of \code{x}
#' 
#' @return matrix
sub_embed <- function(x, width=1, n=1, phase=1){
	emat <- embed(x, dimension=width)
	row_vec <- seq_len(nrow(emat))
	row_ind <- sub_samp(row_vec, n=n, phase=phase)
	return(emat[row_ind, ])
}
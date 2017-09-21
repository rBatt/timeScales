#' AC1
#' Calculate first-order autocorrelation
#' 
#' @param x numeric vector of values of a regular (equally-spaced) time series
#' @param ... unused, here for compatibility with functions that might, by default, pass arguments unused by this function
#' 
#' @return AR(1) coefficient
#' @export
ac1 <- function(x, ...){	
	# # option 1
	l2 <- stats::embed(x, 2)
	ac <- stats::cor(l2[,1], l2[,2], use="na.or.complete")
	#
	# # option 2
	# ac <- stats::ar(x, order.max=1)$ar # returns numeric(0) if nothing fit
	#
	# # option 3
	# l2 <- stats::embed(x, 2)
	# ac <- lm(I(l2[,1]) ~ I(l2[,2]))$coef[2]
	
	# option 4
	# detX <- timeScales::detrend(x)
	# if(all(is.na(x))){return(NA)}
#	# l2 <- stats::embed(x, 2)
#	# ac <- tryCatch(stats::lm(I(l2[,1]) ~ I(l2[,2]))$coef[2], error=function(cond)NA)
	
	# # option 4
	# l2 <- embed(x, 2)
	# out <- RollingWindow::RollingCorr(l2[,2], l2[,1], window=win)
	
	return(ac)
}


#' Subset AC1
#' 
#' Calculate first-order autocorrelation from a subset
#' 
#' @param x numeric vector, time series
#' @param n integer, subset to every n-th observation (every n-th pair of observations)
#' @param phase integer, initial observation (pair) to use when subsetting, passed to \code{\link{sub_samp}}
#' @param ... unused, here for compatability with other functions that might pass arguments by default
#' 
#' To do a subset for autocorrelation, the key is retain the adjacency of observations. For \code{n=2}, for example, we couldn't use every other observation in \code{x}, as this would fundamentally change the time scale of the calculation. Thus, we still want to calculate the correlation between adjacent samples. So in this example, we would use every other adjacent pair. E.g., samples 1-2, (skip 2-3), 3-4, (skip 4-5), ... and so on.
#' 
#' @return AR1 coefficient
#' @export
ac_sub <- function(x, n, phase, ...){
	l2 <- stats::embed(x, 2)
	row_vec <- seq_len(nrow(l2))
	row_ind <- sub_samp(row_vec, n=n, phase=phase)
	l2_sub <- l2[row_ind, ]
	ac <- stats::cor(l2_sub[,1], l2_sub[,2], use="na.or.complete")
	return(ac)
}
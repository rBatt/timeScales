#' AC1
#' Calculate first-order autocorrelation
#' 
#' @param x numeric vector of values of a regular (equally-spaced) time series
#' 
#' @return AR(1) coefficient
#' @export
ac1 <- function(x, ...){	
	# # option 1
	l2 <- embed(x, 2)
	ac <- cor(l2[,1], l2[,2], use="na.or.complete")
	#
	# # option 2
	# ac <- ar(x, order.max=1)$ar # returns numeric(0) if nothing fit
	#
	# # option 3
	# l2 <- embed(x, 2)
	# ac <- lm(I(l2[,1]) ~ I(l2[,2]))$coef[2]
	
	# option 4
	# detX <- detrend(x)
	# if(all(is.na(x))){return(NA)}
#	# l2 <- stats::embed(x, 2)
#	# ac <- tryCatch(stats::lm(I(l2[,1]) ~ I(l2[,2]))$coef[2], error=function(cond)NA)
	
	# # option 4
	# l2 <- embed(x, 2)
	# out <- RollingWindow::RollingCorr(l2[,2], l2[,1], window=win)
	
	return(ac)
}
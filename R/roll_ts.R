#' Roll TS
#' 
#' Apply a rolling window statistic to a time series
#' 
#' @param y numeric vector representing time series
#' @param width integer indicating number of samples per window
#' @param by integer indicating the number of samples by which to increment the window forward between windowed subsamples
#' @param FUN function to be applied to each window
#' @param x optional vector by which the observations in y are ordered (such as date-times)
#' @param DETREND logical, detrend or no? see \code{\link{detrendR}}
#' @param ... additional arguments to be passed to \code{FUN}
#' @param save_output logical, save detrending output to desktop?
#' 
#' The function performs a backward-looking rolling window on a numeric vector. Makes use of the \code{\link{embed}} function, which might be slow for long windows (many observations per window). Internally also makes use of \code{\link{sub_embed}}.
#' 
#' @return a vector with each value being the result of applying \code{FUN} to a window.
#' 
#' @export
#' 
#' @examples
#' # example data
#' x <- 1:50
#' y <- cumsum(rnorm(50))
#' dt <- data.table(x, y)
#' 
#' # do x and y values separately
#' # x value is an ordering value, like time
#' # y is the system state, like chlorophyll
#' xa <- roll_ts(dt[,x], FUN=max, width=5)
#' ya <- roll_ts(dt[,y], width=5)
#' 
#' # do both x and y at same time
#' dta <- dt[,agg_ts(y=y, x=x, width=5)]
#' 
#' # plot two sets of results
#' par(mfrow=c(2,1))
#' plot(xa, ya)
#' dta[,plot(x,y)]
roll_ts <- function(y, width=288, by=1, FUN=mean, x, DETREND=FALSE, save_output=FALSE, ...){
	buff <- rep(NA, width-1)
	if(by > 1){
		mat <- sub_embed(y, width=width, n=by) # sub_embed is for roll win, so subset 'n' is actually window 'by'
	}else{
		mat <- stats::embed(y, width)
	}
	
	# do detrending
	# write a function that will first detrend, then apply the function FUN
	# because detrending might require information about the frequency of the time series y, add that info to the subset z
	if(DETREND){
		stopifnot(stats::is.ts(y))
		FUN2 <- function(z, ...){ # z is the subset (window) from sub_embed() or embed()
			z <- stats::ts(z, freq=stats::frequency(y))
			mp <- 2 #4 # max order of the polynomial; 4 should cover most cases pretty easily
			mf <- floor(min(stats::frequency(y)/2, 2)) # default is fourier order of 6, but this might be too high if the frequency of the time series is high b/c the fourier order has to be limited to 1/2 of the frequency.
			z <- detrendR(z, max_poly=mp, max_fourier=mf, max_interaction=if(mf>0){2}else{0}, save_output=save_output) # detrending
			FUN(z, ...) # then apply FUN to the detrended series
		}
	}else{
		FUN2 <- FUN
	}
	
	agg <- c(buff, apply(X=mat, MARGIN=1, FUN=FUN2, ...))
	if(!missing(x)){
		if(by > 1){
			mat2 <- sub_embed(x, width=width, n=by)
		}else{
			mat2 <- stats::embed(x, width)
		}
		agg2 <- c(buff, apply(mat2, 1, max, na.rm=TRUE))
		data.table(x=agg2, y=agg)
	}else{
		agg
	}
}

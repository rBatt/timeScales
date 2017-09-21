#' Roll TS
#' 
#' Apply a rolling window statistic to a time series
#' 
#' @param y numeric vector representing time series
#' @param width integer indicating number of samples per window
#' @param by integer indicating the number of samples by which to increment the window forward between windowed subsamples
#' @param FUN function to be applied to each window
#' @param x optional vector by which the observations in y are ordered (such as date-times)
#' @param ... additional arguments to be passed to \code{FUN}
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
roll_ts <- function(y, width=288, by=1, FUN=mean, x, ...){
	buff <- rep(NA, width-1)
	if(by > 1){
		mat <- sub_embed(y, width=width, n=by) # sub_embed is for roll win, so subset 'n' is actually window 'by'
	}else{
		mat <- embed(y, width)
	}
	agg <- c(buff, apply(X=mat, MARGIN=1, FUN=FUN, ...))
	if(!missing(x)){
		if(by > 1){
			mat2 <- sub_embed(x, width=width, n=by)
		}else{
			mat2 <- embed(x, width)
		}
		agg2 <- c(buff, apply(mat2, 1, max, na.rm=TRUE))
		data.table(x=agg2, y=agg)
	}else{
		agg
	}
}

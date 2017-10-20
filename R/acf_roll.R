#' ACF Correlation
#' 
#' Helper function to calculate autocorrelation across many time series and then just return the correlations
#' 
#' @param x time series
#' @param ... paramters to pass to \code{stats::acf}
#' 
#' @return numeric vector of correlations
acf_cor <- function(x, ...){
	stats::acf(x, ..., na.action=stats::na.pass, plot=FALSE)$acf
}

#' Rolling Window ACF
#' 
#' Performs rolling window autocorrelation calculation across many time scales
#' 
#' @param x the time series
#' @param width width of rolling window
#' @param by increment window by this many time steps (n observations)
#' @param lag.max integer indicating the maximum lag at which to calculate the ACF (as in code{stats::acf})
#' @param DETREND whether or not to perform AICc-selected detrending; see \code{\link{detrendR}}
#' @param ... unused
#' 
#' @seealso \code{stats::ts} \code{forecast::fourier} \code{\link{trend_xreg}} \code{\link{interaction_xreg}} 
#' @export
acf_roll <- function(x, width=28*24, by=12, lag.max=width/14, DETREND=FALSE, ...){
	
	if(by > 1){
		mat <- sub_embed(x, width=width, n=by) # sub_embed is for roll win, so subset 'n' is actually window 'by'
	}else{
		mat <- stats::embed(x, width)
	}
	
	if(DETREND){
		detrend2 <- function(z){
			z <- stats::ts(z, frequency=stats::frequency(x))
			as.numeric(detrendR(z, max_poly=4, max_fourier=2, max_interaction=2))
		}
		mat <- t(apply(X=mat, MARGIN=1, FUN=detrend2))
	}
	
	out <- apply(X=mat, MARGIN=1, FUN=acf_cor, lag.max=lag.max)
	lag_lab <- 0:(nrow(out)-1)
	if(stats::is.ts(x)){
		obs_lab <- seq(from=stats::tsp(x)[1]+(width/stats::tsp(x)[3]), by=by/stats::tsp(x)[3], length.out=ncol(out))
	}else{
		obs_lab <- seq(from=width, by=by, length.out=ncol(out))
	}
	
	dimnames(out) <- list(lag=lag_lab, obs=obs_lab)
	out <- t(out)
	attr(out, "xlab") <- obs_lab
	attr(out, "ylab") <- lag_lab
	
	return(out)
	
}
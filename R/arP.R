#' AR(p)
#' 
#' Fit an AR(p) Model
#' 
#' @param x time series or numeric vector
#' @param method either 'ar' or 'auto'; 'ar' uses \code{ar} in the 'stats' package, 'auto' uses the \code{auto.arima} in the 'forecast' package.
#' @param stepwise logical, if method is "auto", indicates whether to perform stepwise model selection (fast) or search all models. As in \code{auto.arima}
#' 
#' @return ar coefficients
#' @export
arP <- function(x, method=c("ar", "auto"), stepwise=TRUE){
	method <- match.arg(method)
	if(method=="ar"){
		stats::ar(x, order.max=floor(length(x)/10))$ar
	}
	if(method=="auto"){
		mp <- max(floor(length(x)/20), 2)
		mod <- forecast::auto.arima(y=x, max.p=mp, max.order=mp*2, stepwise=stepwise, seasonal=FALSE, stationary=TRUE)
		co <- stats::coef(mod)
		ar_names <- names(co)[grepl("ar\\d", names(co))] # could just use logic to do subsetting, but leaving here so a bit easier to follow, debug, or extend
		return(co[ar_names])
	}
}
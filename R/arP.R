#' AR(p)
#' 
#' Fit an AR(p) Model
#' 
#' @param x time series or numeric vector
#' @param method either 'ar' or 'auto'; 'ar' uses \code{ar} in the 'stats' package, 'auto' uses the \code{auto.arima} in the 'forecast' package.
#' @param oType output type; default is to return leading eigenvalue; 'eigs' returns all eigenvalues, and 'coefs' returns the AR() coefficients
#' @param stepwise logical, if method is "auto", indicates whether to perform stepwise model selection (fast) or search all models. As in \code{auto.arima}
#' @param ... unused, but present to be compatible with \code{\link{ac1}} (so wrapper functions can accept either)
#' 
#' @details
#' If model selection finds no autocorrelation, then an AR(1) model is fit with \code{stats::arima}.
#' 
#' @return scalar if returning leading eigenvalue, or numeric vector of length p if return eigenvalues or ar() coefficients. For this last case the vector is named.
#' @export
arP <- function(x, method=c("ar", "auto"), oType=c('leadEig', 'eigs', 'coefs'), stepwise=TRUE, ...){
	oType <- match.arg(oType)
	method <- match.arg(method)
	if(method=="ar"){
		ar_coefs <- stats::ar(x, order.max=max(floor(length(x)/10),2), method="yw")$ar
	}
	if(method=="auto"){
		mp <- max(floor(length(x)/20), 2)
		mod <- forecast::auto.arima(y=x, max.p=mp, max.order=mp*2, stepwise=stepwise, seasonal=FALSE, stationary=TRUE)
		co <- stats::coef(mod)
		ar_names <- names(co)[grepl("ar\\d", names(co))] # could just use logic to do subsetting, but leaving here so a bit easier to follow, debug, or extend
		
		ar_coefs <- co[ar_names]
	}
	
	if(!length(ar_coefs)>=1){
		ar_coefs <- arima(x, order=c(1,0,0))$coef["ar1"]
	}
	
	if(oType=="coefs"){
		return(ar_coefs)
	}else{
		eigs <- arEigs(ar_coefs)
		if(oType=="eigs"){
			return(eigs)
		}else{
			return(max(Mod(eigs)))
		}
	}
	
}
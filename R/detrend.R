#' Detrend
#' Detrend time series
#' @param x vector time series
#' @param time vector of time steps
#' 
#' @return detrended x
#' @export
detrend <- function(x, time=1:length(x)){
	stats::residuals(stats::lm(x~time))
}
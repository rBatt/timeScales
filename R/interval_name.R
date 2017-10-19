#' Interval Name
#' 
#' Converts a time interval (integer number of observations) into words, assuming 5 minute data
#' 
#' @param x integer, the number of time steps
#' @param minPerSample time elapsed between samples, in minutes
#' 
#' @examples
#' interval_name(5)
#' interval_name(13)
#' interval_name(288*1.25)
#' 
#' @export
interval_name <- function(x, minPerSample=5){
	# assumes 5 minute data
	interv <- x*minPerSample
	if(interv < 60){
		unit <- "m"
		val <- interv
	}else if(interv >= 60 & interv < 1440){
		unit <- "h"
		val <- round(interv/60, 2)
	}else{
		unit <- "d"
		val <- round(interv/60/24, 2)
	}
	iname <- paste(val, unit, sep="-")
	return(iname)
}

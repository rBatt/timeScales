#' Interval Name
#' 
#' Converts a time interval (integer number of observations) into words, assuming 5 minute data
#' 
#' @param x integer, the number of time steps
#' 
#' @examples
#' interval_name(5)
#' interval_name(13)
#' interval_name(288*1.25)
#' 
#' @export
interval_name <- function(x){
	# assumes 5 minute data
	interv <- x*5
	if(interv < 60){
		unit <- "min"
		val <- interv
	}else if(interv >= 60 & interv < 1440){
		unit <- "hr"
		val <- round(interv/60, 2)
	}else{
		unit <- "day"
		val <- round(interv/60/24, 2)
	}
	iname <- paste(val, unit, sep="-")
	return(iname)
}

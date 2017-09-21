#' Agg TS
#' 
#' Aggregate Time Series
#' 
#' @param x time index vector
#' @param y value vector (numeric)
#' @param width integer indicating number of observations to aggregate at a time
#' @param na.rm logical, remove NA's? passed to \code{FUN}
#' @param FUN aggregation function
#' 
#' @return data.table with 3 columns: x, y, and N. N is number of observation per aggregated value (number of observations for input to \code{FUN})
#' @export
agg_ts <- function(x, y, width=288, na.rm=TRUE, FUN=mean){
	buff <- rep(NA, width-1)
	tot <- round(mean(1/diff(x), na.rm=TRUE), 0)
	frac <- width/tot
	dt <- data.table(x=roundGrid(x, frac), y=y)
	dto <- dt[,list(y=FUN(y, na.rm=na.rm), N=sum(!is.na(y))), by=x]
}
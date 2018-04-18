#' Rolling Autocorrelation, SOS
#' 
#' Calculate autocorrelation in a backwards-looking rolling window; intended for data set collect as part of the SoS Cascade project.
#' 
#' @param X a data.table of lake observations, or possibly list of data.tables. Should have columns named "variable", "lake", "y", and "x". See Details.
#' @param window_elapsed the number of observations that should have elapsed between the first and last observation; is the same as a "window size" when \code{n=1}
#' @param fit_arP logical, if TRUE (default is FALSE), perform model selection to determine the order of the ar() process, then calculate and return the ||eigenvalue|| instead of an AR(1) coeff
#' @param by integer indicating the number of time steps to increment the window as it "rolls". Can be a list or a vector if \code{X} is a list. Can affect rolling window or statistic subsample, see Details.
#' @param n integer indicating to subsample to every n-th observation when performing rolling window. Can be a list or a vector if \code{X} is a list. Can affect rolling window or statistic subsample, see Details.
#' @param phase starting datum to use when subsampling
#' @param vars a character, possibly vector, indicating which variables in "variable" column of \code{X} should be analyzed.
#' @param lakes a character, possibly vector, indicating which lakes in "lake" column of \code{X} should be analyzed.
#' @param subWindow logical, should the number of rolling windows by subset (uses the \code{by} and \code{n} arguments)
#' @param subStat logical, should the statistic be calculated for a subset of the windowed data?
#' @param DETREND logical, detrend or no? Pass to \code{\link{roll_ts}}. Also see \code{\link{detrendR}}.
#' @param save_output logical, save detrending output to desktop?
#' @param mp,mf,mi maximum order of polynomial, Fourier, and interaction
#' 
#' @details I need to fill this in. One thing to remember is that 28-day nominal window size is hard-coded into this function right now.
#' 
#' @return list or list of data.tables containing rolling window statistic
#' @import data.table
#' @export
roll_ac.sos <- function(X, window_elapsed, fit_arP=FALSE, by=1, n=1, phase=1, vars, lakes, subWindow=FALSE, subStat=FALSE, DETREND=FALSE, save_output=FALSE, mp=1, mf=floor(min(stats::frequency(y)/2, 2)), mi=if(mf>0){2}else{0}){
	if(any(unlist(n)!=1) & fit_arP){
		stop("some n>1 and arP is TRUE; Cannot do AR(p) models on downsampled data")
	}
	
	# check vars, set if missing
	if(missing(vars)){
		vars <- X[[1]][,unique(variable)]
	}else{
		vars <- match.arg(vars, choices=X[[1]][,unique(variable)], several.ok=TRUE)
	}
	
	# check lakes, set if missing
	if(missing(lakes)){
		lakes <- X[[1]][,unique(lake)]
	}else{
		lakes <- match.arg(lakes, choices=X[[1]][,unique(lake)], several.ok=TRUE)
	}
		
	# Create n and b based on window_elapsed if subWindow or subStat is TRUE
	# Otherwise, just ensure that by and n are lists
	# Note that lists are for mapply(), and mapply() recycles ... arguments to length of longest arg
	# So no need to check here for same length among X, window_elapsed, by, and n.
	w_e28 <- unlist(window_elapsed)/28
	if(subWindow){
		by <- as.list(w_e28)
	} else if(!is.list(by)){
		by <- as.list(by)
	}
	if(subStat){
		n <- as.list(w_e28)
	} else if(!is.list(n)){
		n <- as.list(n)
	}
	
	# check/ set window elapsed class
	if(!is.list(window_elapsed)){
		window_elapsed <- as.list(window_elapsed)
	}
	
	# The `n` argument is only used by the ac_sub function
	# because this function subsets the series while computing the statistic
	# (By contrast, the `by` argument indicates the size by which roll_ts should  
	# increment the rolling window positionand is unrelated to the statistic used)
	if(all(unlist(n)==1)){
		if(fit_arP){
			funUse <- arP # the defaults of this should do what I want
		}else{
			funUse <- ac1
		}
	}else{
		funUse <- ac_sub
	}
	
	# helper function to apply rolling statistic
	roll_ac <- function(X2, nsteps, by, n){
		X2[variable%in%vars & lake%in%lakes][,j={roll_ts(y=y, x=x, FUN=funUse, width=nsteps, by=by, n=n, phase=phase, DETREND=DETREND, mp=mp, mf=mf, mi=mi, save_output=save_output)}, by=c("lake","variable")]
	}
	out <- mapply(roll_ac, X, window_elapsed, by, n, SIMPLIFY=FALSE)
	out
}

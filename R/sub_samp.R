#' Sub Sample
#' 
#' Sub-sample a vector
#' 
#' @param x a vector to be subsampled
#' @param n an integer indicating to subsample every n-th datum
#' @param phase an integer indicating which sample should be the first to be sampled
#' 
#' @return a sub-sampled vector
#' @export
sub_samp <- function(x, n, phase=1){
	if(phase!=n){
		phase <- max(1, phase%%n) # ensures 1 ≤ phase ≤ n; keeps starting elements when phase > n, unlike min(phase,n)
	}
	vec <- rep(FALSE, n)
	vec[phase] <- TRUE
	x[vec]
}
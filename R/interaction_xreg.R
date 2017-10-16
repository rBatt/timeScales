#' Interaction (X) Regressors
#' 
#' Given two matrices of predictors p and f, get matrix of all columns of f interacting with each column of p
#' 
#' @param p matrix of polynomial trends, such as from \code{trend_xreg}
#' @param f matrix of fourier xregs, such as from \code{forecast::fourier}
#' @param i optional integer indicating the number of columns of p to be used in interaction; if NULL (default), interactions are with all column of p, and output is a list. If \code{i} is not NULL, output is a matrix.
#' 
#' @return a matrix with columns that are interactions between p and f, or a list where each element is all columns of f interacting with one of the columns of p, and each element of that list uses a different column of p
#' @seealso \code{\link{trend_xreg}} \code{forecast::fourier}
#' @export
#' @examples
#' freq <- 24
#' nDays <- 5
#' ampSlope <- 0.25
#' ts_time <- 1:(freq*nDays)
#' sineVec <- sqrt(ts_time)*ampSlope*sin(ts_time*(1/freq)*2*pi)
#' polyVec <- c(trend_xreg(2, ts_time)%*%c(0.75*3.5, -0.78*3))
#' determVec <- sineVec + polyVec
#' noiseVec <- rnorm(freq*nDays, sd=0.25)
#' tsVec <- ts(determVec+noiseVec, frequency=freq)
#' plot(tsVec, lwd=2)
#' 
#' test_p <- trend_xreg(2, tsVec)
#' test_f <- forecast::fourier(tsVec, K=2)
#' fit_mod0 <- lm(tsVec~test_p+test_f)
#' lines(ts(fitted(fit_mod0), freq=freq), col='blue')
#' 
#' test_i_NULL <- interaction_xreg(p=test_p, f=test_f, i=NULL)
#' str(test_i_NULL)
#' 
#' test_i_1 <- interaction_xreg(p=test_p, f=test_f, i=1)
#' str(test_i_1)
#' fit_mod <- lm(tsVec~test_p+test_f+test_i_1)
#' lines(ts(fitted(fit_mod), freq=freq), col='red')
interaction_xreg <- function(p, f, i=NULL){	
	if(is.null(i)){
		# if no i, interaction terms will include all polynomials
		maxPP <- ncol(p)
	}else{
		# if i specificied, don't calculate interactions with all polynomials, limit to i-th polynomial
		maxPP <- min(ncol(p), i)
	}
	
	# calculate interaction terms
	int_list <- structure(vector("list", maxPP), .Names=colnames(p)[1:maxPP])
	for(pp in 1:maxPP){
		PP <- p[,pp] #pick this trend
		pfmat <- matrix(NA, ncol=ncol(f), nrow=nrow(f))
		for(ff in 1:ncol(f)){	
			pfmat[,ff] <- f[,ff]*PP
		}
		colnames(pfmat) <- colnames(f)
		int_list[[pp]] <- pfmat
	}
	
	if(!is.null(i)){
		imm_list <- lapply(int_list[1:maxPP], function(x,f){x[,1:ncol(x)]})
		for(l in 1:length(imm_list)){
			dimnames(imm_list[[l]])[[2]] <- paste(dimnames(imm_list[[l]])[[2]], names(imm_list)[l], sep='-')
		}
		imm <- do.call(cbind, imm_list)
		return(imm)
		# mm <- cbind(intercept=1, pmm, fmm, imm)
	}else{
		return(int_list)
	}
}
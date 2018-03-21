#' Detrend
#' Detrend time series
#' @param x vector time series
#' @param time vector of time steps
#' 
#' @return detrended x
#' @export
detrend <- function(x, time=1:length(x)){
	stats::residuals(stats::lm(x~time, na.action=stats::na.exclude))
}


#' DetrendR 
#' 
#' Detrending of time series using polynomial trend, Fourier series, and their interaction; includes model selection
#' 
#' @param x time series to be detrended
#' @param max_poly integer indicating the order of the polynomial trend
#' @param max_fourier integer indicating the order of the Fourier series
#' @param max_interaction integer indicating the polynomial order that Fourier series should interact with (see \code{\link{interaction_xreg}})
#' @param returnType character indicating the type of output desired. "resid" returns the residuals of the series (the detrended time series), whereas "modelMatrix" returns the matrix of predictor variables. Both outputs are of the AICc-selected model.
#' @param save_output logical, save detrending output to desktop?
#' 
#' @details model to be used for detrending is fit via OLS, the best model is selected by AICc. See \code{\link{interaction_xreg}} for the construction of the interaction matrix. Not all model combinations are considered just all combinations of the orders of polynomial and Fourier series and interactions of the Fourier series with main-effect polynomials; the order of the interaction will never be greater than either of the orders of the polynomial or Fourier. Missing values (NA's) are first linearly interpolated.
#' 
#' @return detrended time series
#' @seealso \code{\link{detrend}}\code{stats::ts} \code{forecast::fourier} \code{\link{trend_xreg}} \code{\link{interaction_xreg}} \code{\link{fill_na}}
#' @export
detrendR <- function(x, max_poly=6, max_fourier=6, max_interaction=3, returnType=c("resid","modelMatrix"), save_output=FALSE){
	returnType <- match.arg(returnType)
	# check if interaction argument makes sense given poly and fourier orders
	if((max_poly==0 | max_fourier==0) & max_interaction!=0){
		warning("max_interaction > 0, but either max_poly or max_fourier == 0. Cannot have interaction w/o both polynomial trend and Fourier series. Proceeding with max_interaction set to 0.")
		max_interaction <- 0
	}
	
	stopifnot(max_poly >= 1)
	
	# polynomial matrix
	full_poly <- trend_xreg(max_poly, x)
	
	# fourier matrix
	if(max_fourier > 0){
		stopifnot(stats::is.ts(x))
	}
	full_fourier <- forecast::fourier(x, max_fourier)
	
	# interaction matrix
	if(max_interaction > 0){
		full_interaction_list <- interaction_xreg(p=full_poly, f=full_fourier)
	}
	
	# function to get model matrix
	get_mm <- function(p, f, i){
		# expects full_poly and full_fourier objects to exist in global/ higher level envir
		
		# polynomial
		if(p > 0){
			pmm <- full_poly[,1:p, drop=FALSE]
		}
		if(f ==0){ # if f ==0, then i should also be 0
			return(cbind(intercept=1, pmm))
		}
		
		# fourier
		get_fmm <- function(x, f){
			x[,1:min(2*f, ncol(x))]
		}
		fmm <- get_fmm(full_fourier, f)
		if(p == 0){
			return(cbind(intercept=1, fmm))
		}
		
		# subset full interaction matrix to correct order
		if(i == 0){
			mm <- cbind(intercept=1, pmm, fmm)
		}else{
			imm_list <- lapply(full_interaction_list[1:i], get_fmm, f=f)
			for(l in 1:length(imm_list)){
				dimnames(imm_list[[l]])[[2]] <- paste(dimnames(imm_list[[l]])[[2]], names(imm_list)[l], sep='-')
			}
			imm <- do.call(cbind, imm_list)
		
			mm <- cbind(intercept=1, pmm, fmm, imm)
		}
		return(mm)
	}
	
	# combinations of model parameters in terms of orders of poly, fourier, and interaction
	pfi_combos0 <- expand.grid(p=1:max_poly, f=0:max_fourier, i=0:max_interaction)
	limit_interactionOrder <- pfi_combos0[,3] <= pfi_combos0[,1] #& pfi_combos0[,3] <= pfi_combos0[,2]
	pfi_combos <- pfi_combos0[limit_interactionOrder,]
	
	# function to get model AICc given a row index (m) of the matrix (pfi_combos) of the possible combinations of model orders 
	mod_aicc <- function(m){
		X <- do.call(get_mm, as.list(pfi_combos[m,]))
		tX <- t(X)
		K <- ncol(X)
		bhat <- solve(tX%*%X)%*%tX%*%Y
		Yhat <- X%*%bhat
		resid <- Y - Yhat
		s2 <- stats::sd(resid)
		NLL <- -sum(stats::dnorm(resid, mean=0, sd=s2, log=TRUE))
		correction <- (2*K*(K+1))/(nrow(X) - K - 1)
		AICc <- 2*K + 2*NLL + correction
		return(AICc)
	}
	Y <- fill_na(as.numeric(x)) # use notation that's easier to recognize in the following
	AICcs <- sapply(1:nrow(pfi_combos), mod_aicc)
	
	# for best model (selected by AICc), re-fit model and get model matrix/ residuals
	X <- do.call(get_mm, as.list(pfi_combos[which.min(AICcs),]))
	if(save_output){
		if(file.exists("~/Desktop/detrend_output.csv")){
			write.csv(pfi_combos[which.min(AICcs),], file="~/Desktop/detrend_output.csv", append=TRUE)
		}else{
			write.csv(pfi_combos[which.min(AICcs),], file="~/Desktop/detrend_output.csv", append=FALSE)
		}
	}
	if(returnType=="modelMatrix"){
		return(X)
	}else if(returnType=="resid"){
		resid <- c(Y - X%*%(solve(t(X)%*%X)%*%t(X)%*%Y))
		resid <- stats::ts(resid, freq=stats::frequency(x))
		return(resid)
	}
	
}
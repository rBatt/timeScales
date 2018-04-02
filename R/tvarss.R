#' TVARSS
#' 
#' Fit a time-varying autoregressive state-space model in a Bayesian framework using JAGS
#' 
#' @param Y vector of observed values of time series
#' @param nP scalar, integer; order of the AR process
#' @param niter scalar, integer; number of MCMC iterations 
#' @param thinout scalar, the number of posterior samples you want; the thinning rate (n.thin) passed to JAGS is max(1, floor(niter/thinout))=n.thin.
#' @param tvMean logical; if TRUE, include a time-varying 'constant' parameter indicating that the mean of the time series changes over time
#' @param parallel logical; use parallel computing?
#' @param oType character indicating type of output; 'jags' is default, outputs an 'rjags' class object. If 'tvarss', output type is class 'tvarss', which is just a subset of the original 'rjags' object (just includes posterior samples of parameters as a list)
#' 
#' @details
#' The model is as follows.  
#'   
#' Observation:\cr
#' \eqn{y_t = z_t + v_t}
#' 
#' Process:\cr
#' \eqn{z_t = C_t + \Phi_t*B^p(Z_t - C_t) + w_t}
#' 
#' \eqn{B^p} is backshift of order p; if \eqn{p=2}\cr
#' \eqn{B^pZ = [z_{t-1}, z_{t-2}]^{-1}} ; note that the \code{[} are used here to indicate a matrix
#' 
#' Parameters are time-varying:\cr
#' \eqn{\Phi_t = \Phi_{t-1} + \epsilon_t}\cr
#' \eqn{C_t = C_{t-1} + \xi_t}\cr
#' 
#' @return a model of class 'rjags' or 'tvarss'
#' @export
tvarss <- function(Y, nP=1, niter=1E3, thinout=500, tvMean=FALSE, parallel=FALSE, oType=c("jags","tvarss")){
	requireNamespace("R2jags", quietly=TRUE)
	oType <- match.arg(oType)
	
	n <- length(Y)
	nP1 <- nP+1
	
	param_names <- c("Z", "Phi", "C", "tau_v", "tau_eps", "tau_w", "tau_xi")
	
	if(tvMean){
		# model_file <- "~/Documents/School&Work/epaPost/timeScales/inst/jags/tv_arp_mean.jags"
		model_file <- file.path(system.file(package="timeScales"), "jags/tv_arp_mean.jags")
	}else{
		# model_file <- "~/Documents/School&Work/epaPost/timeScales/inst/jags/tv_arp_noMean.jags"
		model_file <- file.path(system.file(package="timeScales"), "jags/tv_arp_noMean.jags")
	}
	
	inputData <- list(nP=nP, Y=Y, n=n, nP1=nP1)
	
	# initialValues <- rep(list(list(Z=rep(0.1, n), Phi=matrix(rep(0.1, n*nP),ncol=nP, nrow=n), C=rep(0.1,n), tau_v=1, tau_eps=1, tau_w=1, tau_xi=1)), 4)
	# out <- R2jags::jags.parallel(
	# 	data=names(inputData), inits=initialValues[1], parameters.to.save=param_names, model.file=model_file, n.chains=4, n.iter=3E3, export_obj_names=c("param_names","model_file", "initialValues")
	# )
	
	nthin <- max(1, floor(niter/thinout))
	if(parallel){
		out <- R2jags::jags.parallel(
			data=names(inputData), parameters.to.save=param_names, model.file=model_file, n.chains=4, n.iter=niter, n.thin=nthin, export_obj_names=c("param_names","model_file", "niter", "nthin"), envir = environment()
		)
	}else{
		out <- R2jags::jags(data=inputData, parameters.to.save=param_names, model.file=model_file, n.chains=4, n.iter=niter, n.thin=nthin)
	}
	
	if(oType=="tvarss"){
		out <- as.tvarss(out)
	}
	
	return(out)
}

#' Convert to Class 'tvarss'
#' 
#' Converts an object of class 'rjags' to class 'tvarss', which is just s subset of 'rjags'.
#' 
#' @param x an object of class 'rjags'
#' 
#' @return an object of class 'tvarss'
#' @export
as.tvarss <- function(x){
	if("tvarss"%in%class(x)){return(x)}
	stopifnot(class(x)%in%c("rjags", "rjags.parallel"))
	o <- x$BUGSoutput$sims.list
	class(o) <- "tvarss"
	return(o)
}

#' Summarize a 'tvarss' Posterior
#' 
#' Apply a summary statistic across the posterior samples of each parameter in a 'tvarss' object.
#' 
#' @param x a 'tvarss'  or 'rjags' object
#' @param FUN a function, or character string indicating a function, to be applied to posterior samples
#' 
#' @return a data.table with columns for each parameter, and rows for indices of that parameter. Scalar parameters will be repeated.
#' 
#' @note JAGS and \code{R2jags::jags()} accepts non-scalar parameter, such as vectors or matrices. In the context of time-varying autoregressive state space models, some parameters are going to vary over time. E.g., a parameter "Phi" might be an N-length vector containing parameters \code{phi[n=1], phi[n=2], ..., phi[n=N]}. As such, the output of this function will have 1 column for all "Phi", with each row representing a different phi[n]. The data.table will have N rows.
#' 
#' @export
summarize.tvarss <- function(x, FUN="mean"){
	requireNamespace("fields", quietly=TRUE)
	
	if(any(class(x)%in%c("rjags", "rjags.parallel"))){x <- as.tvarss(x)}
		
	FUN <- match.fun(FUN)
	
	ldim <- function(x){
		ds <- lapply(x, dim) # dimensions for each parameter
		ts_params <- sapply(ds, function(xx)xx[2] > 1)# parameters that are likely time series
		return(ts_params)
	}
	sw <- function(x, D, .FUN=FUN){
		apply(x, MARGIN=D, FUN=.FUN)
	}
	
	ts_params <- ldim(x)
	if(any(ts_params)){
		ts_sw <- lapply(x[ts_params], sw, D=c(2))
		nt <- length(ts_sw[[1]])
		results <- data.table(time=1:nt)
		results[,names(ts_params[ts_params]):=ts_sw]
	}else{
		nt <- 1
		results <- data.table(time=1:nt)
	}

	nonTS_sw <- sapply(x[!ts_params], FUN)
	results[,names(ts_params[!ts_params]):=as.list(nonTS_sw)]
	
	return(results)	
	
}

#' Plot the posterior of a 'tvarss' parameter
#' 
#' Plot a heat map representing the probability density of a 'tvarss' parameter that varies over time
#' 
#' @param x a 'tvarss' or 'rjags' object
#' @param varName NULL (default) or character indicating names of parameters in \code{x}. If NULL, plots all time-varying parameters
#' @param relative logical, if TRUE (default), at each time step the probability density of all values is divided by the maximum.
#' @param main,xlab,ylab character vector indicating the plotting labels to use. If a length-1 vector is supplied, the same character will be used for that label for all parameters
#' @param xvals optional numeric vector of values to be associated with the horizontal (x) axis. If NULL (default), x-axis values are assign a sequence of integers that increase by 1, and the length of this sequence is the number of parameters within a named parameter obejct. E.g., default behavior is to give the Phi parameters, of which there are N, the x-values of 1:N, such as for a time series of Phi that is of length N.
#' 
#' @return invisibly returns a named list of (possibly relative) probability densities
#' @export
plotPost.tvarss <- function(x, varName=NULL, relative=TRUE, main, xlab, ylab, xvals=NULL){
	if(any(class(x)%in%c("rjags", "rjags.parallel"))){x <- as.tvarss(x)}
	ldim <- function(x){
		ds <- lapply(x, dim) # dimensions for each parameter
		ts_params <- sapply(ds, function(xx)xx[2] > 1)# parameters that are likely time series
		return(ts_params)
	}
	# sw <- function(x, D, FUN=FUN){
	# 	apply(x, MARGIN=D, FUN=FUN)
	# }
	
	ts_params <- ldim(x)
	
	if(!is.null(varName)){
		stopifnot(all(varName%in%names(ts_params[ts_params])))
		ind <- ts_params & names(ts_params)%in%varName
		ts_params[!ind] <- FALSE
	}
	
	if(missing(main)){
		if(!is.null(varName)){
			main <- varName
		}else{
			main <- names(ts_params[ts_params])
		}
	}
	if(missing(xlab)){
		xlab <- rep("Time", sum(ts_params))
	}else{if(length(xlab)<sum(ts_params)&length(xlab)==1){
		xlab <- rep(xlab, sum(ts_params))
	}}
	if(missing(ylab)){
		ylab <- rep("Relative probability density", sum(ts_params))
	}else{if(length(ylab)<sum(ts_params)&length(ylab)==1){
		ylab <- rep(ylab, sum(ts_params))
	}}
	
	densL <- list()
	for(j in 1:sum(ts_params)){
		tx <- x[ts_params][[j]]
		ft <- range(c(tx)[-(1:nrow(tx))])
		if(is.null(xvals)){
			xvals <- 1:ncol(dens)
		}
		yvals <- do.call(seq, args=list(from=ft[1], to=ft[2], length.out=512))
		
		if(!relative){
			dens <- apply(tx, 2, function(xx){d<-density(xx, from=ft[1], to=ft[2])$y;return(d)})
		}else{
			# dens <- apply(tx, 2, function(xx){d<-density(xx, from=ft[1], to=ft[2])$y;d<-c(scale(d));return(d)})
			dens <- apply(tx, 2, function(xx){d<-density(xx, from=ft[1], to=ft[2])$y;d<-d/max(d);return(d)})
		}
		densL[[names(ts_params)[ts_params][j]]] <- dens
		image(xvals, yvals, t(dens), col=fields::tim.colors(256), main=main[j], ylab=ylab[j], xlab=xlab[j])		
	}
	invisible(densL)
}


#' TVARSS Wrapper Function
#' 
#' A wrapper function for tvarss() that handles 1) applying tvarss to a list of data.tables/ matrices containing time series to be analyzed; 2) detrending; 3) finding 'p' for AR(p) models; 4) calculating and extracting eigenvalues of the AR(p) models (and doing so in a way that handles them just like another 'parameter' from the model)
#' 
#' @param x an object conaining a list of data.tables to analyze. Each data.table should have a column named "variable", of which some elements need to be "X". When variable=="X", the a column 'y' of the data.table will be taken to be the time series to analyze. The column 'y' should be a "ts" object, with a specified \code{stats::frequency}
#' @param det logical, to detrend?
#' @param fit_arP logical, if TRUE, will fit an TVAR(p)SS model
#' @param niter number of iterations
#' @param thinout number of samples of the posterior to save in ouput (reduced from niter via 'thinning')
#' @param ... additional arguments to be passed to \code{\link{tvarss}}; note that the arguments \code{niter=3E3}, \code{thinout=500}, \code{parallel}, \code{nP}, and \code{oType="tvarss"} are already provided either as arguments to this wrapper function, or specified directly within this function.
#' 
#' @return a list of class "tvarss" objects; the names of the elements of this list will be the same as the names of \code{x}
#' @export
tvarss_wrapper <- function(x, det=FALSE, fit_arP=FALSE, niter=3E3, thinout=500, ...){
	
	ylist <- lapply(x, function(x){x[variable=="X", y]})
	if(det){
		mp <- 2
		mf <- floor(pmin(sapply(ylist, stats::frequency)/2, 2)) # set fourier order to 2 if sufficient samples per cycle, otherwise make 0
		mi <- c(0,mp)[(mf>0)+1]
		ylist <- mapply(detrendR, ylist, max_fourier=mf, max_interaction=mi, MoreArgs=list(max_poly=mp))
	}
	
	find_nP <- function(x){
		arp_mod_full <- ar(x, order.max=ceiling(stats::frequency(x))*5)
		nP <- length(arp_mod_full$ar)
		nP <- max(nP, 1)
		return(nP)
	}
	
	if(fit_arP){
		nPs <- sapply(ylist, find_nP)
	}else{
		nPs <- rep(1, length(ylist))
	}
	# tvarss_list <- lapply(ylist, tvarss, nP=nP, niter=2E3, parallel=TRUE, oType="tvarss", ...)
	# tvarss_list1 <- tvarss(ylist[[1]], nP=nPs[1], niter=2E3, parallel=TRUE, oType="tvarss") #mapply(tvarss, ylist, nP=nPs, MoreArgs=list(niter=2E3, parallel=TRUE, oType="tvarss"))
	tvarss_list <- mapply(tvarss, ylist, nP=nPs, MoreArgs=list(niter=niter, thinout=thinout, parallel=TRUE, oType="tvarss", ...), SIMPLIFY=FALSE)
	# tvarss_list <- mapply(tvarss, ylist, nP=nPs, MoreArgs=list(niter=niter, thinout=thinout, parallel=TRUE, oType="tvarss"), SIMPLIFY=FALSE) # for debugging
	# names(tvarss_list) <- names(x) # not needed, will keep names
	
	if(fit_arP){
		if(any(nPs>1)){
			requireNamespace("doParallel", quiety=TRUE)
			requireNamespace("foreach", quiety=TRUE)
			
			doParallel::registerDoParallel(cores=4)
			eigs <- foreach::foreach(j=1:length(tvarss_list), .combine=list, .multicombine=TRUE) %dopar% {
				apply(tvarss_list[[j]]$Phi, c(1,2), function(x)max(Mod(arEigs(x))))
			}
			for(j in 1:length(tvarss_list)){
				tvarss_list[[j]]$Eigen <- eigs[[j]]
			}
		}else{
			for(j in 1:length(tvarss_list)){
				tvarss_list[[j]]$Eigen <- tvarss_list[[j]]$Phi
			}
		}
	}
	return(tvarss_list)
}
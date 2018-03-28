#' TVARSS
#' 
#' Fit a time-varying autoregressive state-space model in a Bayesian framework using JAGS
#' 
#' @param Y vector of observed values of time series
#' @param nP scalar, integer; order of the AR process
#' @param niter scalar, integer; number of MCMC iterations 
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
tvarss <- function(Y, nP=1, niter=1E3, tvMean=FALSE, parallel=FALSE, oType=c("jags","tvarss")){
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
	
	if(parallel){
		out <- R2jags::jags.parallel(
			data=names(inputData), parameters.to.save=param_names, model.file=model_file, n.chains=4, n.iter=niter, export_obj_names=c("param_names","model_file")
		)
	}else{
		out <- R2jags::jags(data=inputData, parameters.to.save=param_names, model.file=model_file, n.chains=4, n.iter=niter)
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
	if(class(x)=="tvarss"){return(x)}
	stopifnot(class(x)%in%c("rjags"))
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
	
	if(class(x)=="rjags"){x <- as.tvarss(x)}
		
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
#' 
#' @return invisibly returns a named list of (possibly relative) probability densities
#' @export
plotPost.tvarss <- function(x, varName=NULL, relative=TRUE, main, xlab, ylab){
	if(class(x)=="rjags"){x <- as.tvarss(x)}
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
		xvals <- do.call(seq, args=list(from=ft[1], to=ft[2], length.out=512))
		if(!relative){
			dens <- apply(tx, 2, function(xx){d<-density(xx, from=ft[1], to=ft[2])$y;return(d)})
		}else{
			# dens <- apply(tx, 2, function(xx){d<-density(xx, from=ft[1], to=ft[2])$y;d<-c(scale(d));return(d)})
			dens <- apply(tx, 2, function(xx){d<-density(xx, from=ft[1], to=ft[2])$y;d<-d/max(d);return(d)})
		}
		densL[[names(ts_params)[ts_params][j]]] <- dens
		image(1:ncol(dens), xvals, t(dens), col=fields::tim.colors(256), main=main[j], ylab=ylab[j], xlab=xlab[j])		
	}
	invisible(densL)
}
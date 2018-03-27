library(data.table)
#phi[j-1] + 0.95*cos(2*pi*dt*j)*(dt*2*pi)
# dt <- 1/160 # period of the sine wave [so not really a 'dt', don't be confused by this]
 

y <- c(0) # observed state
z <- c(0) # true state
n <- 300 # number of time steps
sigma_w <- 0.5 # sd of the process error
sigma_v <- 0.01 # sd of the observation error
sigma_xi <- 0.0001 # sd of the mean parameter's random walk
sigma_eps <- 0.025 # sd of the AR(1) coefficient's random walk
phi0 <- 0.05 # starting value for the AR(1) coefficient
phi <- c(phi0) # will be a vector of AR(1) coefficients
C <- rep(0, n) # vector of the parameter governing the time series mean; will be filled in with non-0 later

# simulate a time-varying AR(1) model
# both the AR(1) coefficient and the mean parameter vary.
# based on Dakos & Ives 2012 Ecosphere
for(j in 2:n){
	phi[j] <- 1*phi[j-1] + rnorm(1, mean=0, sd=sigma_eps) #AR(1) coefficient
	C[j] <- C[j-1] + rnorm(1, mean=0, sd=sigma_xi) # mean parameter
	z[j] <- (z[j-1]-C[j])*phi[j] + C[j] + rnorm(1, mean=0, sd=sigma_w) # true state
	y[j] <- z[j] + rnorm(1, mean=0, sd=sigma_v) # observed state
}
# set vectors of interest to time series objects for convenient plotting
z <- ts(z)
y <- ts(y)
C <- ts(C)
phi <- ts(phi)


# ---- 'dlm' ----
# initial attempt at 'dlm' package
# # this package confused the hell out of me
# # steve said he had a hard time getting it to work too
# library(dlm)
# buildFun <- function(x){
# 	dlmModARMA(ar=x[1], sigma2=exp(x[2]))
# }
# dlmMLE(y, parm=c(0,0), build=buildFun)

# ---- 'tvReg' ----
# initial attempt at 'tvReg' package
# this package basically just does rolling window ac(), but instead of doing a regression of
# x_t vs x_{t-1}, it puts 'weights' on x_{t-1} (weighted least squares) where the weights are determined by the kernel
# library(tvReg) # https://cran.r-project.org/web/packages/tvReg/vignettes/tvReg-vignette.html
# tvmod <- tvAR(z, p=1, type='none', bw=0.1)
# plot(tvmod$tvcoef[,1], phi[-1])
# plot(tvmod)


# ---- Ryan + 'JAGS' ----
library(R2jags)
nP <- 1 # an AR(p) model
Y <- y # will supply JAGS with Y, in case I want to modify it from y (y simulated above)
# Y[sample(10:n, 0.1*(length(Y)))] <- NA # note that JAGS v4.0 and later can now handle missing observations in this type of model; i tested, works great
n <- length(Y) # n is the time series length
nP1 <- nP+1 # convenient
param_names <- c("Z", "Phi", "C", "tau_v", "tau_eps", "tau_w", "tau_xi") # parameters to track and save output for
model_file <- "~/Documents/School&Work/epaPost/timeScales/inst/jags/tv_arp_mean.jags" # where the tvarp.jags file is located on my (your) computer
inputData <- list(nP=nP, Y=Y, n=n, nP1=nP1) # input values; constants and observations, for JAGS

# run model
out <- R2jags::jags.parallel(
	data=names(inputData), parameters.to.save=param_names, model.file=model_file, n.chains=4, n.iter=5E3, export_obj_names=c("param_names","model_file")
)

# run model after setting initial values (seems to not matter)
# initialValues <- rep(list(list(Z=rep(0.1, n), Phi=matrix(rep(0.1, n*nP),ncol=nP, nrow=n), C=rep(0.1,n), tau_v=1, tau_eps=1, tau_w=1, tau_xi=1)), 4)
# out <- R2jags::jags.parallel(
# 	data=names(inputData), inits=initialValues[1], parameters.to.save=param_names, model.file=model_file, n.chains=4, n.iter=3E3, export_obj_names=c("param_names","model_file", "initialValues")
# )

out_mean <- as.data.table(out$BUGSoutput$mean[c("Z","Phi","C")]) # posterior means of time varying parameters
out_mean_consts <- out$BUGSoutput$mean[c("tau_eps","tau_v","tau_w","tau_xi")] # posterior means of static parameters
results <- data.table(y=y, Y=Y, z=z, C_true=C, phi_true=phi, out_mean) # I use the data.table package; install it if you don't have it; syntactically concise and computationally fast

# plot summary results as time series
par(mfrow=c(2,2))
# panel 1: observations, true state, and estimated state
results[,plot(y, ylim=range(y,z,Z), col='lightgray', lwd=3, type='l', ylab="State")]
results[,lines(z, col='black')]
results[,lines(Z, col='red')]

# panel 2: true AR(1) coefficient and estimated AR(1) coefficient
results[,plot(phi_true, ylim=range(c(phi_true, Phi.V1)), ylab="AR(1) Coefficient")]
results[,lines(Phi.V1, col='blue')]

# panel 3: true mean parmeter and estimated mean parameter
results[,plot(C_true, ylim=range(c(C_true, C)), ylab="Constant Parameter")]
results[,lines(C, col='blue')]

# print out a table comparing true and estimated static parameters
# note that these should only match where the simulated and true processes are identical
# e.g., if Phi is being simulated as a sine wave, don't compare sigma_eps
consts_table <- rbind(c(sigma_w, sigma_v, sigma_eps, sigma_xi),1/sqrt(c(out_mean_consts$tau_w, out_mean_consts$tau_v, out_mean_consts$tau_eps, out_mean_consts$tau_xi)))
colnames(consts_table) <- c("sigma_w","sigma_v","sigma_eps","sigma_xi")
rownames(consts_table) <- c("true", "estimated")
print(consts_table)

# assess estimates of true state 
# interesting to do when I replaced random values of the observation vector with 
results[is.na(Y), cor(z, Z)]
results[!is.na(Y), cor(z, Z)] 


# ===============================
# = Function to Estimate TVARSS =
# ===============================
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
#' The model is as follows.
#' Observation:
#' \eqn{y_t = z_t + v_t}
#' Process:
#' \eqn{z_t = C_t + \Phi_t*B^p(Z_t - C_t) + w_t}
#' \eqn{B^p} is backshift of order p; if \eqn{p=2}
#' \eqn{B^pZ = [z_{t-1}, z_{t-2}]^{-1}} ; note that the [ are used here to indicate a matrix
#' Parameters are time-varying:
#' \eqn{\Phi_t = \Phi_{t-1} + \epsilon_t}
#' \eqn{C_t = C_{t-1} + \xi_t}
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
		model_file <- "~/Documents/School&Work/epaPost/timeScales/inst/jags/tv_arp_mean.jags"
		# model_file <- file.path(system.file(package="timeScales"), "/inst/jags/tv_arp_mean.jags")
	}else{
		model_file <- "~/Documents/School&Work/epaPost/timeScales/inst/jags/tv_arp_noMean.jags"
		# model_file <- file.path(system.file(package="timeScales"), "/inst/jags/tv_arp_noMean.jags")
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
summary.tvarss <- function(x, FUN="mean"){
	if(class(x)=="rjags"){x <- as.tvarss(x)}
		
	FUN <- match.fun(FUN)
	
	ldim <- function(x){
		ds <- lapply(x, dim) # dimensions for each parameter
		ts_params <- sapply(ds, function(xx)xx[2] > 1)# parameters that are likely time series
		return(ts_params)
	}
	sw <- function(x, D){
		apply(x, MARGIN=D, FUN=FUN)
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
#' @param varName, NULL (default) or character indicating names of parameters in \code{x}. If NULL, plots all time-varying parameters
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
	sw <- function(x, D){
		apply(x, MARGIN=D, FUN=FUN)
	}
	
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




# ---- Fit Model w/o Mean ----
out_noMean <- tvarss(Y, niter=2E3)
summary(as.tvarss(out_noMean))[]
dev.new()
par(mfrow=c(2,2), mar=c(2.5,2.5,1,1), mgp=c(1.5,0.5,0), tcl=-0.2)
plot(phi, type='l')
plot(C, type='l')
plotPost.tvarss(test, varName=c('Phi'))

# ---- Fit Model w/ Time-Varying Mean ----
out_mean <- tvarss(Y, niter=2E3, tvMean=TRUE)
summary(as.tvarss(out_noMean))[]
dev.new()
par(mfrow=c(2,2), mar=c(2.5,2.5,1,1), mgp=c(1.5,0.5,0), tcl=-0.2)
plot(phi, type='l')
plot(C, type='l')
plotPost.tvarss(test, varName=c('Phi','C'))




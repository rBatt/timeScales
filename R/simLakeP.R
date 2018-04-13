# =====================
# = Packages Required =
# =====================
# You need to install.packages("rootSolve")


# ===================================
# = Draft of Function Documentation =
# ===================================
#' Eutrophication Bifurcation
#' Model of lake eutrophication via elevated phosphorus loading. When as phosphorus loading increases, a fold bifurcation is reached, and the system switches into a eutrophic state from a oligotrophic state. 
#' @param state state variables, a vector of length 3 with names "X", "M", and "U"
#' @param pars a named vector of parameters; currently only option is "C", a coefficient controlling the fraction of soil P that washes into the lake (units of 'per year'); set up this way so as to be compatible with functions in package "rootSolve" (C is a control parameter affecting water quality).
#' @param F input rate of phosporus to soil
#' @param C coefficient for transfer of soil phosphorus to the lake
# @param H noise for input to the lake
#' @param s sedimentation loss
#' @param h hydrologic loss (outflow)
#' @param r recycling coefficient
#' @param sigma standard deviation of recycling noise
#' @param b permanent burial rate of phosphorus in sediments
# @param dW a white noise process with mean zero and variance dt
#' @param m in the recycling function, the value of X at which recycling is half the maximum rate
#' @param q in the recycling function, the exponent q determines the slope of R(X) near m
#' @return change of state variables
#' @export

# =============================
# = Notes on Parameter Values =
# =============================
# # parameter table from Carpenter 2005 PNAS
# b=0.001 # per year
# C=0.00115 # per year
# F=31.6 # grams per meter squared per year
# h=0.15 # per year
# H=18.6 # grams per meter squared per year; different from carpenter and Brock 2006; in Carpenter 2005 PNAS this is "annual export of P from the watershed in farm products, per unit lake area". In Carpenter & Brock 2006 this is "noise for input to the lake"
# m=2.4 # grams per meter squared
# r=0.019 # per year
# q=8 # no units
# s=0.7 # per year
# W=0.147 # grams per meter squared per year, non-agriculture P input to watershed soil prior to disturbance
# WD=1.55 # grams per meter squared per year, non-agriculture P input to watershed soil after disturbance


# # parameters needed for carpenter & brock 2006 model
# b=0.001 # per year
# C=0.00115 # per year
# F= 14.6 # grams per meter squared per year # is 14.6 in the copy of Table S1 that Steve sent me; 31.6 is the value in the 2005 PNAS paper
# h=0.15 # per year
# m=2.4 # grams per meter squared
# r=0.019 # per year
# q=8 # no units
# s=0.7 # per year
# sigma=0.01 # scaling for dW_dt
# lambda=0.35

# =======================================
# = Model for 1 Deterministic Time Step =
# =======================================
modelDeterministic <- function(state, pars=c(C=0.00115), F=14.6, b=0.001, h=0.15, m=2.4, r=0.019, q=8, s=0.7, sigma=0.01){
	with(as.list(c(state, pars)),{ # allows referring to names as variables; convenient
		R <- (X^q)/(m^q + X^q) # recycling
		H <- 1 # set to 1 as per appendix
		MRX <- M*R
		rMRX <- r*MRX
		dU_dt <- F - C*U*H # soil P
		dX_dt <- C*U*H - (s+h)*X + rMRX # + MRX # water P
		dM_dt <- s*X - b*M - rMRX # - MRX # mud P
		c(dU_dt=dU_dt, dX_dt=dX_dt, dM_dt=dM_dt) # rate of change in UXM
	}) # end with
}

# ================================================================
# = Model Holding Soil P Constant (to solve for water and mud P) =
# ================================================================
#' Deterministic Eutrophication Model
#' A 2D Eutrophication model with state variables for water P (X), mud P (M), and P loading as bifurcation parameter (I)
#' @param state vector of length 2 with names X and M (water and mud P, respectively)
#' @param pars named vector of length 1 with name I (P loading)
#' @param b P burial coefficient
#' @param h outflow coefficient
#' @param m the value of water P when recycling from mud to water is half of its maximum value
#' @param r a scaling term for recycling between mud and water
#' @param q the steepness of the recycling function at \code{m}
#' @param s sedimentation rate (turns water P into mud P)
#' @return named vector (X, M) indicating rate of change for water and mud P
#' @references 
#' SR Carpenter and WA Brock 2006 Rising variance: a leading indicator of ecological transition. Ecology Letters 9: 311-318
#' @export
modelDeterministicXM <- function(state, pars=c(I=0.00115*750), b=0.001, h=0.15, m=2.4, r=0.019, q=8, s=0.7){
	with(as.list(c(state, pars)),{ # use with() so that I can refer to column names in state and names in pars as variables
		R <- (X^q)/(m^q + X^q) # recycling
		MRX <- M*R
		rMRX <- r*MRX
		dX_dt <- I - (s+h)*X + rMRX # water P dynamics
		dM_dt <- s*X - b*M - rMRX # sediment P dynamics
		c(dX_dt=dX_dt, dM_dt=dM_dt) # output rate of change for X and M
	}) # end with
}
#' Wrapper for model when finding Jacobian matrix
#' 
#' Wrapper for model when finding Jacobian matrix using \code{rootSolve} package
#' 
#' @param t time, for compatibility with rootSolve::jacobian.full, but not used otherwise
#' @param state state variables as in \code{\link{modelDeterministicXM}}
#' @param pars paramters as in \code{\link{modelDeterministicXM}}
#' @param ... additional arguments to pass to \code{\link{modelDeterministicXM}}
#' @return rate of change of state variables
#' @export
mDXM_jac <- function(t, state, pars, ...){
	modelDeterministicXM(state, pars, ...)
}

#' Wrapper for model when finding stability properties
#' 
#' Wrapper for model finding stability properties using \code{phaseR} package
#' 
#' @param t time, for compatibility with \code{phaseR} and \code{deSolve} and \code{rootSolve}
#' @param y state variables as in \code{\link{modelDeterministicXM}}
#' @param parameters same as \code{pars} in \code{\link{modelDeterministicXM}}; had to be called 'parameters' for compatibility with \code{phaseR}
#' @param ... other arguments to be passed to \code{\link{modelDeterministicXM}}
#' @return rate of change of state variables
#' @export
mDXM <- function(t, y, parameters, ...){
	names(y) <- c("X","M")
	list(modelDeterministicXM(y, pars=parameters, ...))
}

#' Get Initial Values for Eutrohpication Model
#' Gets a grid/ mesh of values for eventual input to model
#' @param gridN number of values for I and X
#' @param Irange range of values for I, input rate of P to lake
#' @param Xrange range of values for X, water phosphorus
#' @param Mrange range of values for M, mud P
#' @param In number of values for I; default is set to gridN
#' @param Xn number of values for X; default is set to gridN
#' @param Mn number of values for M; default is set to gridN
#' @return a data.frame of initial values
#' @export
getInit <- function(gridN=100, Irange=c(0.25, 1.75), Xrange=c(0.05, 8), Mrange=c(400, 400), In=gridN, Xn=gridN, Mn=gridN){
	# gridN <- 100
	Igrid <- seq(Irange[1], Irange[2], length.out=In)
	Xgrid <- seq(Xrange[1], Xrange[2], length.out=Xn)
	Mgrid <- seq(Mrange[1], Mrange[2], length.out=Mn)

	# sequences <- data.frame(I=Igrid, X=Xgrid, M=Mgrid)
	initialValues <- expand.grid(I=Igrid, X=Xgrid, M=Mgrid)
	# startMudP <-
	initialValues <- initialValues[order(initialValues[,"I"]),]
	rownames(initialValues) <- NULL
	initialValues
}

#' Get the Root of the Eutrophication Model
#' Finds the roots of the eutrophication model given starting values (water, mud P) and parameter (P loading)
#' @param x a vector of length 3 with names X (water P), M (mud P), and I (P loading)
#' @param pars vector of parameters to pass to P model
#' @param maxiter maximum number of iterations to use when finding root
#' @param ... additional arguments to pass to \code{\link{modelDeterministicXM}}
#' @return numeric vector with roots
#' @export
getRoot <- function(x, pars, maxiter=1E3, ...){
	if(missing(pars) | is.null(pars)){
		pars <- c(I=unname(x["I"]))
	}else{
		if(!"I"%in%names(pars)){
			pars <- c(I=unname(x["I"]), pars)
		}else{
			pars <- c(pars)
		}
	}
	tryCatch({
		rootSolve::multiroot(f=modelDeterministicXM, start=x[c("X","M")], parms=pars, maxiter=maxiter, ctol=1E-9, rtol=1E-9, atol=1E-9, ...)$root
	}, warning=function(w)c(X=NA,M=NA))
}

#' Get the Root of the Eutrophication Model from a Data Frame
#' Finds the roots of the eutrophication model given starting values (water, mud P) and parameter (P loading)
#' @param initialValues a \code{data.frame} with 3 columns named X (water P), M (mud P), and I (P loading)
#' @param maxiter maximum number of iterations
#' @param ... additional arguments to pass to \code{\link{modelDeterministicXM}}
#' @return a matrix with columns for initial values of X and M, roots (of X and M), and the bifurcation parameter I. Each row of output corresponds to a row of the input (though output reorded to ascending I). Original columns for state variables (X and M) will be renamed to "init.X" and "init.M". The equilibria that would be approached from these initial values, i.e. the 'roots', will take on the column names of "X" and "M". Thus, the columns X & M will not have the same values in the input as in the output (output are equilibrium values), unless the input states were already at equilibrium.
#' @export
getRoot_df <- function(initialValues, maxiter=1E3, ...){
	rootGrid <- t(apply(initialValues, 1, getRoot, maxiter=maxiter, ...))
	rootGrid_I <- data.matrix(data.frame(I=initialValues[,"I"], init=initialValues[,c("X","M")], rootGrid))
	rootGrid_I <- rootGrid_I[order(rootGrid_I[,"I"]),] # i should really remove this ...
	return(rootGrid_I)
}

#' Get Eigenvalues for Eutrophication Model
#' Given the rate of change of state variables (water and mud P) and the input P, returns the eigenvalues
#' @param x a named vector of length 3 (X, M, I), as in \code{\link{modelDeterministicXM}}
#' @param pars vector of parameters to pass to P model
#' @param ... additional arguments to pass to \code{\link{modelDeterministicXM}}
#' @return eigenvalues
#' @export
getEigs <- function(x, pars, ...){
	if(missing(pars) | is.null(pars)){
		pars <- c(I=unname(x["I"]))
	}else{
		if(!"I"%in%names(pars)){
			pars <- c(I=unname(x["I"]), pars)
		}else{
			pars <- c(pars)
		}
	}
	
	if(any(is.na(x))){return(c(NA,NA))}
	jac <- rootSolve::jacobian.full(y=x[c("X","M")], func=mDXM_jac, parms=x["I"], ...)
	eigen(jac)$values
}
#' Get Eigenvalues for Eutrophication Model (from data frame of states)
#' Given the rate of change of state variables (water and mud P) and the input P, returns the eigenvalues
#' @param rootGrid_I a matrix with columns named for X, M, I, as in \code{\link{modelDeterministicXM}}; these columns are the equilibria (roots) of the state variables, and I, the P input parameter.
#' @param ... additional arguments to pass to \code{\link{modelDeterministicXM}}
#' @note the "Grid" name just comes from legacy coding where I was trying a "grid" (as in, e.g., \code{expand.grid}) of X, M, and I values as initial values for which I was calculating equilibria.
#' @seealso \code{\link{getEigs}}
#' @return data.frame of eigenvalues
#' @export
getEigs_df <- function(rootGrid_I, ...){
	rootGrid_I_comp <- data.matrix(rootGrid_I) #data.matrix(rootGrid_I[complete.cases(rootGrid_I),])
	eigGrid <- t(apply(rootGrid_I_comp,1,getEigs, ...))
	eigGrid_I <- data.frame(rootGrid_I_comp, eigGrid)
	return(eigGrid_I)
}


#' Stability Classification
#' 
#' Find equilibria for a model system of lake eutrophication and classify their stability
#' 
#' @param I numeric scalar representing P input; reasonable values might be between 0 and 2
#' @param pars a named numeric vector of optional additional parameter values to be passed to \code{func}
#' @param func function whose stability is to be classified; nominally this is \code{mDXM}, a wrapper for \code{\link{modelDeterministicXM}}; currently this function only works with the default b/c \code{\link{getRoot_df}} hasn't been set up to accept other functions as an argument
#' @param Xvals,Mvals numeric vector of X (water P) or M (sediment P) values over which to search for equilibria
#' @return a data.frame with a row for each equilibrium point found, and columns for the X,M coordinates of those points, a character describing its stability classification, the trace (tr) and determinant (Delta) of the Jacobian matrix at that point, the 'discriminant' value (tr^2 - 4*delta), and the parameter values supplied through I and pars.
#' @export
stabClass <- function(I, pars, func=mDXM, Xvals=seq(0,15,length.out=15), Mvals=seq(0,1E3,length.out=15)){
	requireNamespace("phaseR", quietly=TRUE)
	if(missing(pars)){pars <- NULL}
	if("I"%in%names(pars)){pars <- pars[names(pars)!="I"]}
	gridVals <- cbind(expand.grid(X=Xvals, M=Mvals), I=I)
	rs <- getRoot_df(gridVals, pars=pars)
	rs <- rs[complete.cases(rs),,drop=FALSE]
	urs <- rs[!duplicated(paste0(round(rs[,"X"],4),round(rs[,"M"],4))), c("I","X","M"), drop=FALSE]

	stab <- function(x){
		st <- phaseR::stability(mDXM, y.star=x[-1,drop=FALSE], parameters=c(I=I,pars), summary=FALSE)
		# if(is.null(names(st$parameters))){names(st$parameters) <- "I"}
		o <- cbind(data.frame(X=st$y.star[1], M=st$y.star[2], classification=st$classification, Delta=st$Delta, discriminant=st$discriminant, tr=st$tr),as.list(st$parameters))
		rownames(o) <- NULL
		o$classification <- as.character(o$classification)
		return(o)
	}
	do.call('rbind', apply(urs, 1, stab))
}

#' Simulate Lake Phosphorus
#' 
#' Wrapper function for simulating a bifurcation involving lake phosphorus
#' 
#' @param nYears integer, number of years for simulation (not necessarily n time steps)
#' @param I_range numeric vector of length 2, the starting and stopping values of P inputs to be used in the simulation; will increment linearly over the course of the simulation (including changing within a time step)
#' @param pars named numeric vector of possible additional parameters to pass to \code{\link{modelDeterministicXM}}
# @param agg_steps numeric vector; should be at least length 1. Currently doesn't do anything to the simulation, but is recorded in output, wich can be convenient later (should probably remove this argument)
# @param steps_per_day numeric vector; default is to calculate as 24/agg_steps; when the first value of agg_steps is 1, then the first default value of steps_per_day will be 24, which means that 1 day (year) of the simulation will have 24 observations
#' @param dt numeric vector of length 1; should be <= 1. The fraction of a year at which a simulation should be incremented. For example, 1/24 (the default) would indicate 24 samples per year.
#' @param add_sin logical; if TRUE (default), a sine wave is imposed on the process of the simulation
#' @param sin_amp numeric vector of length 2; amplitude of the sine wave; first value is the amplitude for water P, second is sediment P
#' @param noise_coeff numeric vector of length 2; standard deviation of the process error for water P (first value) and sediment P (second value)
#' 
#' @return a list of length 2; first element is a matrix of state variables, second is the dt value
#' @export
# simP_ts <- function(nYears, I_range=c(1.0, 1.33), agg_steps=c(1,4,24,48), steps_per_day=24/agg_steps, dt=1/steps_per_day[1], add_sin=TRUE, sin_amp=c(0.0075, 0.005), noise_coeff=c(0.005, 0.1)){
simP_ts <- function(nYears, I_range=c(1.0, 1.33), pars=NULL, dt=1/24, add_sin=TRUE, sin_amp=c(0.0075, 0.005), noise_coeff=c(0.005, 0.1)){
	stateMat <- matrix(NA, nrow=nYears/dt, ncol=2, dimnames=list(NULL, c("X","M"))) # empty matrix to hold state variables
	I <- seq(from=I_range[1], to=I_range[2], length.out=nYears/dt) #0.25 #1.5 is a good value to show for simulation, maybe # nominal fraction of soil P washed into the lake; i think the critical point is at just over 1.336
	
	# stateMat[1,] <- c(1.5, 200) # set initial values for time series simulation
# 	stateMat[1,] <- getRoot(c(I=I[1], stateMat[1,"X"], stateMat[1,"M"]), pars) # make the initial values at equilibrium
# 	if(any(is.na(stateMat[1,])) | any(stateMat[1,]<0)){ # if getRoot returns NA, it didn't find a root in max.iter steps; probably because are parameters are such that it'd be faster to start at higher values of the state parameter
# 		stateMat[1,] <- c(5, 600) # set initial values for time series simulation
# 		stateMat[1,] <- getRoot(c(I=I[1], stateMat[1,"X"], stateMat[1,"M"]), pars) # make the initial values at equilibrium
# 	}
	stab <- stabClass(I=I[1], pars=pars)
	stab <- stab[stab$X>0 & stab$M>0,]
	stable <- grepl("Stable", stab$classification)
	if(any(stable)){
		init <- stab[stable,]
		init <- init[which.min(init$X),]
	}else{
		init <- stab[which.min(stab$X),]
	}
	stateMat[1,] <- unlist(init[1,c("X","M")])
	
	
	for(j in 2:nrow(stateMat)){ # iterate through time steps
		state <- stateMat[j-1,]
		# if((j%%(288))==0){state <- getRoot(c(state, I=I[j-1]))} # this was intended to tie the system to the equilibrium, but don't do this, because it messes up the statistics. It'd be okay if set to equilibrium, simulated at constant I, and calculated ar() for that period of constant I. But changing I within a window of data for which ar() is calculated, while also abruptly setting to equilibrium within that same window ... doing that really messes stuff up. It creates decreasing autocorrelation at some time scales ... I think because of the additional jumps (depending how often it is tied back to root).
		dState_dt <- modelDeterministicXM(state=state, pars=c(I=(I[j]), pars))
	
		# Runge-Kutta Approximation
		# Not sure if I did this correctly
		# k1 <- dState_dt
		# k2 <- modelDeterministic(state=state+dt*k1/2)
		# k3 <- modelDeterministic(state=state+dt*k2/2)
		# k4 <- modelDeterministic(state=state+dt*k3)
		# rkState <- stateMat[j-1,] + dt/6*(k1+2*k2+2*k3+k4)
		# stateMat[j,] <- rkState # runge kutta
	
		# Euler Method Approximation
		dState <- dState_dt*dt + rnorm(2, sd=c(noise_coeff[1], noise_coeff[2]))*dt + c(sin_amp[1], sin_amp[2])*cos(2*pi*dt*j)*(dt*2*pi)
		eulerState <- state + dState
		stateMat[j,] <- eulerState # euler
	
	}
	stateMat <- cbind(time=seq(0, nYears-dt, by=dt), I=I, stateMat)
	# return(list(stateMat=stateMat, agg_steps=agg_steps, steps_per_day=steps_per_day, dt=dt))
	return(list(stateMat=stateMat, dt=dt))
}


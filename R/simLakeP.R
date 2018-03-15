# =====================
# = Packages Required =
# =====================
# You need to install.packages("rootSolve")


# ===================================
# = Draft of Function Documentation =
# ===================================
#' Eutrophication Bifurcation
#' Model of lake eutrophication via elevated phosphorus loading. When as phosphorus loading increases, a fold bifurcation is reached, and the system switches into a eutrophic state from a oligotrophic state. 
#' @param F input rate of phosporus to soil
#' @param C coefficient for transfer of soil phosphorus to the lake
#' @param H noise for input to the lake
#' @param s sedimentation loss
#' @param h hydrologic loss (outflow)
#' @param r recycling coefficient
#' @param sigma standard deviation of recycling noise
#' @param b permanent burial rate of phosphorus in sediments
#' @param dW a white noise process with mean zero and variance dt
#' @param m in the recycling function, the value of X at which recycling is half the maximum rate
#' @param q in the recycling function, the exponent q determines the slope of R(X) near m


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
modelDeterministic <- function(state, pars=c(C=0.00115), F=14.6, b=0.001, h=0.15, m=2.4, r=0.019, q=8, s=0.7, sigma=0.01, lambda=0.35){
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
modelDeterministicXM <- function(state, pars=c(I=0.00115*750), b=0.001, h=0.252, m=2.4, r=0.019, q=10, s=0.748){
	with(as.list(c(state, pars)),{ # use with() so that I can refer to column names in state and names in pars as variables
		R <- (X^q)/(m^q + X^q) # recycling
		MRX <- M*R
		rMRX <- r*MRX
		dX_dt <- I - (s+h)*X + rMRX # water P dynamics
		dM_dt <- s*X - b*M - rMRX # sediment P dynamics
		c(dX_dt=dX_dt, dM_dt=dM_dt) # output rate of change for X and M
	}) # end with
}
#' Jacobian matrix for deterministic eutrophication model
#' Jacboian for the 2D eutrophication model involving water, mud, and input phosphorus
#' @param t time, for compatability with rootSolve::jacobian.full, but not used otherwise
#' @param state state variables as in \code{\link{modelDeterministicXM}}
#' @param pars paramters as in \code{\link{modelDeterministicXM}}
#' @return Jacobian matrix
#' @export
mDXM_jac <- function(t, state, pars){
	modelDeterministicXM(state, pars)
}


# ========================
# = Simulate Time Series =
# ========================
nYears <- 300 # number of years for simulation
dt <- 1/36 # number of time steps per year
# stateMat <- matrix(NA, nrow=nYears/dt, ncol=3, dimnames=list(NULL, c("U","X","M"))) # empty matrix to hold state variables
stateMat <- matrix(NA, nrow=nYears/dt, ncol=2, dimnames=list(NULL, c("X","M"))) # empty matrix to hold state variables

# C <- 0.00115 # nominal fraction of soil P washed into the lake
I <- 1.75 #1.5 is a good value to show for simulation, maybe # nominal fraction of soil P washed into the lake
stateMat[1,] <- c(1.2, 100) # set initial values for time series simulation

for(j in 2:nrow(stateMat)){ # iterate through time steps
	state <- stateMat[j-1,]
	dState_dt <- modelDeterministicXM(state=state, pars=c(I=I))
	
	# Runge-Kutta Approximation
	# Not sure if I did this correctly
	# k1 <- dState_dt
	# k2 <- modelDeterministic(state=state+dt*k1/2)
	# k3 <- modelDeterministic(state=state+dt*k2/2)
	# k4 <- modelDeterministic(state=state+dt*k3)
	# rkState <- stateMat[j-1,] + dt/6*(k1+2*k2+2*k3+k4)
	# stateMat[j,] <- rkState # runge kutta
	
	# Euler Method Approximation
	dState <- dState_dt*dt # + rnorm(2)*dt
	eulerState <- state + dState
	stateMat[j,] <- eulerState # euler
	
}
par(mfrow=c(2,1), mar=c(1.85,1.85,0.5,0.5), mgp=c(1,0.25,0), tcl=-0.25, ps=8)
plot(ts(stateMat[,"X"], freq=1/dt), ylab="X, water P")
plot(ts(stateMat[,"M"], freq=1/dt), ylab="M, sediment P")
stateMat[nrow(stateMat),]


#' Get Initial Values for Eutrohpication Model
#' Gets a grid/ mesh of values for eventual input to model
#' @param gridN number of values for I and X
#' @param Irange range of values for I, input rate of P to lake
#' @param Xrange range of values for X, water phosphorus
#' @param Mrange range of values for M, mud P
#' @param In number of values for I; default is set to gridN
#' @param Xn number of values for X; default is set to gridN
#' @param Mn number of values for M; default is set to gridN
#' @return a data.frame of 
getInit <- function(gridN=100, Irange=c(0.25, 1.75), Xrange=c(0.05, 8), Mrange=c(400, 400), In=gridN, Xn=gridN, Mn=gridN){
	# gridN <- 100
	Igrid <- seq(Irange[1], Irange[2], length.out=In)
	Xgrid <- seq(Xrange[1], Xrange[2], length.out=Xn)
	Mgrid <- seq(Mrange[1], Mrange[2], length.out=Mn)

	sequences <- data.frame(I=Igrid, X=Xgrid, M=Mgrid)
	initialValues <- expand.grid(I=Igrid, X=Xgrid, M=Mgrid)
	# startMudP <-
	initialValues <- initialValues[order(initialValues[,"I"]),]
	rownames(initialValues) <- NULL
	initialValues
}

#' Get the Root of the Eutrophication Model
#' Finds the roots of the eutrophication model given starting values (water, mud P) and parameter (P loading)
#' @param x a vector of length 3 with names X (water P), M (mud P), and I (P loading)
#' @return numeric vector with roots
getRoot <- function(x){
	tryCatch({
		rootSolve::multiroot(f=modelDeterministicXM, start=x[c("X","M")], parms=list(I=x["I"]), maxiter=1E2)$root
	}, warning=function(w)c(X=NA,M=NA))
}

#' Get the Root of the Eutrophication Model from a Data Frame
#' Finds the roots of the eutrophication model given starting values (water, mud P) and parameter (P loading)
#' @param x a \code{data.frame} with 3 columns named X (water P), M (mud P), and I (P loading)
#' @return a matrix with columns for initial values of X and M, roots (of X and M), and the bifurcation parameter I. Each row of output corresponds to a row of the input (though output reorded to ascending I)
#' @export
getRoot_df <- function(initalValues){
	rootGrid <- t(apply(initialValues, 1, getRoot))
	rootGrid_I <- data.matrix(data.frame(I=initialValues[,"I"], init=initialValues[,c("X","M")], rootGrid))
	rootGrid_I <- rootGrid_I[order(rootGrid_I[,"I"]),]
	return(rootGrid_I)
}

#' Get Eigenvalues for Eutrophication Model
#' Given the rate of change of state variables (water and mud P) and the input P, returns the eigenvalues
#' @param x a named vector of length 3 (X, M, I), as in \code{\link{modelDeterministicXM}}
#' @return eigenvalues
getEigs <- function(x){
	if(any(is.na(x))){return(c(NA,NA))}
	jac <- rootSolve::jacobian.full(y=x[c("X","M")], func=mDXM_jac, parms=x["I"])
	eigen(jac)$values
}
#' Get Eigenvalues for Eutrophication Model (from data frame of states)
#' Given the rate of change of state variables (water and mud P) and the input P, returns the eigenvalues
#' @param x a matrix with columns named for X, M, I, as in \code{\link{modelDeterministicXM}}
#' @seealso \code{\link{getEigs}}
#' @return data.frame of eigenvalues
#' @export
getEigs_df <- function(rootGrid_I){
	rootGrid_I_comp <- data.matrix(rootGrid_I) #data.matrix(rootGrid_I[complete.cases(rootGrid_I),])
	eigGrid <- t(apply(rootGrid_I_comp,1,getEigs))
	eigGrid_I <- data.frame(rootGrid_I_comp, eigGrid)
	return(eigGrid_I)
}


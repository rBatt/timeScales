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

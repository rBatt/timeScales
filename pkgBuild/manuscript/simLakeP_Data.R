

# ========================
# = Simulate Time Series =
# ========================
nYears <- 200 # number of years for simulation
dt <- 1/(288) # number of time steps per year
# stateMat <- matrix(NA, nrow=nYears/dt, ncol=3, dimnames=list(NULL, c("U","X","M"))) # empty matrix to hold state variables
stateMat <- matrix(NA, nrow=nYears/dt, ncol=2, dimnames=list(NULL, c("X","M"))) # empty matrix to hold state variables

# C <- 0.00115 # nominal fraction of soil P washed into the lake
I <- seq(from=1, to=1.42, length.out=nYears/dt) #0.25 #1.5 is a good value to show for simulation, maybe # nominal fraction of soil P washed into the lake; i think the critical point is ~1.34
stateMat[1,] <- c(1.2, 200) # set initial values for time series simulation
stateMat[1,] <- getRoot(c(I=I[1], stateMat[1,"X"], stateMat[1,"M"])) # make the initial values at equilibrium

for(j in 2:nrow(stateMat)){ # iterate through time steps
	state <- stateMat[j-1,]
	dState_dt <- modelDeterministicXM(state=state, pars=c(I=I[j]))
	
	# Runge-Kutta Approximation
	# Not sure if I did this correctly
	# k1 <- dState_dt
	# k2 <- modelDeterministic(state=state+dt*k1/2)
	# k3 <- modelDeterministic(state=state+dt*k2/2)
	# k4 <- modelDeterministic(state=state+dt*k3)
	# rkState <- stateMat[j-1,] + dt/6*(k1+2*k2+2*k3+k4)
	# stateMat[j,] <- rkState # runge kutta
	
	# Euler Method Approximation
	dState <- dState_dt*dt + rnorm(2, sd=c(0.025, 0.05))*dt
	eulerState <- state + dState
	stateMat[j,] <- eulerState # euler
	
}
par(mfrow=c(2,1), mar=c(1.85,1.85,0.5,0.5), mgp=c(1,0.25,0), tcl=-0.25, ps=8)
plot(ts(stateMat[,"X"], freq=1/dt), ylab="X, water P")
plot(ts(stateMat[,"M"], freq=1/dt), ylab="M, sediment P")
stateMat[nrow(stateMat),]

# Record like a data set
lakeP <- data.table(lake="Sim", doy=time(ts(stateMat[,"X"], freq=1/dt)), I=I, X=stateMat[,"X"], M=stateMat[,"M"])
lakeP_df <- data.frame(lakeP[,list(I,X,M)]) # for finding eigenvalues
lakePm <- melt(lakeP, id.vars=c("lake","doy"))#[variable%in%vars & lake%in%lakes]


# ==========================
# = Test Stats on Data Set =
# ==========================
# set_ts <- function(y, x, freq=288){
# 	ts(y, freq=288, start=x)
# }
lakePm[, value:=set_ts(y=log(value), x=doy[1]), by=c("lake","variable")]
agg_sos2 <- function(aggsteps){
	out <- lakePm[,j={agg_ts(y=value, x=doy, width=aggsteps)},by=c("lake","variable")]
	out
}
lakeP_agg <- lapply(agg_steps, FUN=agg_sos2)
names(lakeP_agg) <- paste0("agg", agg_steps)
plotac_simp <- function(X, ...){
	X <- copy(X)
	ylim <- X[,range(y, na.rm=TRUE)]
	
	ul <- X[,unique(lake)]
	for(l in 1:length(ul)){
		dud <- X[lake==ul[l],j={
			tcol <- c("Sim"="forestgreen", "Paul"="blue","Peter"="red", "zdiff"="black")[lake[1]]
			plot(x,y, type='l', col=tcol, ylim=ylim, ...)
			NULL
		}]
	}
	
	invisible()
}

simAC <- roll_ac.sos(lakeP_agg, window_elapsed=steps_per_window, vars="X", lakes="Sim", DETREND=TRUE, by=window_by)
png("~/Desktop/simP_rollingAC_Detrend.png", width=3.5, height=5.5, res=150, units='in')
par(mfrow=c(4,1), mar=c(2,2,0.5,0.5), ps=8, mgp=c(1,0.25, 0), tcl=-0.15)
mapply(plotac_simp, simAC)
dev.off()

simAC_noDetrend <- roll_ac.sos(lakeP_agg, window_elapsed=steps_per_window, vars="X", lakes="Sim", DETREND=FALSE, by=window_by)
png("~/Desktop/simP_rollingAC_noDetrend.png", width=3.5, height=5.5, res=150, units='in')
par(mfrow=c(4,1), mar=c(2,2,0.5,0.5), ps=8, mgp=c(1,0.25, 0), tcl=-0.15)
mapply(plotac_simp, simAC_noDetrend)
dev.off()




simP_eigs <- getEigs_df(lakeP_df)
ind <- simP_eigs$I>1.0 & simP_eigs$I<1.38
# plot(simP_eigs$I[ind], Re(simP_eigs$X1)[ind])

# first time eig goes to 0, using observed state variables
# but careful, because theory talks about the Jacobian near equilibrium, which the simulation may not be
firstUnstable <- which(Re(simP_eigs$X1)>=0)[1] # when eigenvalue is 0 based on sim values
lakeP_df[firstUnstable,]
fU_XMroot <- getRoot(unlist(lakeP_df[firstUnstable,]))
firstUnstable_eigs <- getEigs(c(I=simP_eigs$I[firstUnstable], fU_XMroot)) # when sim'd values give eigenvalue of 0, the eigenvalue near the corresponding equilibrium (root of that sim'd state) indicates the system is still stable.

# find the roots at each time step, using simulated values as initial values
# then use these roots to calculate return times
rootOfSim <- getRoot_df(lakeP_df)
rootEigs <- getEigs_df(rootOfSim[,c("I","X","M")])
print(rootEigs[which(Re(rootEigs$X1)>=0)[1],], digits=10) # critical value based on when eig first crosses 0
print(rootEigs[which.max(Re(rootEigs$X1)),], digits=10) # critical value based on maximum eig

png("~/Desktop/simP_eigenvalues.png", width=3.5, height=3.5, res=150, units='in')
par(mar=c(2,2,0.5,0.5), ps=8, mgp=c(1,0.25, 0), tcl=-0.15)
plot(rootEigs$I, Re(rootEigs$X1), xlab="Nutrient Loading", ylab="Eigenvalue", type='l') # plot of eigenvalues across values of I (loading)
dev.off()


lakeP[,c("root.X","root.M", "X1", "X2"):=list(rootEigs$X, rootEigs$M, rootEigs$X1, rootEigs$X2)]


# hmm, the system is very far from equilibrium
# it's especially strange how the simulation P way overshoots the equilibrium value
# in other simulations, i don't recall
png("~/Desktop/simP_simValues_trueEquilibria.png", width=3.5, height=3.5, res=150, units='in')
par(mfrow=c(2,1), mar=c(2,2,0.5,0.5), ps=8, mgp=c(1,0.25, 0), tcl=-0.15)
lakeP[,plot(as.numeric(doy), X, type='l', ylim=range(c(X, root.X), na.rm=TRUE), xlab="Day", ylab="Water P")]
lakeP[,lines(as.numeric(doy), root.X, lty=2, col='forestgreen')]
lakeP[,plot(as.numeric(doy), M, type='l', ylim=range(c(M, root.M), na.rm=TRUE), xlab="Day", ylab="Mud P")]
lakeP[,lines(as.numeric(doy), root.M, lty=2, col='forestgreen')]
dev.off()


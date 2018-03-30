
# ========================
# = Simulate Time Series =
# ========================
library(timeScales)
library(rootSolve)
library(R2jags)
library(foreach)
library(doParallel)

# ================
# = Set RNG Seed =
# ================
set.seed(42)


# ==============================================
# = Set Up Options for Simulation and Analysis =
# ==============================================
agg_steps <- c(1, 4, 24, 48) #c(1, 12, 288, 288*2) # step sizes for aggregation
steps_per_day <- 24/(agg_steps) #60*24/(5 * agg_steps) # obs per day = (60 min / 1 hr) * (24 hrs / 1 day) * (1 obs / 5*n min)

win_days <- 28 #28 # window size in days covered
steps_per_window <- steps_per_day*win_days # steps per window = (n steps / 1 day) * (n days / 1 window)
window_by <- 1 #pmax(1, steps_per_day/(6)) #pmax(1, steps_per_day/(4)) # the denominator is number of window starts per day; if trying to increment window by less than the resolution of the time series, just increment by 1 #c(48, 4, 2, 1)

# out <- simP_ts(nYears=100, I_range=c(1.1, 1.4), agg_steps=agg_steps, steps_per_day=steps_per_day, dt=1/steps_per_day[1])
out <- simP_ts(nYears=100, I_range=c(1.1, 1.4), dt=1/steps_per_day[1])
for(n in names(out)){ # make all the elements of the output list objects in the current env
	assign(n, out[[n]])
}

# Record like a data set
lakeP <- data.table(lake="Sim", doy=time(ts(stateMat[,"X"], freq=1/dt)), I=stateMat[,"I"], X=stateMat[,"X"], M=stateMat[,"M"])
lakeP_df <- data.frame(lakeP[,list(I,X,M)]) # for finding eigenvalues
lakePm <- melt(lakeP, id.vars=c("lake","doy"))#[variable%in%vars & lake%in%lakes]
# save(lakePm, file="~/Documents/School&Work/epaPost/timeScales/data/lakePm.RData", compress="xz")


# ===================================
# = Calculate Eigenvalues and Roots =
# ===================================
lakeP_df <- as.data.frame(dcast(lakePm, formula=lake+doy~variable)[,list(I,X,M)]) # just here to show how to create from what I intend to save
simP_eigs <- getEigs_df(lakeP_df)

# find the roots at each time step, using simulated values as initial values
# then use these roots to calculate return times
rootOfSim <- getRoot_df(lakeP_df)
rootEigs <- getEigs_df(rootOfSim[,c("I","X","M")])
lakeP[,c("root.X","root.M", "X1", "X2"):=list(rootEigs$X, rootEigs$M, rootEigs$X1, rootEigs$X2)]
print(rootEigs[which(Re(rootEigs$X1)>=0)[1],], digits=10) # critical value based on when eig first crosses 0
print(rootEigs[which.max(Re(rootEigs$X1)),], digits=10) # critical value based on maximum eig


# ===================
# = Plot Simulation =
# ===================
dev.new(width=3.5, height=4.5)
par(mfrow=c(2,1), mar=c(1.85,1.85,0.5,0.5), mgp=c(1,0.25,0), tcl=-0.25, ps=8, cex=1)
lakePm[variable=="X", plot(ts(value, freq=1/dt), xlab="X, water P")]
lakePm[variable=="M", plot(ts(value, freq=1/dt), xlab="M, sediment P")]


# ===============================
# = Plot Roots and Observations =
# ===============================
# hmm, the system is very far from equilibrium
# it's especially strange how the simulation P way overshoots the equilibrium value
# in other simulations, i don't recall
dev.new(width=3.5, height=4.5)
# png("~/Desktop/simP_simValues_trueEquilibria.png", width=3.5, height=3.5, res=150, units='in')
par(mfrow=c(2,1), mar=c(2,2,0.5,0.5), ps=8, mgp=c(1,0.25, 0), tcl=-0.15, cex=1)
lakeP[,plot(as.numeric(doy), X, type='l', ylim=range(c(X, root.X), na.rm=TRUE), xlab="Day", ylab="Water P")]
lakeP[,lines(as.numeric(doy), root.X, lty=2, col='forestgreen')]
lakeP[,plot(as.numeric(doy), M, type='l', ylim=range(c(M, root.M), na.rm=TRUE), xlab="Day", ylab="Mud P")]
lakeP[,lines(as.numeric(doy), root.M, lty=2, col='forestgreen')]
# dev.off()


# ===============================
# = Plot Simulation Eigenvalues =
# ===============================
dev.new(width=3.5, height=3.5)
par(mar=c(2,2,2,0.5), ps=8, mgp=c(1,0.25, 0), tcl=-0.15)
plot(simP_eigs$I, Mod(simP_eigs$X1), type='l', xlab="Nutrient Loading", ylab="Eigenvalues around observed states")
par(new=TRUE)
plot(ts(Mod(simP_eigs$X1), freq=1/dt), type='n', xlab="", ylab="", xaxt='n', yaxt='n')
axis(side=3)
mtext("Time", side=3, line=1.25)

# # first time eig goes to 0, using observed state variables
# # but careful, because theory talks about the Jacobian near equilibrium, which the simulation may not be
# firstUnstable <- which(Re(simP_eigs$X1)>=0)[1] # when eigenvalue is 0 based on sim values
# lakeP_df[firstUnstable,]
# fU_XMroot <- getRoot(unlist(lakeP_df[firstUnstable,]))
# firstUnstable_eigs <- getEigs(c(I=simP_eigs$I[firstUnstable], fU_XMroot)) # when sim'd values give eigenvalue of 0, the eigenvalue near the corresponding equilibrium (root of that sim'd state) indicates the system is still stable.


# =================================================
# = Plot 'True' Eigenvalues of Root of Simulation =
# =================================================
# png("~/Desktop/simP_eigenvalues.png", width=3.5, height=3.5, res=150, units='in')
dev.new(width=3.5, height=3.5)
par(mar=c(2,2,2,0.5), ps=8, mgp=c(1,0.25, 0), tcl=-0.15)
plot(rootEigs$I, Re(rootEigs$X1), xlab="Nutrient Loading", ylab="Eigenvalue around Root", type='l') # plot of eigenvalues across values of I (loading)
par(new=TRUE)
plot(ts(Mod(rootEigs$X1), freq=1/dt), type='n', xlab="", ylab="", xaxt='n', yaxt='n')
axis(side=3)
mtext("Time", side=3, line=1.25)
# dev.off()


# =====================================
# = Statistics and Plotting Functions =
# =====================================
set_ts <- function(y, x, freq=288){
	ts(y, freq=freq, start=x)
}
agg_sos2 <- function(aggsteps){
	out <- lakePm[,j={agg_ts(y=value, x=doy, width=aggsteps)},by=c("lake","variable")]
	out
}
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


# ====================================================
# = Prepare Simulation Data for Statistical Analysis =
# ====================================================
lakePm[, value:=set_ts(y=log(value), x=doy[1], freq=1/dt), by=c("lake","variable")]

# ---- aggregate at different time scales ----
lakeP_agg <- lapply(agg_steps, FUN=agg_sos2)
names(lakeP_agg) <- paste0("agg", agg_steps)


# ====================================================================
# = Rolling Window AR(1); 4 Time Scales, Detrended and No Detrending =
# ====================================================================
# ---- Calculate AR(1) for both detrending and not ----
simAC <- roll_ac.sos(lakeP_agg, window_elapsed=steps_per_window, vars="X", lakes="Sim", DETREND=TRUE, by=window_by, save_output=F)
simAC_noDetrend <- roll_ac.sos(lakeP_agg, window_elapsed=steps_per_window, vars="X", lakes="Sim", DETREND=FALSE, by=window_by)

# png("~/Desktop/simP_rollingAC_Detrend.png", width=3.5, height=5.5, res=150, units='in')
ylabs <- paste0("AR(1) Coeff (", gsub("agg", "agg=", names(lakeP_agg)), ")")
dev.new(width=4, height=6)
par(mfcol=c(4,2), mar=c(2,2,0.5,0.5), ps=8, mgp=c(1,0.25, 0), tcl=-0.15, oma=c(0,0,1.0,0), cex=1)
mapply(plotac_simp, simAC, ylab=ylabs)
mtext("Detrended", side=3, outer=TRUE, adj=0.25, line=0)
mapply(plotac_simp, simAC_noDetrend, ylab=ylabs)
mtext("No Detrending", side=3, outer=TRUE, adj=0.85, line=0)
# dev.off()


# ===============================
# = AR(p) Model: Rolling Window =
# ===============================
simARp <- roll_ac.sos(lakeP_agg, fit_arP=TRUE, window_elapsed=steps_per_window, vars="X", lakes="Sim", DETREND=TRUE, by=window_by, save_output=F)
simARp_nod <- roll_ac.sos(lakeP_agg, fit_arP=TRUE, window_elapsed=steps_per_window, vars="X", lakes="Sim", DETREND=FALSE, by=window_by, save_output=F)

ylabs <- paste0("||Eig|| of AR(p) (", gsub("agg", "agg=", names(lakeP_agg)), ")")
dev.new(width=4, height=6)
par(mfcol=c(4,2), mar=c(2,2,0.5,0.5), ps=8, mgp=c(1,0.25, 0), tcl=-0.15, oma=c(0,0,1.0,0), cex=1)
mapply(plotac_simp, simARp, ylab=ylabs)
mtext("Detrended", side=3, outer=TRUE, adj=0.25, line=0)
mapply(plotac_simp, simARp_nod, ylab=ylabs)
mtext("No Detrending", side=3, outer=TRUE, adj=0.75, line=0)


# =================
# = DLM TVAR(1)SS =
# =================
tvarss_xvals <- lapply(lakeP_agg, function(x){x[variable=="X", as.numeric(x)]})


tvdlm_wrapper <- function(x, det=FALSE, fit_arP=FALSE, ...){
	
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
	tvarss_list <- mapply(tvarss, ylist, nP=nPs, MoreArgs=list(niter=2E3, parallel=TRUE, oType="tvarss", ...), SIMPLIFY=FALSE)
	# names(tvarss_list) <- names(x) # not needed, will keep names
	
	if(fit_arP){
		if(any(nPs>1)){
			requireNamespace("doParallel", quiety=TRUE)
			requireNamespace("foreach", quiety=TRUE)
			
			doParallel::registerDoParallel(cores=4)
			eigs <- foreach::foreach(j in 1:length(tvarss_list), .combine=list) %dopar% {
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

# ===================
# = Fit TVARSS DLMs =
# ===================
# ---- TVAR(1)SS ----
tvarss_noMean_list <- tvdlm_wrapper(lakeP_agg, det=TRUE)
tvarss_noMean_nod_list <- tvdlm_wrapper(lakeP_agg)

tvarss_list <- tvdlm_wrapper(lakeP_agg, tvMean=TRUE, det=TRUE)
tvarss_nod_list <- tvdlm_wrapper(lakeP_agg, tvMean=TRUE)

# ---- Plot TVAR(1)SS ----
ylabs <- paste0("TVAR(1) of DLM (", gsub("agg", "agg=", names(lakeP_agg)), ")")

# plot version with a 'constant' mean parameter
dev.new(width=4, height=6)
par(mfcol=c(4,2), mar=c(2,2,0.5,0.5), ps=8, mgp=c(1,0.25, 0), tcl=-0.15, oma=c(0,0,1.5,0), cex=1)
ppout <- mapply(plotPost.tvarss, x=tvarss_noMean_list, xvals=tvarss_xvals, ylab=ylabs, MoreArgs=list(xlab="Time", main="", varName="Phi"))
mtext("Constant Mean, Detrended", side=3, line=0, outer=TRUE, adj=0.85)
ppout <- mapply(plotPost.tvarss, x=tvarss_noMean_nod_list, xvals=tvarss_xvals, ylab=ylabs, MoreArgs=list(xlab="Time", main="", varName="Phi"))
mtext("Constant Mean, No Detrending", side=3, line=0, outer=TRUE, adj=0.15)

# plot version with a time-varying mean parameter
dev.new(width=4, height=6)
par(mfcol=c(4,2), mar=c(2,2,0.5,0.5), ps=8, mgp=c(1,0.25, 0), tcl=-0.15, oma=c(0,0,1.5,0), cex=1)
ppout <- mapply(plotPost.tvarss, x=tvarss_list, xvals=tvarss_xvals, ylab=ylabs, MoreArgs=list(xlab="Time", main="", varName="Phi"))
mtext("Varying Mean, Detrended", side=3, line=0, outer=TRUE, adj=0.85)
ppout <- mapply(plotPost.tvarss, x=tvarss_nod_list, xvals=tvarss_xvals, ylab=ylabs, MoreArgs=list(xlab="Time", main="", varName="Phi"))
mtext("Varying Mean, No Detrending", side=3, line=0, outer=TRUE, adj=0.15)


# ---- TVAR(p)SS ----
tvarss_noMean_nP_list <- tvdlm_wrapper(lakeP_agg, fit_arP=TRUE, det=TRUE)
tvarss_noMean_nP_nod_list <- tvdlm_wrapper(lakeP_agg, fit_arP=TRUE)

tvarss_nP_list <- tvdlm_wrapper(lakeP_agg, tvMean=TRUE, fit_arP=TRUE, det=TRUE)
tvarss_nP_nod_list <- tvdlm_wrapper(lakeP_agg, tvMean=TRUE, fit_arP=TRUE)

# ---- Plot TVAR(p)SS ----
# ylabels for tvar(p)ss
ylabs_nM_nP_nod <- sapply(tvarss_noMean_nP_nod_list, function(x)dim(x$Phi)[3])
ylabs_nP <- sapply(tvarss_nP_list, function(x)dim(x$Phi)[3])
ylabs_nP_nod <- sapply(tvarss_nP_nod_list, function(x)dim(x$Phi)[3])

# plot version with a 'constant' mean parameter
dev.new(width=4, height=6)
par(mfcol=c(4,2), mar=c(1,2,0.5,0.5), ps=8, mgp=c(1,0.15, 0), tcl=-0.15, oma=c(1,0,1.5,0), cex=1)
ylabs <- paste0("Eig of TVAR(p=", ylabs_nM_nP, ");", gsub("agg", "agg=", names(lakeP_agg)))
ppout <- mapply(plotPost.tvarss, x=tvarss_noMean_nP_list, xvals=tvarss_xvals, ylab=ylabs, MoreArgs=list(xlab="", main="", varName="Eigen"))
mtext("Constant Mean, Detrended", side=3, line=0, outer=TRUE, adj=0.17)
mtext("Time", side=1, line=0.75)
ylabs <- paste0("Eig of TVAR(p=", ylabs_nM_nP_nod, ");", gsub("agg", "agg=", names(lakeP_agg)))
ppout <- mapply(plotPost.tvarss, x=tvarss_noMean_nP_nod_list, xvals=tvarss_xvals, ylab=ylabs, MoreArgs=list(xlab="", main="", varName="Eigen"))
mtext("Constant Mean, No Detrending", side=3, line=0, outer=TRUE, adj=0.95)
mtext("Time", side=1, line=0.75)

# plot version with a time-varying mean parameter
dev.new(width=4, height=6)
par(mfcol=c(4,2), mar=c(1,2,0.5,0.5), ps=8, mgp=c(1,0.15, 0), tcl=-0.15, oma=c(1,0,1.5,0), cex=1)
ylabs <- paste0("Eig of TVAR(p=", ylabs_nP, ");", gsub("agg", "agg=", names(lakeP_agg)))
ppout <- mapply(plotPost.tvarss, x=tvarss_nP_list, xvals=tvarss_xvals, ylab=ylabs, MoreArgs=list(xlab="", main="", varName="Eigen"))
mtext("Varying Mean, Detrended", side=3, line=0, outer=TRUE, adj=0.15)
mtext("Time", side=1, line=0.75)
ylabs <- paste0("Eig of TVAR(p=", ylabs_nP_nod, ");", gsub("agg", "agg=", names(lakeP_agg)))
ppout <- mapply(plotPost.tvarss, x=tvarss_nP_nod_list, xvals=tvarss_xvals, ylab=ylabs, MoreArgs=list(xlab="", main="", varName="Eigen"))
mtext("Varying Mean, No Detrending", side=3, line=0, outer=TRUE, adj=0.95)
mtext("Time", side=1, line=0.75)


# =========================================
# = AR(1) Time Varying Regressoin (TVReg) =
# =========================================
library(tvReg)
dX <- detrendR(ts(stateMat[,"X"], frequency=steps_per_day[1]), max_poly=6, max_fourier=3, returnType="resid")
tvmod <- tvAR(dX, type='none', bw=0.6)
tvmod_nod <- tvAR(ts(stateMat[,"X"], frequency=steps_per_day[1]), bw=0.6)

tvreg_ar1 <- c(NA, tvmod$tvcoef)
tvreg_ar1_nod <- c(NA, tvmod_nod$tvcoef[,1])



arp_mod_full <- ar(dX, order.max=steps_per_day[1]*2)
nP <- length(arp_mod_full$ar)
arp_mod_full_nod <- ar(ts(stateMat[,"X"], frequency=steps_per_day[1]), order.max=steps_per_day[1]*2)
nP_nod <- length(arp_mod_full$ar)

tvmodp <- tvAR(dX, type='none', bw=0.1, p=nP)
tvmodp_nod <- tvAR(stateMat[,"X"], type='none', bw=0.1, p=nP_nod)

alleigs <- function(x){
	max(Mod(arEigs(x)))
}
tvmodp_eigs <- apply(tvmodp$tvcoef, 1, alleigs)
tvmodp_nod_eigs <- apply(tvmodp_nod$tvcoef, 1, alleigs)


dev.new(width=5, height=5)
par(mfrow=c(2,2), mar=c(2,2,0.5,0.5), ps=8, mgp=c(1,0.25, 0), tcl=-0.15, oma=c(0,0,1.5,0), cex=1)
plot(ts(tvreg_ar1, freq=1/dt), xlab="Time", ylab="AR(1) Coeff (time-varying reg)")
mtext("Detrended", side=3, line=-0.5, outer=TRUE, adj=0.25)
plot(ts(tvreg_ar1_nod, freq=1/dt), xlab="Time", ylab="AR(1) Coeff (time-varying reg)")
mtext("No Detrending", side=3, line=-0.5, outer=TRUE, adj=0.75)
plot(ts(tvmodp_eigs, frequency=steps_per_day[1]), type='l', xlab="Time", ylab=paste0("Eig of AR(p=",nP,") (time-varying reg)"))
plot(ts(tvmodp_nod_eigs, frequency=steps_per_day[1]), type='l', xlab="Time", ylab=paste0("Eig of AR(p=",nP,") (time-varying reg)"))





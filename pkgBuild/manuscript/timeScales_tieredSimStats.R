#' ---
#' title: "Apply autocorrelation early warning statistic to a 'tiered' simulation"
#' author: "Ryan Batt"
#' date: "2018-04-09"
#' abstract: I have applied autocorrelation, as an early warning statistic, to field data and simulated time series. In the simulation, nutrient loading, the bifurcation parameter, was gradually inreased each time step (of size dt). I then applied DLM or rolling-window style statistics (autocorrelation) to the resultant time series. However, the time series simulation departed from equilibrium values, even when the bifurcationw as changed very slowly and the simulation very long. A different way to do the simulation, which will avoid the potential problem of departing from equilibrium, is to calculate the statistic for different time series initiated at equilibrium and only vary the bifurcation parameter across time series, not within each time series. I am referring to this approach as a "tiered simulation". The goal of this experiment involvign the tiered simulation is to determine what the autocorrelation (AR(1)) or eigenvalue of an AR(p) should look like at several time scales. This objective is identical to that of the analyses of the field data and the previous simulation experiment, but with the advantage the the interpretation of the results should be easier, reasons aforementioned.
#'       .
#' output:
#'   html_document:
#'     toc: true
#'     toc_depth: 4
#'     fig_caption: true
#'     theme: "readable"
#'     template: default
#'   pdf_document:
#'     toc: true
#'     toc_depth: 4
#'     template: latex-ryan.template
#'     fig_caption: yes
#' geometry: margin=1.0in
#' lineno: true
#' lineSpacing: false
#' titlesec: true
#' documentclass: article
#' placeins: true
#' ---


#+ setup, include=FALSE, echo=FALSE, cache=FALSE
# ==========================================
# = Record Time to Get Elapsed Time at End =
# ==========================================
t1 <- Sys.time()


# =================
# = Load Packages =
# =================
library("data.table")
library("zoo")
library("forecast")
library("timeScales")
library(rootSolve)
library(R2jags)
library(foreach)
library(doParallel)
library(tvReg)

# Report
library(knitr)
library(rmarkdown)


# ================
# = Report Setup =
# ================
doc_type <- c("html", "pdf")[1]
table_type <- c("html"="html", "pdf"="latex")[doc_type]
options("digits"=3) # rounding output to 4 in kable() (non-regression tables)
o_f <- paste(doc_type, "document", sep="_")

# problem with pdflatex in El Capitan? It might be directories. Check http://pages.uoregon.edu/koch/FixLink.pkg

# render!
# rmarkdown::render(
# 	"~/Documents/School&Work/epaPost/timeScales/pkgBuild/manuscript/timeScales_tieredSimStats.R",
# 	output_format=o_f,
# 	output_dir='~/Documents/School&Work/epaPost/timeScales/pkgBuild/manuscript',
# 	clean = TRUE
# )

Sys.setenv(PATH=paste0("/Library/TeX/texbin:",Sys.getenv("PATH")))
opts_chunk$set(
	fig.path = 'timeScales_tieredSimStats/', 
	cache.path='timeScales_tieredSimStats/',
	echo=TRUE, 
	include=TRUE, 
	cache=FALSE,
	autodep=TRUE,
	results='asis',
	warning=FALSE,
	fig.show="hold",
	fig.lp = if(o_f=="html_document"){"**Figure.**"}else{NULL}
)


#+ setSeed
# ================
# = Set RNG Seed =
# ================
set.seed(42)

#' #Simulation Options, Helper Functions
#+ setupSim
# ==============================================
# = Set Up Options for Simulation and Analysis =
# ==============================================
dt <- 1/24
agg_steps <- c(1, 4, 24, 48) #c(1, 12, 288, 288*2) # step sizes for aggregation
steps_per_day <- (1/dt)/(agg_steps) #60*24/(5 * agg_steps) # obs per day = (60 min / 1 hr) * (24 hrs / 1 day) * (1 obs / 5*n min)

win_days <- 28 #28 # window size in days covered
steps_per_window <- steps_per_day*win_days # steps per window = (n steps / 1 day) * (n days / 1 window)
window_by <- 1 #pmax(1, steps_per_day/(6)) #pmax(1, steps_per_day/(4)) # the denominator is number of window starts per day; if trying to increment window by less than the resolution of the time series, just increment by 1 #c(48, 4, 2, 1)

#+ helperFunctions
# =====================================
# = Statistics and Plotting Functions =
# =====================================
set_ts <- function(y, x, freq=288){
	ts(y, freq=freq, start=x)
}
agg_sos3 <- function(aggsteps){
	out <- tieredTSm[,j={agg_ts(y=value, x=timeVec, width=aggsteps)},by=c("lake","variable","I")]
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

# plot_aggStat <- function(x, stats=c("AR1","ARpEig")){
# 	for(a in 1:length(x)){
#
# 	}
# }



#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Time Series of Water P
#' ##Wrapper Function for Simulated the Tiered Time Series
#+ function-tieredTSWrapper
# Simulated a tiered (in P input) time series
# 
# Simulate a time series water P involving stepwise, or "tieried", increases of the P input
# 
# @param nTiers number of tiers, which is number of values of P input (I)
# @param nYears numeric, number of years per tier
# @param dt the per-year size of a time step
# @param ... further arguments to be passed to \code{\link{simP_ts}}; recommended arguments are q and sin_amp
# 
# @return returns a data.table containing time series of X, water phosphorus
# @export
tieredTSWrapper <- function(nTiers, nYears, iranges, dt=1/24, ...){
	ts_x <- seq(0, by=dt, length.out=nYears/dt*nTiers)
	tieredTS <- matrix(NA, nrow=nYears/dt, ncol=nTiers)

	for(j in 1:nTiers){
		tieredTS[,j] <- simP_ts(nYears=nYears, I_range=rep(iranges[j],2), dt=dt, ...)[[1]][,"X"]
	}
	tieredTS2 <- data.table(lake='Sim', timeVec=ts_x, I=rep(iranges, each=nYears/dt), X=c(tieredTS))
	return(tieredTS2)
}

#' ##Simulate Time Series
#+ tieredSimulation
nTiers <- 50
nYears <- 2E2 # every nYears, the nutrient loading is increased
iranges <- seq(0.01, 1.4, length.out=nTiers)
qVals <- c(2,5,8,10)
tiered_list <- lapply(qVals, function(x){tieredTSWrapper(nTiers=nTiers, nYears=nYears, iranges=iranges, dt=1/24, pars=c(q=x), sin_amp=c(0,0))})
names(tiered_list) <- paste0("q",qVals)

# ts_x <- seq(0, by=dt, length.out=nYears/dt*nTiers)
# tieredTS <- matrix(NA, nrow=nYears/dt, ncol=nTiers)
#
# for(j in 1:nTiers){
# 	tieredTS[,j] <- simP_ts(nYears=nYears, I_range=rep(iranges[j],2), dt=dt, sin_amp=c(0,0), q=11)[[1]][,"X"]
# }
# tieredTS2 <- data.table(lake='Sim', timeVec=ts_x, I=rep(iranges, each=nYears/dt), X=c(tieredTS))

#' ##Plot Time Series of I
#+ figure1-Its, fig.width=3.5, fig.height=3.5, fig.cap="**Figure 1.** Time series of phosphorus input to the lake. Note that the changes in phosphorus are 'tiered' (increase stepwise) over the course of the time series. Simulations of water P will be initiated at equilibrium at the start of each tier (see **Figure 2**). Simulations used 50 tiers of P input, with each tier having a duration of 200 years."
# dev.new()
par(mar=c(2,2,1,0.5), mgp=c(1,0.25,0), tcl=-0.25, cex=1, ps=8)
tiered_list[[1]][,plot(timeVec, I, xlab="Year", ylab="Phosphorus Input", type='l')]
grid()


#' ##Plot Time Series of Water P
#+ figure2-waterPts, fig.width=6, fig.height=6, fig.cap="**Figure 2.** Time series of water phosphorus (X) in the lake as P input (I) is increased in a tiered (stepwise) fashion (see **Figure 1**). Each panel was simulate using a different value for the exponent q. At the start of each tier of I, X (and sediment phosphorus, M, not shown) is initiated at equilibrium. No sinusoidal process has been added to the dynamics in these simulations. The simulations are stochastic, with very low noise. The size of each time step is 1/24 years."
# dev.new()
par(mfrow=c(2,2), mar=c(2,2,1,0.5), mgp=c(1,0.25,0), tcl=-0.25, cex=1, ps=8)
for(j in 1:length(tiered_list)){
	tl <- tiered_list[[j]]
	qv <- names(tiered_list)[j]
	tl[,plot(timeVec, X, xlab="Year", ylab="Water P", main=paste0("q = ",gsub("q","",qv)))]
}


#' ##Prepare Time Series for Analysis
#' Will first select a nominal value of q before calculating on time series.
tieredTSm <- melt(tiered_list[['q8']], measure.vars="X") # undo with: dcast(tieredTSm, ...~variable)
tieredTSm[, value:=set_ts(y=(value), x=timeVec[1], freq=1/dt), by=c("lake","variable","I")]
# tieredTS2[,plot(timeVec, X, type='l')]
tiered_agg <- lapply(agg_steps, agg_sos3)
names(tiered_agg) <- paste0('agg',agg_steps)


#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Statistics
#' ##Calculate Statistics
#+ calculate-statistics
tiered_agg_stat <- vector('list', length(tiered_agg))
names(tiered_agg_stat) <- names(tiered_agg)
for(a in 1:length(tiered_agg)){
	ta <- tiered_agg[[a]]
	tiered_agg_stat[[a]] <- ta[,j={
		ar1coef <- ac1(y)
		arpcoefs <- arP(y, oType="coefs")
		nP <- length(arpcoefs)
		arpeigs <- arEigs(arpcoefs)
		data.table(x=mean(x), AR1=ar1coef, nP=nP, ARpEig=max(Mod(arpeigs)))
		
	}, by=c("lake","variable","I")]
}

#' ##Plot AR(1) and ||Eig|| of AR(p) across Time Scales
#+ figure3-ar1-eigARp-timeScales, fig.width=3.5, fig.height=6, fig.cap="**Figure 3.** Time series of the AR(1) coefficient (black line, left vertical axis) and the leading eigenvalue computed from the coefficients of an AR(p) model (blue line, right vertical axis). Statistics were calculated from tie series of water phosphorus (**Figure 2**). The horizontal axis is the value of P input (I), which serves as a control parameter inducing a bifurcation. For each panel, the time series was aggregated to a different time scale before calculating statistics: in the top panel, no aggregation was performed; in the lower panels, the time series was aggregated by averaging to yield the indicated number of samples per year (fractions indicate a sampling frequency < 1 year). Statistics were calulated for each 'tier' of P input (see **Figure 1**); thus, statistics were only applied to sections of the water P time series with constant parameter values, with each section being initated at equilibrium water P (and M, sediment P; not shown)."
# dev.new(width=3.5, height=6)
par(mfrow=c(4, 1), mar=c(1.2,2.5,0.75,2), mgp=c(1,0.25,0), oma=c(1,0.1,0.5,0.1), tcl=-0.25, ps=8, cex=1)
for(a in 1:length(tiered_agg_stat)){
	ta <- tiered_agg_stat[[a]]
	nms <- names(tiered_agg_stat)[a]
	aggVal <- gsub("agg","",nms)
	ylab1 <- "AR(1) Coef"
	minP <- ta[,min(nP)]
	maxP <- ta[,max(nP)]
	ylab2 <- bquote(paste("AR(",paste(.(minP) <= p )<= .(maxP), ") ||Eig||"))
	# ylab2 <- bquote("AR(" .(minP) <= p <= .(maxP) ") ||Eig||") #"AR(p) ||Eig||"
	ta[,plot(I, AR1, type='l', xlab='', ylab=ylab1)]
	par(new=TRUE)
	ta[,plot(I, ARpEig, type='l', col='blue', xlab='', ylab='', xaxt='n', yaxt='n')]
	axis(side=4)
	mtext(ylab2, side=4, line=1)
	mtext(paste0("Scale=", 1/dt/as.numeric(aggVal), " obs. per yr"), side=2, line=1.75, font=2)
	if(a==1){
		legend('top', lty=1, col=c("black","blue"), legend=c("AR(1)","AR(p) ||Eig||"), ncol=2, bty='n', inset=-0.275, xpd=TRUE)
	}
}


# dev.new()
# xlims <- tiered_agg_stat[[1]][,range(I)]
# ylims <- range(sapply(tiered_agg_stat, function(x)x[,range(AR1)]))
# plot(NA, xlim=xlims, ylim=ylims, xlab="", ylab="")
# for(a in 1:length(tiered_agg_stat)){
# 	ta <- tiered_agg_stat[[a]]
# 	ta[,lines(I, AR1)]
# }






#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Session Information
#+ sessionInfo, results='markup'
difftime(Sys.time(), t1) # how long it took to run these models/ produce this report
Sys.Date()
sessionInfo()
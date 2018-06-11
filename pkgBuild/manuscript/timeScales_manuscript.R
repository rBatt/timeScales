#' ---
#' title: "The time scale of resilience loss: the effect of sampling frequency on an early warning statistic"
#' author: "Ryan Batt"
#' date: "2017-10-17"
#' abstract: |
#'        Many complex systems show abrupt shifts, even when changes drivers change smoothly. Such changes are often referred to as regime shifts, and can be caused by critical transitions. These regime shifts can be very difficult to predict using mechanistic models. However, prior to a critical transition, complex systems may exhibit critical slowing down, a dynamical phenomenon characteristic of a system losing resilience prior to a critical transition. Statistical indicators for detecting the loss of resilience are being developed as a possible tool for predicting regime shifts involving critical transitions. However, these indicators can be sensitive to the observation frequency of the system. Researchers have often relied on expert judgement and intuition for selecting reasonable sampling frequencies for data analysis. However, the issue has not been systematically evaluated for all indicators, especially those based on autocorrelation. Furthermore, modern sensor technologies provide numerous options for sampling frequencies. The potential influence on sampling frequency on statistical indicators, coupled with the availability of high frequency data, presents a problem for the study and application of statistical indicators of resilience loss for researchers and managers alike. Here we analyze the sensitivity of autocorrelation, a statistical indicator of resilience, in a lake undergoing a critical transition. We evaluate how sampling frequency affects the value of this statistic, its use as an early warning indicator, and propose a general method for selecting an appropriate sampling frequency in other systems.
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


#+ setup, include=FALSE, echo=FALSE
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
library(deSolve)
library(phaseR)

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
# 	"~/Documents/School&Work/epaPost/timeScales/pkgBuild/manuscript/timeScales_manuscript.R",
# 	output_format=o_f,
# 	output_dir='~/Documents/School&Work/epaPost/timeScales/pkgBuild/manuscript',
# 	clean = TRUE
# )

Sys.setenv(PATH=paste0("/Library/TeX/texbin:",Sys.getenv("PATH")))
opts_chunk$set(
	fig.path = 'timeScales_manuscript/', 
	cache.path='timeScales_manuscript/',
	echo=TRUE, 
	include=TRUE, 
	cache=F,
	autodep=TRUE,
	results='asis',
	warning=FALSE,
	fig.show="hold",
	fig.lp = if(o_f=="html_document"){"**Figure.**"}else{NULL}
)


#' #Options
#+ options
# ========================================
# = Options for Field and Simulated Data =
# ========================================
win_days <- 28 # window size in days covered


# ==========================
# = Options for Field Data =
# ==========================
agg_steps <- c(1, 12, 288, 288*2) # step sizes for aggregation
lakes <- c("Peter","Paul") # can be vector; lakes to analyze (Paul, Peter)
vars <- c("logChla") #c("logBga") # #"chla" #"bga" # can be vector; variables to analyze (wtr, bga, chla)
field_ylab <- c("chla"="Chlorophyll", "logChla"="Log10 Chlorophyll", "bga"="Phycocyanin", "logBga"="Log10 Phycocyanin")

steps_per_day <- 60*24/(5 * agg_steps) # obs per day = (60 min / 1 hr) * (24 hrs / 1 day) * (1 obs / 5*n min)
steps_per_window <- steps_per_day*win_days # steps per window = (n steps / 1 day) * (n days / 1 window)
acf_lag.max <- steps_per_window[1]/14 # 2 days is /14
window_by <- pmax(1, steps_per_day/(4)) # the denominator is number of window starts per day; if trying to increment window by less than the resolution of the time series, just increment by 1 #c(48, 4, 2, 1)

# detrending
mp <- 1 # max polynomial; linear detrending or none at all
mf <- 2 # max fourier
mi <- 1 # max 'interaction' (polynomial in the amplitude of the fourier series)


# ================
# = Set RNG Seed =
# ================
set.seed(42)


# ==============================================
# = Set Up Options for Simulation and Analysis =
# ==============================================
agg_steps_sim <- c(1, 4, 24, 48) #c(1, 12, 288, 288*2) # step sizes for aggregation
steps_per_day_sim <- 24/(agg_steps_sim) #60*24/(5 * agg_steps) # obs per day = (60 min / 1 hr) * (24 hrs / 1 day) * (1 obs / 5*n min)

steps_per_window_sim <- steps_per_day_sim*win_days # steps per window = (n steps / 1 day) * (n days / 1 window)
window_by_sim <- 1 #pmax(1, steps_per_day/(6)) #pmax(1, steps_per_day/(4)) # the denominator is number of window starts per day; if trying to increment window by less than the resolution of the time series, just increment by 1 #c(48, 4, 2, 1)


# =================================
# = Options for Tiered Simulation =
# =================================
nTiers <- 50
nYears_tiered <- 200


# =================================
# = Options for Linear Simulation =
# =================================
nYears_lin <- 200


# ===================================================
# = Stability Properties Guiding Simulation Options =
# ===================================================
critVals <- sort(findCrit()) # critical values of I
Ivals4 <- sort(c(outer(critVals, c(-0.1, 0.1), FUN="-"))) # 4 values of I to use in examples
Irange <- range(Ivals4) # min and max values of I to use in simulation
Ivals <-  seq(Irange[1], Irange[2], length.out=nTiers)



#'   
#' \FloatBarrier  
#'   
#' ***  
#'   


#' #Simulation
#' ##Mathematical Model
NULL
#' \begin{array}{ll}
#' \frac{dX}{dt} = I - X(s+h) + rMR   &(1) \\[0.5em]
#' \frac{dM}{dt} = sX - bM - rMR   &(2) \\[0.5em]
#' R \equiv \frac{X^q}{m^q + X^q}   &(3) \\[2.0em]
#' \text{solving } \frac{dX}{dt}=0 \text{ for } M: \\[0.5em]
#' M = \frac{-I + X(s+h)}{rR}   &(4)\\[0.5em]
#' \text{substituting for } M \text{ in the } dM/dt \text{ equation, I get}\\[0.5em]
#' dM/dt = sX - \frac{-I + X(s+h)}{rR} (b - rR)   &(5) \\[2.0em]
#' \text{solving } \frac{dM}{dt}=0 \text{ for } M: \\[0.5em]
#' M =  \frac{sX}{b+rR}   &(6)\\[0.5em]
#' \text{substituting for } M \text{ in the } dX/dt \text{ equation, I get}\\[0.5em]
#' dX/dt =  I - X(s+h) + rR \frac{sX}{b+rR}   & (7)\\[2.0em]
#' \text{differentiating Eq. 7 with respect to } X\\[0.5em]
#' f''(X) = \frac{rsX^q(b(m^q(q+1)+X^q)+rX^q)}{(bm^q+bX^q+rX^q)^2} -s -h   & (8)\\
#' \end{array}
#'   
#' It is Eqs 7 and 8 above that I will be working from for the rest of this report. Eq 7 gives me the derivative of X as a function of X (and parameters).  When this equation equal zero, there is a equilibrium. When Eq 8 is 0 AND Eq 7 itself is 0, that's a critical point.  
#'   
#'   
#' ##Simulation Helper Functions
#+ helperFunctions
# =====================================
# = Statistics and Plotting Functions =
# =====================================
#      make measured values of class "ts" with frequency = 288 samples per day ----
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


#' ##Figure 1: Phase Portrait
#+ figure1-phasePortrait, fig.width=6, fig.height=6, fig.cap="**Figure 1.** A phase portrait of the system for varying values of P input when q=8. The vector field (indicating the direction and speed that the system moves through phase space at that point) is represented by gray arrows. Nullclines are represented red and blue lines, indicating where dX/dt and dM/dt are equal to zero, respectively.  Trajectories starting at arbitrary initial points (open diamonds) and continuing the along the accompanying solid black line indicate how the system moves from the initial point through phase space for 50 years. Equilibria are indicated by points: solid filled circle is a stable node, an 'X' is a saddle point. An equilibrium occurs whereever the nullclines cross. The different panels correspond to different values of P loading (I). "
par(mfrow=c(2,2), mar=c(2,2,1,0.5), mgp=c(1,0.25,0), tcl=-0.15, cex=1, ps=9, cex.axis=0.85)
# Is <- c(0.75, 1, 1.25, 1.5)
# Is0 <- diff(sort(critVals))*0.5*c(-1,1)+sort(critVals) # a good range of critical values
# Is <- round(c(Is0[1], critVals, Is0[2])-0.1,2)
for(i in 1:length(Ivals4)){
	timeScales::phasePortrait(pars=c(I=Ivals4[i]), nFlow=10, addLeg=(i==2))
	mtext(paste0("I = ",round(Ivals4[i],2)), side=3, line=0, adj=0, font=2)
}
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   


#' ##Time Series for 'Tiered' P Input Increase
#+ tieredSimulation
agg_sos3 <- function(aggsteps){
	out <- tieredTSm[,j={agg_ts(y=value, x=timeVec, width=aggsteps)},by=c("lake","variable","I")]
	out
}
tieredSim <- tieredTSWrapper(nTiers=nTiers, nYears=nYears_tiered, iranges=Ivals, dt=1/24, sin_amp=c(0,0))

# #' ##Figure: Tiered Time Series
# #+ figure2-tiered-time-series, fig.width=3.5, fig.height=5, fig.cap="**Figure.** Time series of simulated lake phosphorus. Top panel shows the time series of phosphorus inputs to the lake. Bottom panels shows the time series of phosphorus concentration in the water over time. Dashed vertical lines are critical values of phosphorus. At the first critical value new equilibria emerge, but original equilibrium point still exists. At the second critical value, the original equilibrium value is annihilated along with a saddle point as the two collide. For each tier of P input, model state variables were initiated at equilibrium and simulated stochastically for 100 years."
# critTimes_cvInd <- apply(outer(Ivals,critVals, FUN='-'), 2, function(x)which.min(abs(x)))
# critTimes_ind <- tieredSim[,I%in%Ivals[critTimes_cvInd]]
# critTimes <- tieredSim[critTimes_ind, timeVec[1], by="I"][,V1]
# par(mfrow=c(2,1), mar=c(2,2,0.5,0.5), mgp=c(1.1,0.25,0), tcl=-0.15, cex=1, ps=9)
# tieredSim[,plot(timeVec, I, type='l', xlab='Year', ylab='P Input (I)')]
# grid()
# abline(v=critTimes, lty=2)
# tieredSim[,plot(timeVec, X, type='l', xlab="Year", ylab="Water P (X)")]
# grid()
# abline(v=critTimes, lty=2)


#' ##Time Series for 'Continuous' P Input Increase
#+ simTimeSeries
# ========================
# = Simulate Time Series =
# ========================
# out <- simP_ts(nYears=100, I_range=c(1.1, 1.4), agg_steps=agg_steps, steps_per_day=steps_per_day, dt=1/steps_per_day[1])
contSim <- simP_ts(nYears=nYears_lin, I_range=Irange, dt=1/steps_per_day_sim[1])
# for(n in names(contSim)){ # make all the elements of the output list objects in the current env
# 	assign(n, contSim[[n]])
# }
stateMat <- contSim$stateMat
dt <- contSim$dt

# ====================
# = Reshape Data Set =
# ====================
# Record like a data set
lakeP <- data.table(lake="Sim", doy=time(ts(stateMat[,"X"], freq=1/dt)), I=stateMat[,"I"], X=stateMat[,"X"], M=stateMat[,"M"])
lakeP_df <- data.frame(lakeP[,list(I,X,M)]) # for finding eigenvalues
lakePm <- melt(lakeP, id.vars=c("lake","doy"))#[variable%in%vars & lake%in%lakes]
# save(lakePm, file="~/Documents/School&Work/epaPost/timeScales/data/lakePm.RData", compress="xz")

# ====================================================
# = Prepare Simulation Data for Statistical Analysis =
# ====================================================
# lakePm[, value:=set_ts(y=log(value), x=doy[1], freq=1/dt), by=c("lake","variable")]
lakePm[, value:=set_ts(y=(value), x=doy[1], freq=1/dt), by=c("lake","variable")]

# aggregate at different time scales ----
lakeP_agg <- lapply(agg_steps_sim, FUN=agg_sos2)
names(lakeP_agg) <- paste0("agg", agg_steps_sim)

# =============================================
# = Xvalues to be used in all TVARSS plotting =
# =============================================
# not to be used until much later in report
# tvarss_xvals <- lapply(lakeP_agg, function(x){x[variable=="X", as.numeric(x)]})


# #' ##Figure: Continous Time Series
# #+ figure-contTimeSeries
# cont_cv_ind <- lakePm[variable=="I",(apply(outer(value,critVals, FUN='-'),2,function(x)which.min(abs(x))))] #apply(outer(Ivals,critVals, FUN='-'), 2, function(x)which.min(abs(x)))
# cont_critTimes <- lakePm[cont_cv_ind,as.numeric(doy)]
# par(mfrow=c(2,1), mar=c(1.85,1.85,0.5,0.5), mgp=c(1,0.25,0), tcl=-0.25, ps=8, cex=1)
# lakePm[variable=="I", plot(ts(value, freq=1/dt), xlab='Year', ylab="P input (I)")]
# abline(v=cont_critTimes, lty=2)
# lakePm[variable=="X", plot(ts(value, freq=1/dt), xlab='Year', ylab="Water P (X)")]
# abline(v=cont_critTimes, lty=2)


#' ##Figure 2: Tiered and Continuous Time Series, with ACF
#+ figure2-tieredContinuous-TS-ACF, fig.width=3, fig.height=5, fig.cap="**Figure 2.** Time series of simulated lake phosphorus. Simulations included 24 observations per year. The first row shows time series of phosphorus input (I), the second row water phosphorus (X), and the third row the autocorrelation function. The left column of panels contains results from a simulation in which I was increased in a stepwise fashion; the model was initiated at equilibrium at the start of each tier of P input, each of which lasted 100 years. The maximum number of lags in the bottom-left panel is equal to the number of observations in 1 tier of P input (100*24). The right column of panels contain results from a simulation in which P input is increased smoothly and linearly over time; each time step has slightly higher I than the previous. The number of lags in the bottom-right panel are equal to the number of observations in 28 years. The vertical dashed lines indicate the times at which alternate equilibria emerge (a stable node and a saddle point), and when the original equilibrium (a stable node) is annihilated as it collides with a saddle point (see Figure 1). Both simulations are stochastic, but the simulation on the right includes an additional sinusoidal process. The simulation with smooth and linear increases in I (right) does not show an abrupt increase in X because the system is lagging behind its equilibrium value, which changes with I."
critTimes_cvInd <- apply(outer(Ivals,critVals, FUN='-'), 2, function(x)which.min(abs(x)))
critTimes_ind <- tieredSim[,I%in%Ivals[critTimes_cvInd]]
critTimes <- tieredSim[critTimes_ind, timeVec[1], by="I"][,V1]
cont_cv_ind <- lakePm[variable=="I",(apply(outer(value,critVals, FUN='-'),2,function(x)which.min(abs(x))))] #apply(outer(Ivals,critVals, FUN='-'), 2, function(x)which.min(abs(x)))
cont_critTimes <- lakePm[cont_cv_ind,as.numeric(doy)]


par(mfcol=c(3,2), mar=c(1.85,1.5,0.5,0.25), oma=c(0.1,0.35,0.3,0.1), mgp=c(0.85,0.15,0), tcl=-0.15, cex=1, ps=9, ylbias=0.3, cex.axis=0.85)

tieredSim[,plot(timeVec, I, type='l', xlab='Year', ylab='')]
abline(v=critTimes, lty=2)
mtext("P Input (I)", side=2, line=1)
mtext("Tiered", side=3, line=-0.2, font=2)
tieredSim[,plot(timeVec, X, type='l', xlab="Year", ylab="")]
mtext("Water P (X)", side=2, line=1)
abline(v=critTimes, lty=2)
acf(tieredSim[,X], lag.max=nYears_tiered/dt, ylab="")
mtext("ACF", side=2, line=1)

lakePm[variable=="I", plot(ts(value, freq=1/dt), xlab='Year', ylab="")]
abline(v=cont_critTimes, lty=2)
mtext("Linear", side=3, line=-0.2, font=2)
lakePm[variable=="X", plot(ts(value, freq=1/dt), xlab='Year', ylab="")]
abline(v=cont_critTimes, lty=2)
acf(lakePm[variable=="X",as.numeric(value)], lag.max=win_days/dt, ylab="")

#'   
#' \FloatBarrier  
#'   
#' ***  
#'   

#' ##Statistics for Tiered Simulations
#+ tieredStatistics
tieredTSm <- melt(tieredSim, measure.vars="X") # undo with: dcast(tieredTSm, ...~variable)
tieredTSm[, value:=set_ts(y=(value), x=timeVec[1], freq=1/dt), by=c("lake","variable","I")]
# tieredTS2[,plot(timeVec, X, type='l')]
tiered_agg <- lapply(agg_steps_sim, agg_sos3)
names(tiered_agg) <- paste0('agg',agg_steps_sim)

tiered_agg_stat <- vector('list', length(tiered_agg))
names(tiered_agg_stat) <- names(tiered_agg)
for(a in 1:length(tiered_agg)){
	ta <- tiered_agg[[a]]
	tiered_agg_stat[[a]] <- ta[,j={
		ar1coef <- ac1(y, trend.rm=FALSE)
		arpcoefs <- arP(y, oType="coefs")
		nP <- length(arpcoefs)
		arpeigs <- arEigs(arpcoefs)
		data.table(x=mean(x), AR1=ar1coef, nP=nP, ARpEig=max(Mod(arpeigs)), sdVal=sd(y, na.rm=TRUE))
		
	}, by=c("lake","variable","I")]
}


#' ##Statistics for Continous Simulations
#+ continousStatistics
simAC <- roll_ac.sos(lakeP_agg, window_elapsed=steps_per_window_sim, vars="X", lakes="Sim", DETREND=TRUE, by=window_by_sim, save_output=F, mp=mp, mf=mf, mi=mi)
simARp <- roll_ac.sos(lakeP_agg, fit_arP=TRUE, window_elapsed=steps_per_window_sim, vars="X", lakes="Sim", DETREND=TRUE, by=window_by_sim, save_output=F, mp=mp, mf=mf, mi=mi)

rollSD_contSim <- function(d, w){
	d[variable%in%"X" & lake%in%"Sim"][,j={roll_ts(y=y, x=x, FUN=sd, width=w, DETREND=TRUE, save_output=FALSE, mp=mp, mf=mf, mi=mi)}, by=c("lake","variable")]
}
simSD <- mapply(rollSD_contSim, d=lakeP_agg, w=steps_per_window_sim, SIMPLIFY=FALSE)


#' ##Figure 3: Statistics for Simulations
#+ figure3-simulationStatistics, fig.width=3, fig.height=5, fig.cap="**Figure 3.** Autocorrelation calculated at different time scales for the two simulated scenarios of tiered increases in P input (left column) and continuous linear increases in P input (right column). For the tiered increases, autocorrelation was calculated within a given tier, never across tiers. For the continuous linear increases, autocorrelation was calculated in a backwards-looking rolling window. Each row of panels corresponds to a different level of pre-analysis aggregation: in the top row no aggregation was performed, in the second row groups of 4 observations were averaged, in the third row the 24 annual observations were averaged, and for the fourth row 2 'years' of observations were averaged. Vertical dashed lines correspond to the critical points described in Figs1&2."
critTimes_cvInd <- apply(outer(Ivals,critVals, FUN='-'), 2, function(x)which.min(abs(x)))
critTimes_ind <- tieredSim[,I%in%Ivals[critTimes_cvInd]]
critTimes <- tieredSim[critTimes_ind, timeVec[1], by="I"][,V1]
cont_cv_ind <- lakePm[variable=="I",(apply(outer(value,critVals, FUN='-'),2,function(x)which.min(abs(x))))] #apply(outer(Ivals,critVals, FUN='-'), 2, function(x)which.min(abs(x)))
cont_critTimes <- lakePm[cont_cv_ind,as.numeric(doy)]

par(mfcol=c(4,2), mar=c(1.5,1.8,0.25,0.25), oma=c(0.5,0.8,0.5,0.1), mgp=c(1,0.25,0), tcl=-0.15, cex=1, ps=9)

for(i in 1:length(tiered_agg_stat)){
	ta <- tiered_agg_stat[[i]]
	ty <- paste0(1/dt/as.numeric(gsub("agg","",names(tiered_agg_stat)[i])), " obs per yr")
	ta[,plot(x, AR1, type='l', ylab=ty, xlab="")]
	if(i==1){
		mtext("Tiered", side=3, line=-0.2, xpd=TRUE, font=2)
	}
	abline(v=critTimes, lty=2)
}
mtext("AR(1) Coef", side=2, outer=TRUE, line=-0.1, font=2)
mtext("Year", side=1, line=1)


for(i in 1:length(simAC)){
	ta <- simAC[[i]]
	# ty <- paste0(1/dt/as.numeric(gsub("agg","",names(simAC)[i])), " obs per yr")
	ta[,plot(x, y, type='l', ylab="", xlab="")]
	if(i==1){
		mtext("Linear", side=3, line=-0.2, xpd=TRUE, font=2)
	}
	abline(v=cont_critTimes, lty=2)
}
mtext("Year", side=1, line=1)


# plot standard deviation for the 'linear' simulation
# par(mfrow=c(4,1), mar=c(2,2,0.5,0.5), mgp=c(1,0.25,0), tcl=-0.15, cex=1, ps=9)
# for(i in 1:length(simSD)){
# 	simSD[[i]][,plot(x,y, type='l')]
# }


#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Field Experiment
#' ##Data Functions
#+ data-prep-aggregation
agg_sos <- function(aggsteps){
	out <- sosm[,j={agg_ts(y=value, x=doy, width=aggsteps)},by=c("lake","variable")]
	out
}

#' ##Plotting Functions
#' ###Plotting ACF
#+ functions-plottingACF
plot_acf <- function(ln=c("Paul","Peter"), v=c("chla", "bga", "logChla", "logBga"), na.action=na.exclude, lag.max=288*12, ...){
	ln <- match.arg(ln)
	v <- match.arg(v)
	dots <- list(...)
	if(is.null(dots$main)){main <- paste(ln, v, 'acf')}else{main <- dots$main}
	d <- sos[lake==ln, get(v)]
	acf2 <- function(..., main=main){acf(...)}
	o <- acf2(d, lag.max=lag.max, na.action=na.action, ...)
	mtext(main, side=3, line=0.1, font=2)
	invisible(NULL)
}

#' ###Plotting ACF Heat Maps
#+ acf-map-functions
acf_map <- function(out, ...){
	obs_lab <- attr(out, "xlab")
	lag_lab <- attr(out, "ylab")
	rwbCols <- colorRampPalette(c("blue","white","red"))(256) 
	image(x=obs_lab, y=lag_lab[-1], z=out[,2:ncol(out)], col=rwbCols, ...)
	invisible(NULL)
}

add_legend <- function(out, legend.mar=2, col, axis.args){
	if(missing(col)){
		col <- colorRampPalette(c("blue","white","red"))(256)
	}
	if(missing(axis.args)){
		axis.args <- list(mgp=c(0.5, 0.15, 0), tcl=-0.1)
	}
	fields::image.plot(out, legend.only=TRUE, col=col, legend.mar=legend.mar, axis.args=axis.args)
}

add_axis <- function(out, side=1){
	map_colors <- colorRampPalette(c("blue","white","red"))(256)
	ylab_pretty <- pretty(attr(out, "ylab")/288)
	ylab_pretty[1] <- 1/288 #attr(out_L_sub, "ylab")[3]/288
	ylab_names <- sapply(ylab_pretty, interval_name, minPerSample=1440) #interval_name(attr(out_L_sub, "ylab"))
	axis(side=side, at=ylab_pretty*288, labels=ylab_names)
	invisible()
}

sub_out <- function(out, ind=list(1,1), type=c("sub", "thin")){
	type <- match.arg(type)
	xlab <- attr(out, "xlab")
	ylab <- attr(out, "ylab")
	nr <- nrow(out)
	nc <- ncol(out)
	if(type=="thin"){
		rVec <- ((1:nr)%%ind[[1]])==0
		cVec <- ((1:nc)%%ind[[2]])==0
	}
	if(type=="sub"){
		rVec <- ind[[1]]
		cVec <- ind[[2]]
	}
	out2 <- out[rVec, cVec]
	attr(out2, "xlab") <- xlab[rVec]
	attr(out2, "ylab") <- ylab[cVec]
	return(out2)
}

#' ###Plotting Peter Paul AC, Diff
#+ function-plotLakeAC
plotac <- function(X, ...){
	X <- copy(X)
	ylim <- X[,range(y, na.rm=TRUE)]
	
	ydiff <- X[lake=="Peter", y] - X[lake=="Paul", y]
	zdiff <- data.table(lake="zdiff", variable=X[,variable[1]], x=X[lake==lake[1], x], y=ydiff)
	X <- rbind(X, zdiff)
	
	ul <- X[,unique(lake)]
	for(l in 1:length(ul)){
		dud <- X[lake==ul[l],j={
			tcol <- c("Paul"="blue","Peter"="red", "zdiff"="black")[lake[1]]
			if(lake[1]=="Paul"){
				p1 <- function(..., ylab, xlab, ylab2="", xlab2=""){
					plot(..., ylab=ylab2, xlab=xlab2)
					mtext(ylab, side=2, line=par("mgp")[1])
					mtext(xlab, side=1, line=par("mgp")[1])
				}
				p1(x,y, type='l', col=tcol, ylim=ylim, ...)
			}else if(lake[1]=="Peter"){
				lines(x, y, col=tcol)
			} else if(lake[1]=="zdiff"){
				p2 <- function(..., ylab, ylab2="", xlab, xlab2=""){
					plot(..., ylab=ylab2, xlab=xlab2)
					mtext(xlab, side=1, line=par("mgp")[1])
				}
				p2(x,y, type='l', col=tcol, ...)
			}
		
			NULL
		}]
	}
	
	invisible()
}

#' ##Data Prep
#+ data-prep-basic
#      drop tuesday lake ----
sos <- sos_data[Lake!="Tuesday" & Year==2015]

#      shorter names ----
setnames(sos, 
	old=c("Year", "Lake", "DoY", "Temp_HYLB", "Chl_HYLB", "Chl_logged_HYLB", "BGA_HYLB","BGA_logged_HYLB"),
	new=c("year","lake","doy","wtr","chla", "logChla", "bga", "logBga")
)

#      ensure numeric, re-structure data set ----
sos[,bga:=as.numeric(bga)]
sos[,logBga:=as.numeric(logBga)]
sos[,chla:=as.numeric(chla)]
sos[,logChla:=as.numeric(logChla)]
sosm <- melt(sos, id.vars=c("year","lake","doy"))[variable%in%vars & lake%in%lakes]
sosm[,variable:=as.character(variable)]


sosm[, value:=set_ts(y=(value), x=doy[1]), by=c("year","lake","variable")]

#      grab range limits (primarily for plotting) ----
doy_range <- sos[,range(doy, na.rm=TRUE)]
var_range <- sos[,range(get(vars[1]), na.rm=TRUE)]

# data aggregation
sos_agg <- lapply(agg_steps, agg_sos)
names(sos_agg) <- paste0("agg", agg_steps)


#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##Figure 4: Time Series & ACF
#+ figure4-fieldTimeSeries-ACF, fig.width=3.5, fig.height=3.5, fig.cap="**Figure 4.** Time series and autocorrelation function of high frequency chlorophyll fluorescence measurements in Peter Lake (manipulated) and Paul Lake (reference) in 2015. The red vertical dashed line indicates the day when fertilization was halted in Peter Lake.", results='hide'
par(mfrow=c(2,2), mar=c(1.75, 1.0, 0.25, 0.25), oma=c(0,0.75,0.5,0), mgp=c(1, 0.2, 0), tcl=-0.15, ps=9, cex=1, cex.axis=0.85, ylbias=0.3)
sosm[lake=="Paul" & variable==vars, plot(doy, (value), xlim=doy_range, col="black", type='l', xlab="", ylab="")]
mtext(sosm[,field_ylab[variable[1]]], side=2, line=0.85)
mtext("Day of year", side=1, line=0.85)

mtext("Reference", side=3, line=-0.1, font=2)
sosm[lake=="Peter" & variable==vars, plot(doy, (value), xlim=doy_range, col="black", type='l', xlab="", ylab="")]
abline(v=180, col='red', lty=2)
mtext("Day of year", side=1, line=0.85)
mtext("Manipulated", side=3, line=-0.1, font=2)

plot_acf(v=vars, ylab="", main="", xlab="")
mtext("ACF", side=2, line=0.85)
mtext("Lag", side=1, line=0.85)
plot_acf(ln='Peter', v=vars, ylab="", main="", xlab="")
mtext("Lag", side=1, line=0.85)

#'   
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##Figure 5: Rolling Window Autocorrelation
#+ rollingWindowAC-calculation, cache=FALSE
AC_list <- roll_ac.sos(sos_agg, window_elapsed=steps_per_window, vars=vars, lakes=lakes, DETREND=TRUE, by=window_by, mp=mp, mf=mf, mi=mi)

#+ figure5-rollingWindowAC-PaulPeterDifference, fig.width=3, fig.height=5, fig.cap="**Figure 5.** Rolling windows of first-order autocorrelation from detrended chlorophyll time series. Blue lines are from Paul Lake (reference), red lines are Peter Lake (manipulated). In the second column, the black lines represent the difference (Peter - Paul) between the lines in the first column (positive values indicate that autocorrelation was higher in Peter than in Paul). Time series series in each row were aggregated to a different time scale prior to analysis."
# ylabs <- paste0(c("chla"="Chl-a", "bga"="Phyco", "logChla"="log(Chl-a)", "logBga"="log(Phyco)")[vars], "\nAR(1) (", sapply(agg_steps, interval_name), ")")
ylabs <- paste0("\nAR(1) (", sapply(agg_steps, interval_name), ")")

xlabs <- rep("", length(agg_steps))
xlabs[length(agg_steps)] <- "Day of Year"

par(mfrow=c(length(agg_steps),2), mar=c(1.25,1.25,0.25,0.25), oma=c(0.75,0.75,0.25,0.1), cex=1, tcl=-0.15, mgp=c(1,0.2,0), ps=9, ylbias=0.3, cex.axis=0.85)
invisible(mapply(plotac, X=AC_list, ylab=ylabs, xlab=xlabs))
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##Calculate ACF Map
#+ acf-map-calculate, cache=FALSE
out_L <- acf_roll(x=sosm[lake=="Paul" & variable==vars, value], width=steps_per_window[1], by=window_by[1], lag.max=acf_lag.max, DETREND=TRUE, mp=mp, mf=mf, mi=mi)
out_R <- acf_roll(x=sosm[lake=="Peter" & variable==vars, value], width=steps_per_window[1], by=window_by[1], lag.max=acf_lag.max, DETREND=TRUE, mp=mp, mf=mf, mi=mi)

#' ##Figure 6: ACF Heat Map w/ Time Series Insets
#+ figure6-acf-map-full-tsInsets-figure, fig.width=6, fig.height=7, fig.cap="**Figure 6.** Autocorrelation across many time scales, using the ACF function. Each window is detrended first. Time series in the insets represent subsets of the full heat map at specific time scales."
#      Setup Layout Matrix ----
nMain <- 3 # the number of heat map panels
ts_choices <- c(1, 12*6, 12*24, 12*48)
nScales <- length(ts_choices)
widthFac <- 3 # how much wider the heat maps are relative to the time series
pExpand <- 6 # number of times (*) to expand each time series panel
bExpand <- 1 # number of buffer panels to add above and below each set of nScales time series panels
panelVec <- rep(rep((1:nMain),each=nScales*pExpand+2*bExpand),widthFac)
tsBase <- seq(4, by=nMain, length.out=nScales)+rep(0:2, each=nScales)
tsPos <- rep(tsBase, each=pExpand)
buffer_mat <- matrix(rep(0, nMain*bExpand),ncol=nMain)
ts_buffer_mat <- rbind(buffer_mat, matrix(tsPos, ncol=nMain), buffer_mat)
lay_mat <- matrix(c(c(ts_buffer_mat), panelVec), ncol=widthFac+1)

#      Add Legend -- Custom Function ----
add_legend2 <- function(inputDat){
	mapLegend(x=1.03, y=0.5, w=0.05, h=0.9, zlim=range(inputDat), cols=c("blue","white","red"), horiz=FALSE, axSide=4, lab.cex=par("cex.axis"), lab.sig=2, offset=0.2, xpd=TRUE)
}
add_panel_lab_main <- function(let){mtext(let, side=3, line=-0.85, adj=0.01, font=2, cex=1.2)}

#      Plot Heat Maps ----
layout(lay_mat)
par(mar=c(1.5,2,1,3), oma=c(0.5,0,0,0), mgp=c(1,0.2,0), tcl=-0.15, ps=9, cex=1, las=0, cex.axis=0.85)
acf_map(out_L, xlab="", ylab="Time scale", main="Paul Lake (reference)", yaxt='n')
add_axis(out_L)
# add_legend(out_L)
zrange <- range(c(range(out_L), range(out_R), range(out_R-out_L)))
add_legend2(out_L)
add_panel_lab_main("A")

par(cex=1)
acf_map(out_R, xlab="", ylab="Time scale", main="Peter Lake (manipulated)", yaxt='n')
add_axis(out_R)
# add_legend(out_R)
add_legend2(out_R)
add_panel_lab_main(LETTERS[1+(nScales+1)])

par(cex=1)
out_Diff <- out_R - out_L
acf_map(out_Diff, xlab="", ylab="Time scale", main="Difference (manipulated - reference)", yaxt='n', xpd=TRUE)
mtext("Day of year", side=1, line=1, xpd=TRUE)
add_axis(out_Diff)
# add_legend(out_Diff)
add_legend2(out_R - out_L)
add_panel_lab_main(LETTERS[1+(nScales+1)*2])

#      Plot Time Series ----
suppressWarnings({par(mar=c(0.25,1.25,0.1,0.1), mgp=c(1,0.15,0), tcl=0.15, las=1, ps=9)})
xval <- attr(out_L, "xlab")
for(s in nScales:1){ # iterate through time scales more slowly than throw main plots (paul, peter, diff)
	ts_ind <- list(r=1:nrow(out_L), c=ts_choices[s]+1)
	for(i in 1:nMain){ # iterate through main plots more quickly, plotting same time scale for each of paul, peter, diff
		tout <- switch(i, out_L, out_R, (out_R-out_L))
		ts_temp <- c(sub_out(tout, ind=ts_ind, type="sub"))
		ylim <- range(ts_temp)+c(0, diff(range(ts_temp))*0.2) # add a bit of extra white space at the top of the insets to make room for labels
		plot(xval, ts_temp, xlab="", ylab="", type='l', xaxt='n', yaxt='n', ylim=ylim)
		axis(side=2, labels=TRUE, at=pretty(ts_temp, n=3), xpd=F)
		suppressWarnings({axis(side=1, labels=(s==1), mgp=c(1,-0.2, 0))})
		mtext(interval_name(ts_choices[s]), side=3, adj=0.98, font=2, line=-0.75)
		revS <- (nScales:1)[s]
		panelLab <- LETTERS[(i-1)*nScales+revS+i] # I appologize to my future self if he needs to understand this. Remember that the plots are created in a very jumbled way due to 1) the layout(), and 2) I go nScales:1 not 1:nScales, and 3) the s loop is outside the i loop [this point interacts with #1, such that I think they cancel each other]. Also, I had to fiddle a bit so I'm not even sure I understand the pattern, so don't be confused by those two points, especially the first, it might not be relevant. The +i at the end is just b/c I want the heat maps to be labeled A, E, and I, such that all the Paul Lake panels can be summarized as A-D, all the Peter Lake panels as E-H, and all the 'difference' panels as I-L.
		mtext(panelLab, side=3, adj=0.04, font=2, line=if(panelLab=='B'){-0.75}else{-0.75}, cex=1)
		# mtext(paste(panelLab, interval_name(ts_choices[s]), sep=", "), side=3, adj=0.99, font=2, line=-0.67)
		# A
	}
}



#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Supplement or Extra Figures
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##Figure: Correlation of ACF across Time Scales
#+ crossScaleCorrelationACF, fig.width=3, fig.height=6, fig.cap="**Figure.** Color is the correlation between time scales."
add_legend3 <- function(inputDat){
	mapLegend(x=1.03, y=0.5, w=0.05, h=0.9, zlim=range(inputDat), cols=fields::tim.colors(), horiz=FALSE, axSide=4, lab.cex=par("cex.axis"), lab.sig=2, offset=0.2, xpd=TRUE)
}

library(fields)
corL <- cor(out_L[,-1])
corR <- cor(out_R[,-1])
corDiff <- cor((out_R-out_L)[,-1])
y <- attributes(out_L)$ylab

par(mfrow=c(3,1), mar=c(2,2,0.5,3), mgp=c(1,0.2,0), tcl=-0.15, ps=9, cex=1, cex.axis=0.85)
image(x=y[-1], y=y[-1], z=corL[-1,-1], col=tim.colors(), xlab="", ylab="", xaxt='n', yaxt='n')
add_axis(out_L, side=1)
add_axis(out_L, side=2)
add_legend3(corL)

image(x=y[-1], y=y[-1], z=corR[-1,-1], col=tim.colors(), xlab="", ylab="", xaxt='n', yaxt='n')
add_axis(out_R, side=1)
add_axis(out_R, side=2)
add_legend3(corR)

image(x=y[-1], y=y[-1], z=corDiff[-1,-1], col=tim.colors(), xlab="", ylab="", xaxt='n', yaxt='n')
add_axis(out_R, side=1)
add_axis(out_R, side=2)
add_legend3(corDiff)

#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##Figure: Significance of Correlation of ACF across Time Scales
#+ signif-crossScaleCorrelationACF, fig.width=3, fig.height=5, fig.cap="**Figure.** Color is the correlation between time scales."
n <- nrow(out_L)
tstat  <- function(r, n){r*sqrt((n-2)/(1-r^2))}
tp <- function(tstat,df){pt(tstat, df=n-2, lower.tail=FALSE)}
ttest_pval <- function(corMat, n, adjMethod="BH"){
	out <- tp(tstat(corMat, n), n)#[-1,-1]
	out[] <- p.adjust(out, method='BH')
	return(out)
}

ttest_pval_L <- ttest_pval(corL, n) #tp(tstat(corL, n), n)#[-1,-1]
ttest_pval_R <- ttest_pval(corR, n)
ttest_pval_Diff <- ttest_pval(corDiff, n)

par(mfrow=c(3,1), mar=c(2,2,0.5,0.5), mgp=c(1,0.2,0), tcl=-0.15, ps=9, cex=1, cex.axis=0.85)

image(x=y[-1], y=y[-1], z=(ttest_pval_L<0.05), xlab="", ylab="", xaxt='n', yaxt='n')
add_axis(out_L, side=1)
add_axis(out_L, side=2)

image(x=y[-1], y=y[-1], z=(ttest_pval_R<0.05), xlab="", ylab="", xaxt='n', yaxt='n')
add_axis(out_L, side=1)
add_axis(out_L, side=2)

image(x=y[-1], y=y[-1], z=(ttest_pval_Diff<0.05), xlab="", ylab="", xaxt='n', yaxt='n')
add_axis(out_L, side=1)
add_axis(out_L, side=2)



#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##Figure: Full ACF Heat Map
#+ figureS-acf-map-full-figure, fig.width=3.5, fig.height=6, fig.cap="**Figure** Autocorrelation across many time scales, using the ACF function. Each window is detrended first.", fig.show='hold', include=TRUE
#      Thin-out for Fast/ Lighter Plotting ----
out_L_sub <- sub_out(out_L, ind=list(r=8, c=4), type='thin')
out_R_sub <- sub_out(out_R, ind=list(r=8, c=4), type='thin')
out_Diff_sub <- out_R_sub - out_L_sub

#      Begin Plotting ----
xlimL <- c(min(attr(out_L_sub, "xlab")), 240)
xlimR <- c(min(attr(out_R_sub, "xlab")), 240)

par(mfrow=c(3,1))
par(mar=c(2,2,1,3), mgp=c(1,0.2,0), tcl=-0.15, ps=8, cex=1)
acf_map(out_L_sub, xlab="", ylab="Time scale", main="Paul Lake (reference)", xlim=xlimL, yaxt='n')
add_axis(out_L_sub)
add_legend(out_L_sub)

par(cex=1)
acf_map(out_R_sub, xlab="", ylab="Time scale", main="Peter Lake (manipulated)", xlim=xlimR, yaxt='n')
add_axis(out_R_sub)
add_legend(out_R_sub)

par(cex=1)
acf_map(out_Diff_sub, xlab="Day of year", ylab="Time scale", main="Difference", xlim=xlimR, yaxt='n')
add_axis(out_Diff_sub)
add_legend(out_Diff_sub)
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   

#' ##Figure: Subset ACF Heat Map
#+ figureS-acf-map-subset-figure, fig.width=3.5, fig.height=6, fig.cap="**Figure** Autocorrelation at a across many time scales, using the ACF function. Each window is detrended first. Subset of full data set (zoom on high frequencies and early part of the time series)."
#      Subset to Zoom in on Relevant Bits ----
rInd <- attr(out_L, "xlab") <= 190
cInd <- attr(out_L, "ylab") <= 144
out_L_sub2 <- sub_out(out_L, ind=list(r=rInd, c=cInd), type='sub')
out_R_sub2 <- sub_out(out_R, ind=list(r=rInd, c=cInd), type='sub')
out_Diff_sub2 <- out_R_sub2 - out_L_sub2

#      Begin Plotting ----
par(mfrow=c(3,1))
par(mar=c(2,2,1,3), mgp=c(1,0.2,0), tcl=-0.15, ps=8, cex=1)
acf_map(out_L_sub2, xlab="", ylab="Time scale", main="Paul Lake (reference)", yaxt='n', zlim=range(out_L_sub2))
add_axis(out_L_sub2)
add_legend(out_L_sub2)

par(cex=1)
acf_map(out_R_sub2, xlab="", ylab="Time scale", main="Peter Lake (manipulated)", yaxt='n', zlim=range(out_R_sub2))
add_axis(out_R_sub2)
add_legend(out_R_sub2)

par(cex=1)
acf_map(out_Diff_sub2, xlab="Day of year", ylab="Time scale", main="Difference", yaxt='n', zlim=range(out_Diff_sub2))
add_axis(out_Diff_sub2)
add_legend(out_Diff_sub2)
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#+ figureS-acf-map-full-tsInsets-figure, fig.width=6, fig.height=6, fig.cap="**Figure** Autocorrelation at a across many time scales, using the ACF function. Each window is detrended first. Time series in the insets represent subsets of the full heat map at specific time scales. Each time scale is standardized, with its mean removed and divded by its standard deviation."
#      Plot Heat Maps ----
layout(lay_mat)
par(mar=c(1.5,2,1,3), oma=c(0.5,0,0,0), mgp=c(1,0.2,0), tcl=-0.15, ps=8, cex=1, las=0)
acf_map(scale(out_L), xlab="", ylab="Time scale", main="Paul Lake (reference)", yaxt='n')
add_axis(scale(out_L))
# add_legend(out_L)
zrange <- range(c(range(scale(out_L)), range(scale(out_R)), range(scale(out_R-out_L))))
add_legend2(scale(out_L))
add_panel_lab_main("A")


par(cex=1)
acf_map(scale(out_R), xlab="", ylab="Time scale", main="Peter Lake (manipulated)", yaxt='n')
add_axis(scale(out_R))
# add_legend(out_R)
add_legend2(scale(out_R))
add_panel_lab_main(LETTERS[1+(nScales+1)])

par(cex=1)
out_Diff <- scale(out_R - out_L)
acf_map(out_Diff, xlab="", ylab="Time scale", main="Difference (manipulated - reference)", yaxt='n', xpd=TRUE)
mtext("Day of year", side=1, line=1, xpd=TRUE)
add_axis(out_Diff)
# add_legend(out_Diff)
add_legend2(scale(out_R - out_L))
add_panel_lab_main(LETTERS[1+(nScales+1)*2])

#      Plot Time Series ----
suppressWarnings({par(mar=c(0.25,1.25,0.1,0.1), mgp=c(1,0.15,0), tcl=0.15, las=1, ps=8)})
xval <- attr(out_L, "xlab")
for(s in nScales:1){ # iterate through time scales more slowly than throw main plots (paul, peter, diff)
	ts_ind <- list(r=1:nrow(out_L), c=ts_choices[s]+1)
	for(i in 1:nMain){ # iterate through main plots more quickly, plotting same time scale for each of paul, peter, diff
		tout <- switch(i, scale(out_L), scale(out_R), scale(out_R-out_L))
		ts_temp <- c(sub_out(tout, ind=ts_ind, type="sub"))
		plot(xval, ts_temp, xlab="", ylab="", type='l', xaxt='n', yaxt='n')
		axis(side=2, labels=TRUE, at=pretty(ts_temp, n=3), xpd=F)
		suppressWarnings({axis(side=1, labels=(s==1), mgp=c(1,-0.2, 0))})
		mtext(interval_name(ts_choices[s]), side=3, adj=0.98, font=2, line=-0.75)
		revS <- (nScales:1)[s]
		panelLab <- LETTERS[(i-1)*nScales+revS+i] # I appologize to my future self if he needs to understand this. Remember that the plots are created in a very jumbled way due to 1) the layout(), and 2) I go nScales:1 not 1:nScales, and 3) the s loop is outside the i loop [this point interacts with #1, such that I think they cancel each other]. Also, I had to fiddle a bit so I'm not even sure I understand the pattern, so don't be confused by those two points, especially the first, it might not be relevant. The +i at the end is just b/c I want the heat maps to be labeled A, E, and I, such that all the Paul Lake panels can be summarized as A-D, all the Peter Lake panels as E-H, and all the 'difference' panels as I-L.
		mtext(panelLab, side=3, adj=0.04, font=2, line=if(panelLab=='B'){-0.75}else{-0.75}, cex=1)
		# mtext(paste(panelLab, interval_name(ts_choices[s]), sep=", "), side=3, adj=0.99, font=2, line=-0.67)
		# A
	}
}

#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' # Info
#+ Info, results='markup'
difftime(Sys.time(), t1) # how long it took to run these models/ produce this report
Sys.time()
sessionInfo()







#' ---
#' title: "Applying early warning statistics to simulated lake P: a test of detrending, DLMs, and AR(p) models"
#' author: "Ryan Batt"
#' date: "2018-04-01"
#' abstract: |
#'        I had previously analyzed field data for early warning signals. Particularly, I was testing if the early warning signals varied across the time scale of the analysis/ data collection. My analysis indicated that the strength of these signals, which took the form of increasing autocorrelation, varied across scales. However, questions remained whether the detrending method I used muted the signals at some frequencies, and whether other statistics might be more robust to time scale. These other statistics include the leading eigenvalue of an AR(p) model, which might be robust to time scale because higher order AR terms reach further back in time, and so might help to bridge the gap across time scales. AR(p) models can be fit in several contexts, one of which is a dynamic linear model (DLM). DLM's could also be convenient for circumventing the issue of rolling window size. Here I simulate a critical transition in lake P, and use the simulated time series to test the effects of time scales and detrending on the performance of DLM and rolling window models, as well as on AR(1) vs AR(p) models.
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
# 	"~/Documents/School&Work/epaPost/timeScales/pkgBuild/manuscript/timeScales_dlm_arp_sim.R",
# 	output_format=o_f,
# 	output_dir='~/Documents/School&Work/epaPost/timeScales/pkgBuild/manuscript',
# 	clean = TRUE
# )

Sys.setenv(PATH=paste0("/Library/TeX/texbin:",Sys.getenv("PATH")))
opts_chunk$set(
	fig.path = 'timeScales_dlm_arp_sim/', 
	cache.path='timeScales_dlm_arp_sim/',
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
agg_steps <- c(1, 4, 24, 48) #c(1, 12, 288, 288*2) # step sizes for aggregation
steps_per_day <- 24/(agg_steps) #60*24/(5 * agg_steps) # obs per day = (60 min / 1 hr) * (24 hrs / 1 day) * (1 obs / 5*n min)

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

#' #Simulate Time Series
#+ simTimeSeries
# ========================
# = Simulate Time Series =
# ========================
# out <- simP_ts(nYears=100, I_range=c(1.1, 1.4), agg_steps=agg_steps, steps_per_day=steps_per_day, dt=1/steps_per_day[1])
out <- simP_ts(nYears=100, I_range=c(1.1, 1.4), dt=1/steps_per_day[1])
for(n in names(out)){ # make all the elements of the output list objects in the current env
	assign(n, out[[n]])
}


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
lakeP_agg <- lapply(agg_steps, FUN=agg_sos2)
names(lakeP_agg) <- paste0("agg", agg_steps)


# =============================================
# = Xvalues to be used in all TVARSS plotting =
# =============================================
# not to be used until much later in report
tvarss_xvals <- lapply(lakeP_agg, function(x){x[variable=="X", as.numeric(x)]})



#' #Section 1: Time Series Simulation and Equilibria
#' ##Calculate Eigenvalues and Roots
#+ eigsRoots, results='markup'
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

#' ##Figure 1a: Time Series
#+ figure1a-timeSeries, fig.width=3.5, fig.height=4.5, fig.cap="**Figure 1a.** Time series of water P and sediment P.  P Input to the lake (not shown) is increasing linearly over time. The tipping point is reach somewhere around time step 80, and the system shifts shortly thereafter. The simulated lines do not represent equilibria, but are a stochastic realization. There are cycles of period = 1 unit 'Time' of very small amplitude (small relative to the size of shift at the end of the time series).", results='hide'
# ===================
# = Plot Simulation =
# ===================
# dev.new(width=3.5, height=4.5)
par(mfrow=c(2,1), mar=c(1.85,1.85,0.5,0.5), mgp=c(1,0.25,0), tcl=-0.25, ps=8, cex=1)
lakePm[variable=="X", plot(ts(value, freq=1/dt), ylab="X, water P")]
lakePm[variable=="M", plot(ts(value, freq=1/dt), ylab="M, sediment P")]


#' ##Figure 1b: Time Series w/ Roots
#+ figure1b-timeSeriesRoots, fig.width=3.5, fig.height=4.5, fig.cap="**Figure 1b.** Time series as shown previously, but now with the 'roots' calculated. For each value in the black lines, the corresponding equilibrium value was calculated with the black line's value as a starting point. This root is plotted in forest green (intended to be a dashed line, btu the plotting messes up, sorry). It might be important to note that the simulated state is often pretty far off from the equilibrium value. This behavior persists even when the parameter is changed more slowly (say, over 3,000 years instead of just 100 years) or w/ a small delta-t. I see no important consequences of this yet, but might be worth noting. ", results='hide'
# ===============================
# = Plot Roots and Observations =
# ===============================
# hmm, the system is very far from equilibrium
# it's especially strange how the simulation P way overshoots the equilibrium value
# in other simulations, i don't recall
# dev.new(width=3.5, height=4.5)
# png("~/Desktop/simP_simValues_trueEquilibria.png", width=3.5, height=3.5, res=150, units='in')
par(mfrow=c(2,1), mar=c(2,2,0.5,0.5), ps=8, mgp=c(1,0.25, 0), tcl=-0.15, cex=1)
lakeP[,plot(as.numeric(doy), X, type='l', ylim=range(c(X, root.X), na.rm=TRUE), xlab="Day", ylab="Water P")]
lakeP[,lines(as.numeric(doy), root.X, lty=2, col='forestgreen')]
lakeP[,plot(as.numeric(doy), M, type='l', ylim=range(c(M, root.M), na.rm=TRUE), xlab="Day", ylab="Mud P")]
lakeP[,lines(as.numeric(doy), root.M, lty=2, col='forestgreen')]
# dev.off()

#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##Figure 2: Eigenvalues at Simulated Values and Roots
#+ figure2-eigenvalSim, fig.width=3.5, fig.height=5.5, fig.cap="**Figure 2.** Eigenvalues of the simulation (indicating stability) across time (top horizontal axis) and across nutrient loading (bottom horizontal axis; in the simulation, nutrient loading increased linearly with time). Eigenvalues were calculated from the Jacobian matrix (square matrix of partial derivatives of state variables w/ respect to each other). The Jacobian was computed numerically by calculating the rate of change of the system around a mulutivariate system state (mud P and water P). I calculated the Jacobian, and therefore the eigenvalues, using 2 sets of system state: 1) the system state in the simulated time series (top panel); 2) the system state as the equilibrium value given starting points (and parameters) in the simulated time series. In other words, the top panel corresponds to the eigenvalues in the black lines in **Figure 1b**, and the bottom panel corresponds to the dark green lines in **Figure 1b**. Plotted eigenvalues are the maximum of the modulus of the Jacobian's eigenvalues. The bottom panel shows the canonical exponential increase in eigenvalue as the tipping point is approached (bifurcation at I= ~1.3356, which is where the ||eigenvalue|| = 0)."
# ===============================
# = Plot Simulation Eigenvalues =
# ===============================
# dev.new(width=3.5, height=3.5)
par(mfrow=c(2,1), mar=c(2,2,2,0.5), ps=8, mgp=c(1,0.25, 0), tcl=-0.15)
plot(simP_eigs$I, Mod(simP_eigs$X1), type='l', xlab="Nutrient Loading", ylab="Eigenvalues around observed states")
par(new=TRUE)
plot(ts(Mod(simP_eigs$X1), freq=1/dt), type='n', xlab="", ylab="", xaxt='n', yaxt='n')
axis(side=3)
mtext("Time", side=3, line=1)


# =================================================
# = Plot 'True' Eigenvalues of Root of Simulation =
# =================================================
# png("~/Desktop/simP_eigenvalues.png", width=3.5, height=3.5, res=150, units='in')
# dev.new(width=3.5, height=3.5)
# par(mar=c(2,2,2,0.5), ps=8, mgp=c(1,0.25, 0), tcl=-0.15)
plot(rootEigs$I, Re(rootEigs$X1), xlab="Nutrient Loading", ylab="Eigenvalue around Root", type='l') # plot of eigenvalues across values of I (loading)
par(new=TRUE)
plot(ts(Mod(rootEigs$X1), freq=1/dt), type='n', xlab="", ylab="", xaxt='n', yaxt='n')
axis(side=3)
mtext("Time", side=3, line=1.0)
# dev.off()


#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Section 2: Rolling Window AR(1) and AR(p)
#' ##Calculate Rolling Window AR(1)
#+ rollingWindowAR1

# ====================================================================
# = Rolling Window AR(1); 4 Time Scales, Detrended and No Detrending =
# ====================================================================
# Calculate AR(1) for both detrending and not ----
simAC <- roll_ac.sos(lakeP_agg, window_elapsed=steps_per_window, vars="X", lakes="Sim", DETREND=TRUE, by=window_by, save_output=F)
simAC_noDetrend <- roll_ac.sos(lakeP_agg, window_elapsed=steps_per_window, vars="X", lakes="Sim", DETREND=FALSE, by=window_by)

#' ##Figure 3: Rolling Window AR(1)
#+ figure3-rollingWindowAR1, fig.width=4, fig.height=6, fig.cap="**Figure 3.** Time series of rolling window first-order autocorrelation (AR(1)). Windows are 28 'years' long, with 24 time steps per year (dt=1/24). The left column of panels are calculated by detrending the time series in each 28-day window; detrending has a polynomial component that has up to a linear and quadratic term (up to order 2), a Fourier series component with up to 2 sin-cos pairs of period=24, and a possible interaction of the Fourier pairs and the polynomial term(s). The order of detrending was selected by AICc. After detrending, the AR(1) coefficient is calculated as the correlation coefficient between x_t and x_{t-1}. In the right column of panels, detrending is not used before calculating AR(1). The rows differ by the degree to which the time series were **agg**regated: the top panel had no aggregation, the second panel was aggregated by taking the average of every 2 conseuctive observations (not overlapping, such that when agg=2 dt effectively becomes 1/24/2=1/12), the third row was aggregated taking the average of each group of 24 observations (annual average), and the fourth row was aggregated by taking the average of 48 consecutive observations (yielding a biennial time series). **Note:** This convention of detrending and aggregation, and this arrangement of panels, will be re-used in upcoming plots."
# png("~/Desktop/simP_rollingAC_Detrend.png", width=3.5, height=5.5, res=150, units='in')
ylabs <- paste0("AR(1) Coeff (", gsub("agg", "agg=", names(lakeP_agg)), ")")
# dev.new(width=4, height=6)
par(mfcol=c(4,2), mar=c(1,2,0.5,0.5), ps=8, mgp=c(1,0.25, 0), tcl=-0.15, oma=c(1,0,1.0,0), cex=1)
mapply(plotac_simp, simAC, ylab=ylabs, xlab="")
mtext("Time", side=1, line=1)
mtext("Detrended", side=3, outer=TRUE, adj=0.25, line=0)
mapply(plotac_simp, simAC_noDetrend, ylab=ylabs, xlab="")
mtext("No Detrending", side=3, outer=TRUE, adj=0.85, line=0)
mtext("Time", side=1, line=1)
# dev.off()


#' ##Calculate Rolling Window AR(p)
#+ rollingwindowARp
# ===============================
# = AR(p) Model: Rolling Window =
# ===============================
simARp <- roll_ac.sos(lakeP_agg, fit_arP=TRUE, window_elapsed=steps_per_window, vars="X", lakes="Sim", DETREND=TRUE, by=window_by, save_output=F)
simARp_nod <- roll_ac.sos(lakeP_agg, fit_arP=TRUE, window_elapsed=steps_per_window, vars="X", lakes="Sim", DETREND=FALSE, by=window_by, save_output=F)


#' ##Figure 4: Rolling Window AR(p)
#+ figure4-rollingWindowARp, fig.width=4, fig.height=6, fig.cap="**Figure 4.** The modulus of the leading eigenvalue of AR(p) models fit to 28-year rolling window time series. AR(p) models were fit using R's ar() function, with the fitting method set to argument corresponding to 'Yule-Walker'; the maximum order considered for model selection was set to the length of the windowed time series (28 years times 24 observations per year) divided by 10; i.e., a very high limit. Selected orders were small (e.g., less than 10). The order was selected for each window separately, and so cannot be displayed for the time series as a whole."
ylabs <- paste0("||Eig|| of AR(p) (", gsub("agg", "agg=", names(lakeP_agg)), ")")
# dev.new(width=4, height=6)
par(mfcol=c(4,2), mar=c(1,2,0.5,0.5), ps=8, mgp=c(1,0.25, 0), tcl=-0.15, oma=c(1,0,1.0,0), cex=1)
mapply(plotac_simp, simARp, ylab=ylabs, xlab="")
mtext("Time", side=1, line=1)
mtext("Detrended", side=3, outer=TRUE, adj=0.25, line=0)
mapply(plotac_simp, simARp_nod, ylab=ylabs, xlab="")
mtext("Time", side=1, line=1)
mtext("No Detrending", side=3, outer=TRUE, adj=0.75, line=0)

#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Section 3: Use Dynamic Linear Models to Fit Time-Varying AR(1) and AR(p)
#' The general form of the models is as follows:  
#' $$\begin{array}{ll}
#' y_t = z_t + v_t     &v\sim \mathcal{N}(0, 1/\tau_v) \\
#' z_t = C_t + \Phi_t\bold{B}^p(Z-C) + \omega_t     &\omega \sim \mathcal{N}(0,1/\tau_{\omega}) \\
#' C_t=C_{t-1}+\xi_t & \xi \sim \mathcal{N}(0,1/\tau_\xi)\\
#' \Phi_t = \Phi_{t-1} + \varepsilon_t &\varepsilon \sim \mathcal{N}([0,\dots,0], \Sigma_\Phi)\\[0.5em]
#' &\Sigma_\Phi = \begin{bmatrix} \sigma^2_{\varphi_1}& \dots& 0\\ \vdots& \ddots& \vdots\\ 0& \dots& \sigma^2_{\varphi_p} \end{bmatrix}\\[2em]
#' \text{if }  p=2: \begin{cases}
#' \Phi_t=\begin{bmatrix} \varphi_{1,t}, \varphi_{2,t} \end{bmatrix}\\[0.5em]
#' \bold{B}^p(Z)=\begin{bmatrix}Z_{t-1}\\Z_{t-2} \end{bmatrix}\\[1em]
#' \bold{B}^p(C)=\begin{bmatrix}C_{t-1}\\C_{t-2} \end{bmatrix}\\
#' \end{cases}\\
#' \end{array}$$
#'   
#' Note that I will explore the application of this model with and without the following elements: 1) detrending (for whole time series at once); 2) time-varying mean (the parameter C above); 3) an AR(1) vs an AR(p) model. When exploring AR(p) models, the **order of the model was selected** using R's `ar()` function (as above). A more appropriate approach, albeit computationally intensive, would be to fit many time-varying models of different orders and compare their deviance information criteria (DIC). However, for now I'll be sticking to determining model order ahead of time b/c it is computationally cheap.  
#' 
#' ##Fit TVAR(1)SS Models
#+ fitTVAR1SSmodels
# ===================
# = Fit TVARSS DLMs =
# ===================
# TVAR(1)SS
tvarss_noMean_list <- tvarss_wrapper(lakeP_agg, det=TRUE) # constant mean, detrended
tvarss_noMean_nod_list <- tvarss_wrapper(lakeP_agg) # constant mean, no detrending

tvarss_list <- tvarss_wrapper(lakeP_agg, tvMean=TRUE, det=TRUE) # varying mean, detrended
tvarss_nod_list <- tvarss_wrapper(lakeP_agg, tvMean=TRUE) # varying mean, no detrending

#' ##Figure 5: DLM TVAR(1)SS w/ Constant Mean
#+ prepare-xylabs
# prepare xylabs
# Plot TVAR(1)SS ----
ylabs <- paste0("TVAR(1) of DLM (", gsub("agg", "agg=", names(lakeP_agg)), ")")
xlabs <- c(rep("", (length(lakeP_agg)-1)), "Time")

#+ figure5-TVAR1SS-noMean, fig.width=4, fig.height=6, fig.cap="**Figure 5.** Time series of the relative density of the posterior of the first-order autocorrelation coefficient (AR(1)) as fit by a Bayesian state space model (dynamic linear model) in which the mean parameter is held constant. Results are shown for the case where the entire time series was first detrended, and the case where no detrending was performed. The color shown is the relative desnity of the posterior; each posterior (i.e., at each time step) was relativized by dividing the density at all values (all values of the coefficient) by the maximum density. Therefore, at each time step, the (y-axis) value at which the maximum posterior density for that time step occurred will have a dark red color, the minimum value will have a blue color. "
# plot version with a 'constant' mean parameter
# dev.new(width=4, height=6)
par(mfcol=c(4,2), mar=c(1,2,0.5,0.5), ps=8, mgp=c(1,0.25, 0), tcl=-0.15, oma=c(1,0,1.5,0), cex=1)
ppout <- mapply(plotPost.tvarss, x=tvarss_noMean_list, xvals=tvarss_xvals, ylab=ylabs, xlab=xlabs, MoreArgs=list(main="", varName="Phi"))
mtext("Constant Mean, Detrended", side=3, line=0, outer=TRUE, adj=0.15)
ppout <- mapply(plotPost.tvarss, x=tvarss_noMean_nod_list, xvals=tvarss_xvals, ylab=ylabs, xlab=xlabs, MoreArgs=list(main="", varName="Phi"))
mtext("Constant Mean, No Detrending", side=3, line=0, outer=TRUE, adj=0.85)

#' ##Figure 6: DLM TVAR(1)SS w/ Time-Varying Mean
#+ figure6-TVAR1SS, fig.width=4, fig.height=6, fig.cap="**Figure 6.** Time series of posterior samples of first-order autocorrelation (AR(1)) fit by a DLM. Like Figure 5, but here the mean parameter varies over time. Note that a time-varying mean parameter should effectively detrend the time series. Colors are as in Figure 5."
# plot version with a time-varying mean parameter
# dev.new(width=4, height=6)
par(mfcol=c(4,2), mar=c(1,2,0.5,0.5), ps=8, mgp=c(1,0.25, 0), tcl=-0.15, oma=c(1,0,1.5,0), cex=1)
ppout <- mapply(plotPost.tvarss, x=tvarss_list, xvals=tvarss_xvals, ylab=ylabs, xlab=xlabs, MoreArgs=list(main="", varName="Phi"))
mtext("Varying Mean, Detrended", side=3, line=0, outer=TRUE, adj=0.15)
ppout <- mapply(plotPost.tvarss, x=tvarss_nod_list, xvals=tvarss_xvals, ylab=ylabs, xlab=xlabs, MoreArgs=list(main="", varName="Phi"))
mtext("Varying Mean, No Detrending", side=3, line=0, outer=TRUE, adj=0.85)

#' ##Fit TVAR(p)SS Models
#+ fitTVARpSSmodels
# TVAR(p)SS ----
tvarss_noMean_nP_list <- tvarss_wrapper(lakeP_agg, fit_arP=TRUE, det=TRUE)
tvarss_noMean_nP_nod_list <- tvarss_wrapper(lakeP_agg, fit_arP=TRUE)

tvarss_nP_list <- tvarss_wrapper(lakeP_agg, tvMean=TRUE, fit_arP=TRUE, det=TRUE)
tvarss_nP_nod_list <- tvarss_wrapper(lakeP_agg, tvMean=TRUE, fit_arP=TRUE)

# Plot TVAR(p)SS ----
# part of y-labels for tvar(p)ss
# different than for tvar(1)ss b/c can include the model order (p)
# also note that for the rolling window ar(p) this could not be done b/c
# p can differ among windows; in DLM, p is constant
# Note that these aren't the full y-labels, just the order (p)
ylabs_nM_nP <- sapply(tvarss_noMean_nP_list, function(x)dim(x$Phi)[3])
ylabs_nM_nP_nod <- sapply(tvarss_noMean_nP_nod_list, function(x)dim(x$Phi)[3])
ylabs_nP <- sapply(tvarss_nP_list, function(x)dim(x$Phi)[3])
ylabs_nP_nod <- sapply(tvarss_nP_nod_list, function(x)dim(x$Phi)[3])

#' ##Figure 7: DLM TVAR(p)SS w/ Constant Mean
#+ figure7-TVARpSS, fig.width=4, fig.height=6, fig.cap="**Figure 7.** Time series of posterior samples of the modulus of the leading eigenvalue of the AR coefficients, of which there are p. The mean parameter was held constant over time. Similar to the interpretation of an AR(1) coefficient, the system loses stability as this value approaches 1 from below. The left column is detrended, right column was not detrended. Colors are as in Figure 5."
# plot version with a 'constant' mean parameter
# dev.new(width=4, height=6)
par(mfcol=c(4,2), mar=c(1,2,0.5,0.5), ps=8, mgp=c(1,0.15, 0), tcl=-0.15, oma=c(1,0,1.5,0), cex=1)
ylabs <- paste0("Eig of TVAR(p=", ylabs_nM_nP, ");", gsub("agg", "agg=", names(lakeP_agg)))
ppout <- mapply(plotPost.tvarss, x=tvarss_noMean_nP_list, xvals=tvarss_xvals, ylab=ylabs, xlab=xlabs, MoreArgs=list(main="", varName="Eigen"))
mtext("Constant Mean, Detrended", side=3, line=0, outer=TRUE, adj=0.17)
# mtext("Time", side=1, line=0.75)
ylabs <- paste0("Eig of TVAR(p=", ylabs_nM_nP_nod, ");", gsub("agg", "agg=", names(lakeP_agg)))
ppout <- mapply(plotPost.tvarss, x=tvarss_noMean_nP_nod_list, xvals=tvarss_xvals, ylab=ylabs, xlab=xlabs, MoreArgs=list(main="", varName="Eigen"))
mtext("Constant Mean, No Detrending", side=3, line=0, outer=TRUE, adj=0.95)
# mtext("Time", side=1, line=0.75)



#' ##Figure 8: DLM TVAR(p)SS w/ Time-Varying Mean
#+ figure8-TVARpSS, fig.width=4, fig.height=6, fig.cap="**Figure 8.** Time series of posterior samples of the modulus of the leading eigenvalue of the AR coefficients, of which there are p. The mean parameter varied over time. Colors are as in Figure 5."
# plot version with a time-varying mean parameter
# dev.new(width=4, height=6)
par(mfcol=c(4,2), mar=c(1,2,0.5,0.5), ps=8, mgp=c(1,0.15, 0), tcl=-0.15, oma=c(1,0,1.5,0), cex=1)
ylabs <- paste0("Eig of TVAR(p=", ylabs_nP, ");", gsub("agg", "agg=", names(lakeP_agg)))
ppout <- mapply(plotPost.tvarss, x=tvarss_nP_list, xvals=tvarss_xvals, ylab=ylabs, xlab=xlabs, MoreArgs=list(main="", varName="Eigen"))
mtext("Varying Mean, Detrended", side=3, line=0, outer=TRUE, adj=0.15)
# mtext("Time", side=1, line=0.75)
ylabs <- paste0("Eig of TVAR(p=", ylabs_nP_nod, ");", gsub("agg", "agg=", names(lakeP_agg)))
ppout <- mapply(plotPost.tvarss, x=tvarss_nP_nod_list, xvals=tvarss_xvals, ylab=ylabs, xlab=xlabs, MoreArgs=list(main="", varName="Eigen"))
mtext("Varying Mean, No Detrending", side=3, line=0, outer=TRUE, adj=0.95)
# mtext("Time", side=1, line=0.75)

#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Section 4: Time-Varying Regression
#' ##Fit Time-Varying Regression Model
#+ tvRegFit
# =========================================
# = AR(1) Time Varying Regressoin (TVReg) =
# =========================================
dX <- detrendR(ts(stateMat[,"X"], frequency=steps_per_day[1]), max_poly=2, max_fourier=2, returnType="resid")
tvmod <- tvAR(dX, type='none', bw=0.6)
tvmod_nod <- tvAR(ts(stateMat[,"X"], frequency=steps_per_day[1]), bw=0.6)

tvreg_ar1 <- c(NA, tvmod$tvcoef)
tvreg_ar1_nod <- c(NA, tvmod_nod$tvcoef[,1])

arp_mod_full <- ar(dX, order.max=steps_per_day[1]*2)
nP <- length(arp_mod_full$ar)
arp_mod_full_nod <- ar(ts(stateMat[,"X"], frequency=steps_per_day[1]), order.max=steps_per_day[1]*2)
nP_nod <- length(arp_mod_full$ar)

tvmodp <- tvAR(dX, type='none', bw=0.1, p=nP)
tvmodp_nod <- tvAR(stateMat[,"X"], bw=0.1, p=nP_nod)

alleigs <- function(x){
	max(Mod(arEigs(x)))
}
tvmodp_eigs <- apply(tvmodp$tvcoef, 1, alleigs)
tvmodp_nod_eigs <- apply(tvmodp_nod$tvcoef, 1, alleigs)


#' ##Figure 9: Time Varying Regression
#+ figure9-tvReg, fig.height=5, fig.width=5, fig.cap="**Figure 9.** Time series of results from time-varying regression. Left column was detrended, right column was not detrended. First row are time series of AR(1) coefficient (time varying), and second row are time series of the modulus of leading eigenvalue of the AR(p) coefficients. I allowed the mean parameter to vary when the time series was not first detrended, and held the mean constant when the time series was detrended prior to model fitting. This model was only applied to the highest frequency time series (the top panel in previous figures). The method works by fitting a AR coefficients at entire time series at each time step using weighted least squares; the weights are determined from a kernel, which has the effect of down-weighting the influence of observations far from the focal data point (either far in the past or far in the future). [Note: The kernel used here is called 'Epanesnikov' in the help file, but I believe the spelling is 'Epanechnikov'; its shape is similar to cosine (see Wikipedia: https://en.wikipedia.org/wiki/Kernel_(statistics)#Kernel_functions_in_common_use).] Therefore, as in a rolling window analysis, a separate model is fit at each time step, and data points close to that time step influence the statistic. In this method, the influential data points can be in the past or future, and the rate of decay of their influence as distance from the focal point increases is determined by a bandwidth parameter. In this way, the bandwidth parameter serves a similar role as the window size. The bandwidth parameter can be selected by cross validation, but in experimenting, I found that the selected bandwidth was often quite large, resulting in minimal temporal change in the autocorrelation coefficient (and mean) reported. The method is easy to compute (fast), and can produce plausible results consistent with expectations; however, the problems of selecting window size are not really avoided as a bandwidth still needs to be determined (and cross validation, as implemented here, was not always reliable). Furthermore, the kernel approach, as implemented here w/ symmetric kernels, does not make much sense in the applied scenario of not having future data; however, a different shaped kernel might work better to avoid this problem. It should also be noted in the DLM models above (fit in JAGS), all data points are used to compute the statistic at each time step (not an on-line approach), though could still be calculated sequentially (using full available time series up to that date) in the applied scenario."
# dev.new(width=5, height=5)
par(mfrow=c(2,2), mar=c(2,2,0.5,0.5), ps=8, mgp=c(1,0.25, 0), tcl=-0.15, oma=c(0,0,1.5,0), cex=1)
plot(ts(tvreg_ar1, freq=1/dt), xlab="Time", ylab="AR(1) Coeff (time-varying reg)")
mtext("Detrended", side=3, line=-0.5, outer=TRUE, adj=0.25)
plot(ts(tvreg_ar1_nod, freq=1/dt), xlab="Time", ylab="AR(1) Coeff (time-varying reg)")
mtext("No Detrending", side=3, line=-0.5, outer=TRUE, adj=0.75)
plot(ts(tvmodp_eigs, frequency=steps_per_day[1]), type='l', xlab="Time", ylab=paste0("Eig of AR(p=",nP,") (time-varying reg)"))
plot(ts(tvmodp_nod_eigs, frequency=steps_per_day[1]), type='l', xlab="Time", ylab=paste0("Eig of AR(p=",nP,") (time-varying reg)"))


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
#' ---
#' title: "The time scale of resilience loss: the effect of sampling frequency on an early warning statistic"
#' author: "Ryan Batt"
#' date: "2017-08-21"
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
# =================
# = Load Packages =
# =================
library("data.table")
library("zoo")
library("forecast")

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
# 	"~/Documents/School&Work/epaPost/timeScales/pkgBuild/manuscript/timeScales_report.R",
# 	output_format=o_f,
# 	output_dir='~/Documents/School&Work/epaPost/timeScales/pkgBuild/manuscript',
# 	clean = TRUE
# )

Sys.setenv(PATH=paste0("/Library/TeX/texbin:",Sys.getenv("PATH")))
opts_chunk$set(
	fig.path = 'manuscript_report/', 
	cache.path='manuscript_report/',
	echo=TRUE, 
	include=TRUE, 
	cache=F,
	results='asis',
	warning=FALSE,
	fig.show="hold",
	fig.lp = if(o_f=="html_document"){"**Figure.**"}else{NULL}
)

# setwd("~/Documents/School&Work/epaPost/timeScales/pkgBuild/manuscript")
# source("../manuscript/manuscript_figures_functions.R")

#' #Options
#+ options
win_days <- 28 # window size in days covered
agg_steps <- c(1, 12, 144, 288) # step sizes for aggregation
lakes <- "Peter" # can be vector; lakes to analyze (Paul, Peter)
vars <- "chla" # can be vector; variables to analyze (wtr, bga, chla)



#' #Functions in Developement
#' ##Plotting
#+ functions-plotting
plot_acf <- function(ln=c("Paul","Peter"), v=c("chla", "bga"), na.action=na.exclude, lag.max=288*12, ...){
	ln <- match.arg(ln)
	v <- match.arg(v)
	dots <- list(...)
	if(is.null(dots$main)){main <- paste(ln, v, 'acf')}else{main <- dots$main}
	d <- sos[lake==ln, get(v)]
	o <- acf(d, lag.max=lag.max, main="", na.action=na.action)
	mtext(main, side=3, line=0.1, font=2)
	invisible(NULL)
}


#' ##Statistics
#+ functions-statistics
# Detrend
# Detrend time series
# @param x vector time series
# @param time vector of time steps
#
# @return detrended x
# @export
detrend <- function(x, time=1:length(x)){
	stats::residuals(stats::lm(x~time))
}

# AC1
# Calculate first-order autocorrelation
#
# @param x numeric vector of values of a regular (equally-spaced) time series
#
# @return AR(1) coefficient
# @export
ac1 <- function(x, ...){	
	# # option 1
	l2 <- embed(x, 2)
	ac <- cor(l2[,1], l2[,2], use="na.or.complete")
	#
	# # option 2
	# ac <- ar(x, order.max=1)$ar # returns numeric(0) if nothing fit
	#
	# # option 3
	# l2 <- embed(x, 2)
	# ac <- lm(I(l2[,1]) ~ I(l2[,2]))$coef[2]
	
	# option 4
	# detX <- detrend(x)
	# if(all(is.na(x))){return(NA)}
#	# l2 <- stats::embed(x, 2)
#	# ac <- tryCatch(stats::lm(I(l2[,1]) ~ I(l2[,2]))$coef[2], error=function(cond)NA)
	
	# # option 4
	# l2 <- embed(x, 2)
	# out <- RollingWindow::RollingCorr(l2[,2], l2[,1], window=win)
	
	return(ac)
}

ac_sub <- function(x, n, phase, ...){
	l2 <- embed(x, 2)
	row_vec <- seq_len(nrow(l2))
	row_ind <- sub_samp(row_vec, n=n, phase=phase)
	l2_sub <- l2[row_ind, ]
	ac <- cor(l2[,1], l2[,2], use="na.or.complete")
	return(ac)
}


#' ##Rolling Windows and Aggregation
#+ functions-windows-agg
roundGrid <- function(x, frac=1){
	# if frac is 1, then place in a 1º grid
	# if frac is 0.5, then place in the 0.5º grid
	floor(round(x/frac, 6))*frac+frac/2
}

roll_ts <- function(y, width=288, by=1, FUN=mean, x, ...){
	buff <- rep(NA, width-1)
	if(by > 1){
		mat <- sub_embed(y, width=width, n=by) 
	}else{
		mat <- embed(y, width)
	}
	agg <- c(buff, apply(X=mat, MARGIN=1, FUN=FUN, ...))
	if(!missing(x)){
		if(by > 1){
			mat2 <- sub_embed(x, width=width, n=by) 
		}else{
			mat2 <- embed(x, width)
		}
		agg2 <- c(buff, apply(mat2, 1, max, na.rm=TRUE))
		data.table(x=agg2, y=agg)
	}else{
		agg
	}
}
# @examples
# # example data
# x <- 1:50
# y <- cumsum(rnorm(50))
# dt <- data.table(x, y)
#
# # do x and y values separately
# # x value is an ordering value, like time
# # y is the system state, like chlorophyll
# xa <- roll_ts(dt[,x], FUN=max, width=5)
# ya <- roll_ts(dt[,y], width=5)
#
# # do both x and y at same time
# dta <- dt[,agg_ts(y=y, x=x, width=5)]
#
# # plot two sets of results
# par(mfrow=c(2,1))
# plot(xa, ya)
# dta[,plot(x,y)]

agg_ts <- function(x, y, width=288, na.rm=TRUE, FUN=mean){
	buff <- rep(NA, width-1)
	tot <- round(mean(1/diff(x), na.rm=TRUE), 0)
	frac <- width/tot
	dt <- data.table(x=roundGrid(x, frac), y=y)
	dto <- dt[,list(y=FUN(y, na.rm=na.rm), N=sum(!is.na(y))), by=x]
}

sub_samp <- function(x, n, phase=1){
	if(phase!=n){
		phase <- max(1, phase%%n) # ensures 1 ≤ phase ≤ n; keeps starting elements when phase > n, unlike min(phase,n)
	}
	vec <- rep(FALSE, n)
	vec[phase] <- TRUE
	x[vec]
}

#' Embed and Subsample a Time Series
#' 
#' Embed a time series like in \code{\link{embed}}, and then subsample the rows of that resulting matrix.
#' 
#' @param x numeric vector
#' @param width integer size of window/ embedding dimension
#' @param n integer, sample every nth element of \code{x}
#' @param phase integer, start subsampling sequency on the nth element of \code{x}
#' 
#' @return matrix
sub_embed <- function(x, width=1, n=1, phase=1){
	emat <- embed(x, dimension=width)
	row_vec <- seq_len(nrow(emat))
	row_ind <- sub_samp(row_vec, n=n, phase=phase)
	return(emat[row_ind, ])
}
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   


#' #Data Prep
#+ data-prep
sos <- sos_data[Lake!="Tuesday"] # drop tuesday lake
setnames(sos, 
	old=c("Year", "Lake", "DoY", "DateTime", "Temp_HYLB", "Chla_Conc_HYLB", "BGA_Conc_HYLB"),
	new=c("year","lake","doy","datetime","wtr","chla","bga")
) # shorter names
sos[,bga:=as.numeric(bga)]
sosm <- melt(sos, id.vars=c("year","lake","doy","datetime"))
# just a few handy summaries for use in graphing (x and y limits)
doy_range <- sos[,range(doy, na.rm=TRUE)]
chla_range <- sos[,range(chla, na.rm=TRUE)]
bga_range <- sos[,range(bga, na.rm=TRUE)]
wtr_range <- sos[,range(wtr, na.rm=TRUE)]
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   


#' #Fluoresence Time Series
#' ##Plot Chlorophyll & BGA Time Series for each Lake
#+ chla-bga-timeSeries-figure, fig.width=5, fig.height=5, fig.cap="**Figure.** High frequency chlorophyll (chla, micrograms per liter) and blue-green algae (bga cells per liter) time series in Peter (red) and Paul (blue) lakes in 2015.", results='hide'
par(mfcol=c(2,2), mar=c(2.5, 2.0, 0.25, 0.25), mgp=c(1, 0.25, 0), tcl=-0.15, ps=8, cex=1)
sos[lake=="Paul", plot(doy, chla, xlim=doy_range, col="blue", type='l')]
sos[lake=="Peter", plot(doy, chla, xlim=doy_range, col="red", type='l')]
sos[lake=="Paul", plot(doy, bga, xlim=doy_range, col="blue", type='l')]
sos[lake=="Peter", plot(doy, bga, xlim=doy_range, col="red", type='l')]
#'   
#'   
#' Paul Lake is an unmanipulated, reference lake. Peter was being fertilized with nitrogen and phosphorphus every day. The goal was to create a critical transition in Peter Lake (Hopf birfurcation).  
#'   
#' Around day 180, Peter Lake has a blue-green algal bloom. Maximum BGA in Peter is ~4x's higher than maximum BGA in Paul Lake.
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   


#' #Influence of Sampling Frequency on AC
#' ##ACF Plots for Chla and BGA in each Lake
#+ ACF-figure, fig.width=5, fig.height=5, fig.cap="**Figure.** ACF of chlorophyll a (chla, indicator of algal biomass) and blue-green algae (bga, indicator of blue-green biomass; blue-green algae create HABs).", results='hide'
par(mfrow=c(2,2), mar=c(2.5, 2.0, 1, 0.25), mgp=c(1, 0.25, 0), tcl=-0.15, ps=8, cex=1)
plot_acf()
plot_acf(v="bga")
plot_acf(ln='Peter')
plot_acf(ln='Peter', v='bga')
#'   
#'   
#' The figures show the autocorrelation statistic calculated at different time lags. These time series have a 5-minute data observation frequency. Thus, a lag of 1 is a 5-minute lag. There are 288 samples in a day.  
#'   
#' Correlation is a description of the relatedness between two variables. If you take a single variable, and correlated it with itself, the correlation is 1. If you take a variable, and correlated it with itself *after "lagging" it by a certain number of time steps* (shifting the duplicated series relative to the reference series), you get an "autocorrelation" statistic. The AutoCorrelation Function (ACF) gives the autocorrelation statistic. The ACF is 1 when the lag is 0 because this is equivalent to the correlation between a time series and its unlagged copy. By lagging the copy, you get the average relatedness (correlation) between values of a time series with past values of the time series.  
#'   
#' The ACF plots are cyclical. They reach peaks at multiples of 288, because this is comparing chlorophyll concentrations at the same times of day (on different days). Basically, chlorophyll at 12:00 is more correlated with chlorophyll at 12:00 the day before than it is with the concentration at 11:00 (or 1:00) the day before. The troughs in the plots represent perfectly out-of-phase points in the time series --- every 288 points after the first 144 (comparing "noon" to "midnight", or 2:00pm to 2:00am).  
#'   
#' The point of this figure is to show that the time scale of analysis affects the calculation of the AC statistic. AC is a statistic used to gauge changes in ecosystem resilience, which can then be used as an early warning signal of critical transitions like algal blooms.  
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   


#' #Calculate Rolling Averages of Different Scales (hourly, 12-hr, 24-hr)
#'   
#' In the past, I have simply aggregated high frequency data to a lower frequency by taking the average. For example, with 5-min data there are 288 samples per day. So I would average each day's 288 samples into a single value (mean), and then analyze the resultant "daily" time series.  
#'   
#' Below, I replicate this process (24-hr avg), but also aggregate at two other frequencies (1-hr and 12-hr).  
#'   
#+ rollingAvg-chla-3
agg_sos <- function(aggsteps){
	out <- sosm[,j={agg_ts(y=value, x=doy, width=aggsteps)},by=c("lake","variable")]
	out
}
sos_agg <- lapply(agg_steps, agg_sos)
names(sos_agg) <- paste0("agg", agg_steps)



# sos_12 <- sos[lake=="Peter", agg_ts(y=chla, x=doy, width=12)]
# sos_144 <- sos[lake=="Peter", agg_ts(y=chla, x=doy, width=144)]
# sos_288 <- sos[lake=="Peter", agg_ts(y=chla, x=doy, width=288)]

#+ rollingAvg-chla-3-fig, fig.width=3.5, fig.height=6.5, fig.cap="**Figure.** Rolling averages of chlorophyll in Peter Lake. The window size for the rolling average varies among panels.", results='hide'
interval_name <- function(x){
	# assumes 5 minute data
	interv <- x*5
	if(interv < 60){
		unit <- "min"
		val <- interv
	}else if(interv >= 60 & interv < 1440){
		unit <- "hr"
		val <- round(interv/60, 2)
	}else{
		unit <- "day"
		val <- round(interv/60/24, 2)
	}
	iname <- paste(val, unit, sep="-")
	return(iname)
}

plot_agg_ts <- function(steps, ln="Peter", vn="chla"){
	ns <- length(steps)
	par(mfrow=c(ns,1), mar=c(1.75, 1.75, 1, 0.25), oma=c(1,0.1,0.1,0.1), mgp=c(1, 0.25, 0), tcl=-0.15, ps=8, cex=1)
	for(s in 1:ns){
		iname <- interval_name(steps[s])
		tdat <- sos_agg[[paste0("agg",steps[s])]][lake==ln & variable==vn]
		ylab <- paste0(ln, " ", vn, " (", iname, " avg)")
		plot(tdat[,x], tdat[,y], xlab="", ylab=ylab, type='l')
	}
	mtext("Day of Year", side=1, outer=TRUE, line=0.1)
}

plot_agg_ts(agg_steps)
#
# sos_12[,plot(x, y, xlab="DoY", ylab="Peter Chla (1 hr avg)", type='l')]
# sos_144[,plot(x, y, xlab="DoY", ylab="Peter Chla (12 hr avg)", type='l')]
# sos_288[,plot(x, y, xlab="DoY", ylab="Peter Chla (24 hr avg)", type='l')]
#'   
#'   
#' The longer the window for the rolling average, the smoother the time series. There is a lot of high-frequency variability ("noise"). Performing a non-overlapping rolling mean (non-overlapping because a dataum is never used for more than 1 over the averages/ in more than 1of the windows) causes upward and downward fluctuations to "cancel out", producing a smoother time series. The larger the window, the larger the sample size that goes into each average, and the smoother the result.
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   


#' #Influence of Sampling Frequency on AC as EWS
#' Early warning statistics are usually calculated in a "rolling window" fashoin. Unlike the daily average performed above, rolling windows for the EWS are usually overlapping. The idea is that for each point, a statistic (like autocorrelation=AC) is calculated for the past N observations. The window is shifted forward by 1 observation, and calculated for N observations again. Thus, adjacent windows are contain N-1 of the same data points.
#'   
#' In this case, we run into a bit of a conundrum: should windows contain a fixed number of data points (as described above), or should they cover a fixed amount of time? I.e., if a window size is 28 data points, it would cover 28 days for daily data, or 14 days for 12-hr data, or 1 and 1/6 days for hourly data.  
#'   
#' Below, I take both approaches: calculating the AC statistic for a fixed number of points (28 points), and again for a fixed period of time (28 days). For the daily data, these two approaches will be equivalent. The outcome will likely differ for the hour and 12-hr time series.  
#'   
#' ###Fixed Window Size (number observations)
#+ rollingAC-fixedSize
roll_ac.sos <- function(X, window_elapsed, nby=1, vars, lakes){
	# check vars, set if missing
	if(missing(vars)){
		vars <- X[[1]][,unique(variable)]
	}else{
		vars <- match.arg(vars, choices=X[[1]][,unique(variable)], several.ok=TRUE)
	}
	
	# check lakes, set if missing
	if(missing(lakes)){
		lakes <- X[[1]][,unique(lake)]
	}else{
		lakes <- match.arg(lakes, choices=X[[1]][,unique(lake)], several.ok=TRUE)
	}
	
	# check/ set window elapsed class
	if(!is.list(window_elapsed)){
		window_elapsed <- as.list(window_elapsed)
	}
	
	# helper function to apply rolling statistic
	roll_ac <- function(X2, nsteps){
		X2[variable%in%vars & lake%in%lakes][,j={roll_ts(y=y, x=x, FUN=ac1, width=nsteps, by=nby)}, by=c("lake","variable")]
	}
	out <- mapply(roll_ac, X, window_elapsed, SIMPLIFY=FALSE)
	out
}

sos_agg_ac <- roll_ac.sos(X=sos_agg, window_elapsed=win_days*24*60/5/agg_steps, vars=vars, lakes=lakes)

# sos_ac1_fS_12 <- sos_12[,j={roll_ts(y=y, x=x, FUN=ac1, width=28)}]
# sos_ac1_fS_144 <- sos_144[,j={roll_ts(y=y, x=x, FUN=ac1, width=28)}]
# sos_ac1_fS_288 <- sos_288[,j={roll_ts(y=y, x=x, FUN=ac1, width=28)}]

#+ rollingAC-fixedSize-figure, fig.width=3.5, fig.height=6.5, fig.cap="**Figure.** Rolling window AC(1). Panels differ in the window size of the rolling average (applied before calcluating AC statistic). AC statistic is calculated from the same number of points in each panel.", results='hide'
par(mfrow=c(3,1), mar=c(2.5, 2.0, 1, 0.25), mgp=c(1, 0.25, 0), tcl=-0.15, ps=8, cex=1)
sos_ac1_fS_12[,plot(x, y, xlab="DoY", ylab="AC(1) for Peter Chla (1 hr avg)", type='l')]
sos_ac1_fS_144[,plot(x, y, xlab="DoY", ylab="AC(1) Peter Chla (12 hr avg)", type='l')]
sos_ac1_fS_288[,plot(x, y, xlab="DoY", ylab="AC(1) Peter Chla (24 hr avg)", type='l')]
#'   
#'   
#' There is minimal signal in AR(1) for the hourly data, and a very clear increase in AR(1) for the daily data; 12-hr data is somewhere in between. Each point in all three panels is autocorrelation calculated from the same number of data points, but due to differences in sampling frequency, the amount of time covered by those data points differs.  
#'   
#' Using this approach, there is a strong influence of sampling frequency on not only the AR(1) value (note y-axis scaling), but also the change in the statisic over time. Namely, the bottom panel shows a clear increase leading up to day ~195, whereas the top panel shows no trend at all.  
#'   

#' ###Fixed Window Time (amount of time elapsed)
#'   
#'   
#+ rollingAC-fixedTime
sos_ac1_fT_12 <- sos_12[,j={roll_ts(y=y, x=x, FUN=ac1, width=28*2*12)}]
sos_ac1_fT_144 <- sos_144[,j={roll_ts(y=y, x=x, FUN=ac1, width=28*2)}]
sos_ac1_fT_288 <- sos_288[,j={roll_ts(y=y, x=x, FUN=ac1, width=28)}]

#+ rollingAC-fixedTime-fig, fig.width=3.5, fig.height=6.5, fig.cap="**Figure.** Rolling window AC(1). Panels differ in the window size of the rolling average (applied before calcluating AC statistic). AC statistic is calculated over the same duration of time in each panel.", results='hide'
par(mfrow=c(3,1), mar=c(2.5, 2.0, 1, 0.25), mgp=c(1, 0.25, 0), tcl=-0.15, ps=8, cex=1)
sos_ac1_fT_12[,plot(x, y, xlab="DoY", ylab="AC(1) for Peter Chla (1 hr avg)", type='l')]
sos_ac1_fT_144[,plot(x, y, xlab="DoY", ylab="AC(1) Peter Chla (12 hr avg)", type='l')]
sos_ac1_fT_288[,plot(x, y, xlab="DoY", ylab="AC(1) Peter Chla (24 hr avg)", type='l')]
#'   
#'   
#' In this case, point plotted in all three panels are from an AR(1) statistic calculated from the same period of time. However, the first 2 panels have more data points per time period than the bottom panel, so more observations go into each AR(1) point for these panels.  
#'   
#' Unlike the "fixed-size" example previosly-shown, in this "fixed-time" approach, all 3 panels show increases in AC leading up to day ~195. One could even argue that the clearest increase is in the top panel, which had no discernable signal at all when using the previous method!  
#'   
#' Also worth noting is the difference in scaling in the y-axis. Between days 160 and 195, the top panel increases from ~0.8 to ~1.0, middle panel from 0.4 to 1.0, and the bottom panel from ~0.4 to ~1.4. (Note that usually AC is a statistic constrained to be between -1 and 1, but if there is a trend in the data, this can inflate the statistic b/c the time series is not 'stationary'; for now, I'm ignoring this.) Even though the top panel's AC has less of a range, the trend is at least as clear as for the panels showing a y-axis range that is 3-4x's bigger!  
#'   
#' I need to think more about these results, because this second outcome was not quite expected. Ideally, all 3 panels would cover the same period of time **and** have the same number of data points. A better approach to try next would be to use down-sampling (instead of averaging, just sample every Nth data point) to calculate the AC statistic. This approach could be done in a bootstrapping type approach, where every Nth sample is used, and the time series resampled N times (but not randomly, so that all data are used). E.g., sample every N points starting with the Mth observation, then every N points starting with M+1, and so on.  
#'   
#' Thus far, these results are very interesting. **There are two important messages so far.** 1) if you can collect high-frequency (which is becoming more and more common due to advances in sensor technologies) do it, because you can just increase the window size when calculating your AC EWS statistic and get a have more power for that calculation (more data points); 2) You don't have to worry too much about the effects of data aggregation when calculating your statistic --- you can use the points separately or in aggregate, and you apparently get very similar results. This is good news if you are worried about computation time, or simply don't want to worry about having to fine-tune another user-defined parameter when computing EWS. This result earns EWS a point for user-friendliness!  
#'   

#' ###Fixed Window Time & Size (subsampling)
#'   
#'   
#' When sampling a variable over a particular window of time (say, August 2017), varying the sampling resolution changes affects both time scale and sample size. It affects time scale because adjacent measurements are closer together in high resolution data. If affects sample size because there will be more samples in high resolution data.  
#'   
#' If my goal is to understand how time scale influences autocorrelation (the correlation between adjacent samples, aka AR(1) aka AC), I need to control for the effect of sample size.  
#'   
#' Aggregating data to a coarser time scale (e.g., taking a daily average) still makes use of the larger sample size by incorporating the additional samples into the aggregated statistic. The solution here is to subsample the high-frequency time series instead of aggregating it. By subsampling, we can use an original high-frequency time series to create additional time series of varying resolutions w/o incorporating variable amounts of information into each element of the coarser series.  
#'   
#' When we go to calculate a rollowing window statistic, we immediate run into an issue where it is not immediately obvious whether the "window size" should be interpreted as maining a constant number of observations (e.g., 24 data points) across time series of different resolutions, or if it should be interpreted to reflect a constant period of elapsed time (e.g., 24 hours). If we set the window to reflect a constant period of time (24 hr), then the high-resolution series will have more data point per window, inflating the sample size of the statistic applied to that window. If we set the window to reflect a constant number of observations (e.g., 24 points), then the high resolution series will not cover as much time, and therefore summary statistics will not be exposed to the same overarching dynamics when they are computed.  
#'   
#' To resolve these issues of window size across time series of varying resolution, I used a fixed-time-elapsed window, and then I reduced (subsampled) the number of data points used in the calculation of the windowed summary statistic. For example, for a time series with twice the resolution of a reference series, I would halve the number of data points used in calculating the summary statistic in each window. This process becomes extra tricky for autocorrelation --- in the previous example, one cannot simply use every other observation, because this would be akin to calculating lag-2 autocorrelation, not calculating lag-1 autocorrelation for half as many points. In order to calculate lag-1 autocorrelation I put all pairs of temporally adjacent observations in a 2-column matrix (one column for each observation forming the 'pair'), and subsampled to every other (or every nth) row. I think calculated the (auto)correlation between those two columns. I now have a lag-1 AC statistic that was claculated from half as many "pairs" of observations. Note that in this particular example all points in the subsampled time series would contribute to the AC statistic, but if 2-column matrix wasn't again subsampled then each point would contribute twice. E.g., if I had a series 1,2,3,4,5,6,7,8,9,10,11,12, I would subsample it to 1,3,5,7,9,11. Then, instead of calculating the correlation between the columns of $\begin{matrix} 1 & 3 \\ 3 & 5 \\ 5 & 7 \\ 7 & 9\\ \end{matrix}$, I'd calculate it between the columns of  $\begin{matrix} 1 & 3 \\ 5 & 7 \\ \end{matrix}$ (of course, there'd hopefully be more than 2 rows in the resultant matrix, this is just an example). Further note that in the general case of subsample to every nth point, not all points are typically included. E.g., if we were subsampling to every 3rd point b/c the high-resolution time series was sampled 3 x's as often, the subsampled matrix of the previous example would be $\begin{matrix} 1 & 3 \\ 7 & 9\\ \end{matrix}$ (again note that the number of rows depends on the length of time series).  
#' 
#' In addition controling for sample size in the summary statistic, I also needed to control for the number of windows. Even if a window covers a fixed period of time, and if the statistic applied to that window uses a fixed number of observations, a high-resolution time series will still contain more potential windows. As a result, there will be a larger number of AC summary statistics in the resultant AC time series. If trend statistics like a linear regression of Kendall's tau are applied to this high-res AC series, then it'll have more power (more AC estimates) than an AC series based on a lower-res time series. This solution is easy: when applying the rolling window, increment the window by the same interval used in the autocorrelation subsampling (normally windows are incremented by 1 time step). I.e., if the series has twice the resolution, increment the window by 2.   
#'   
#' In summary, I took the following 3 steps to control for sample size when investigating the influence of time scales on time series of different resolutions:  
#'     1. I created time series of different resolutions by using subsampling, not aggregating.  
#'     2. When calculating rolling window statistics, I subsampled again within each window, being sure to preserve the adjacency structure of the subsampled time series.  
#'     3. Lastly, I instead of incrementing the rolling window by 1, I incremented by the same interval used in subsampling.  
#'   


#+ subSamp-chla-3
sos_ss_12 <- sos[lake=="Peter", list(x=sub_samp(x=doy, n=12), y=sub_samp(x=chla, n=12))]
sos_ss_144 <- sos[lake=="Peter", list(x=sub_samp(x=doy, n=144), y=sub_samp(x=chla, n=144))]
sos_ss_288 <- sos[lake=="Peter", list(x=sub_samp(x=doy, n=288), y=sub_samp(x=chla, n=288))]

sos_ss_ac1_12 <- list()
sos_ss_ac1_144 <- list()
sos_ss_ac1_288 <- list()

#+ rollingAC-subData-fullStat-fullWindow
sos_ss_ac1_12[['fSfW']] <- sos_ss_12[,j={roll_ts(y=y, x=x, FUN=ac1, width=28*12*2)}]
sos_ss_ac1_144[['fSfW']] <- sos_ss_144[,j={roll_ts(y=y, x=x, FUN=ac1, width=28*2)}]
sos_ss_ac1_288[['fSfW']] <- sos_ss_288[,j={roll_ts(y=y, x=x, FUN=ac1, width=28)}]

#+ rollingAC-subData-subStat-fullWindow
sos_ss_ac1_12[['sSfW']] <- sos_ss_12[,j={roll_ts(y=y, x=x, FUN=ac_sub, width=28*12*2, by=1, n=24, phase=1)}]
sos_ss_ac1_144[['sSfW']] <- sos_ss_144[,j={roll_ts(y=y, x=x, FUN=ac_sub, width=28*2, by=1, n=2, phase=1)}]
sos_ss_ac1_288[['sSfW']] <- sos_ss_288[,j={roll_ts(y=y, x=x, FUN=ac1, width=28)}]

#+ rollingAC-subData-fullStat-subWindow
sos_ss_ac1_12[['fSsW']] <- sos_ss_12[,j={roll_ts(y=y, x=x, FUN=ac1, width=28*12*2, by=24)}]
sos_ss_ac1_144[['fSsW']] <- sos_ss_144[,j={roll_ts(y=y, x=x, FUN=ac1, width=28*2, by=2)}]
sos_ss_ac1_288[['fSsW']] <- sos_ss_288[,j={roll_ts(y=y, x=x, FUN=ac1, width=28)}]

#+ rollingAC-subSample-subStat-subWindow
sos_ss_ac1_12[['sSsW']] <- sos_ss_12[,j={roll_ts(y=y, x=x, FUN=ac_sub, width=28*12*2, by=24, n=24, phase=1)}]
sos_ss_ac1_144[['sSsW']] <- sos_ss_144[,j={roll_ts(y=y, x=x, FUN=ac_sub, width=28*2, by=2, n=2, phase=1)}]
sos_ss_ac1_288[['sSsW']] <- sos_ss_288[,j={roll_ts(y=y, x=x, FUN=ac1, width=28)}]

#+ rollingAC-subSample-fig, fig.width=6.5, fig.height=6, fig.cap="**Figure.** Rolling window AC(1). Panels differ in the amount of subsampling applied before calculating statistic. AC statistic is calculated over the same duration of time in each panel.", results='hide'
par(mfcol=c(3,4), mar=c(2.5, 1.75, 1, 0.25), mgp=c(1, 0.25, 0), tcl=-0.15, ps=8, cex=1)
sos_ss_ac1_12[['fSfW']][,plot(x, y, xlab="", ylab="AC(1) for Peter Chla (1-hr res)", type='l')]
mtext("Full Stat, Full Win", side=3, line=0, font=2)
sos_ss_ac1_144[['fSfW']][,plot(x, y, xlab="", ylab="AC(1) Peter Chla (12-hr res)", type='l')]
sos_ss_ac1_288[['fSfW']][,plot(x, y, xlab="DoY", ylab="AC(1) Peter Chla (24-hr res)", type='l')]

sos_ss_ac1_12[['sSfW']][,plot(x, y, xlab="", ylab="", type='l')]
mtext("Sub Stat, Full Win", side=3, line=0, font=2)
sos_ss_ac1_144[['sSfW']][,plot(x, y, xlab="", ylab="", type='l')]
sos_ss_ac1_288[['sSfW']][,plot(x, y, xlab="DoY", ylab="", type='l')]

sos_ss_ac1_12[['fSsW']][,plot(x, y, xlab="", ylab="", type='l')]
mtext("Full Stat, Sub Win", side=3, line=0, font=2)
sos_ss_ac1_144[['fSsW']][,plot(x, y, xlab="", ylab="", type='l')]
sos_ss_ac1_288[['fSsW']][,plot(x, y, xlab="DoY", ylab="", type='l')]

sos_ss_ac1_12[['sSsW']][,plot(x, y, xlab="", ylab="", type='l')]
mtext("Sub Stat, Sub Win", side=3, line=0, font=2)
sos_ss_ac1_144[['sSsW']][,plot(x, y, xlab="", ylab="", type='l')]
sos_ss_ac1_288[['sSsW']][,plot(x, y, xlab="DoY", ylab="", type='l')]

#' Of all the changes I made in this section relative to previous sections, I think the alteration with the biggest impact was not aggregating the original time series (e.g., daily average), and instead subsampling the original time series (e.g., one sample per day). In this round of results, I am most surprised that, when calculating the windowed AC statistic, using the full windowed series (much larger sample size for high-res time series) was not very different from using the subsampled series. Previously, my best guess as to why the high-res time series performed so well -- arguably *better* than the daily average -- in the fixed-time-elapsed + hourly-average scenario was because each calculation of AC had so many more observations available. This does not appear to be the case any longer. In the "Full Stat" columns of the above figure, the top row (hourly) has 24-times as many data points per AC calculation relative to the bottom row (daily). But apparently sample size supplied to the AC statistic is **not** the what was giving the high-frequency data such good performance in the early analyses.  
#'   
#' For next steps I need to 1) check my code; 2) try higher- and lower- resolutions than what I'm using here in order to "break" the AC indicator (get it to not show a warning); 3) If I can't "break" the indicator, I should apply & read about DFA (I'm thinking that if I can't break AC, then it is scale-free, though I'm pretty sure that DFA as a warning just means that the scale-free property changes [decays?] as the tipping point is approach, and in that case, DFA wouldn't tell me what I want to know); 4) after 2&3, plot a heat map of the ACF (`acf()`), because that will explicitly show me how AC is changing over time and across time scales. My guess is that this is essentially going to be identical to a heat map based on `spec.ar()`.  


#' \FloatBarrier  
#'   
#' ***  
#'   


#' #A General Method for Determining Time Scale
#' For this I think I need to check out `fractal::timeLag`! (https://www.rdocumentation.org/packages/fractal/versions/2.0-1/topics/timeLag)

#'   
#' \FloatBarrier  
#'   
#' ***  
#'   


#' #General AC as EWS

#'   
#' \FloatBarrier  
#'   
#' ***  
#'   





#' #Info
#+ systemSettings, results='markup'
Sys.time()
sessionInfo()
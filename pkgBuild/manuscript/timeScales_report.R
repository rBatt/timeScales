#' ---
#' title: "The influence of sampling frequency on statistical indicators of declining resilience"
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


#' #Functions in Developement
#' ##Plotting
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
	# l2 <- embed(x, 2)
	# ac <- cor(l2[,1], l2[,2])
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
	l2 <- stats::embed(x, 2)
	ac <- tryCatch(stats::lm(I(l2[,1]) ~ I(l2[,2]))$coef[2], error=function(cond)NA)
	
	# # option 4
	# l2 <- embed(x, 2)
	# out <- RollingWindow::RollingCorr(l2[,2], l2[,1], window=win)
	
	return(ac)
}


#' ##Rolling Windows and Aggregation
roundGrid <- function(x, frac=1){
	# if frac is 1, then place in a 1º grid
	# if frac is 0.5, then place in the 0.5º grid
	floor(x/frac)*frac+frac/2
}

roll_ts <- function(y, width=288, na.rm=TRUE, FUN=mean, x){
	buff <- rep(NA, width-1)
	mat <- embed(y, width)
	agg <- c(buff, apply(mat, 1, FUN, na.rm=na.rm))
	if(!missing(x)){
		mat2 <- embed(x, width)
		agg2 <- c(buff, apply(mat2, 1, max, na.rm=na.rm))
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
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   


#' #Data Prep
sos <- sos_data[Lake!="Tuesday"] # drop tuesday lake
setnames(sos, 
	old=c("Year", "Lake", "DoY", "DateTime", "Temp_HYLB", "Chla_Conc_HYLB", "BGA_Conc_HYLB"),
	new=c("year","lake","doy","datetime","wtr","chla","bga")
) # shorter names

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
sos_12 <- sos[lake=="Peter", agg_ts(y=chla, x=doy, width=12)]
sos_144 <- sos[lake=="Peter", agg_ts(y=chla, x=doy, width=144)]
sos_288 <- sos[lake=="Peter", agg_ts(y=chla, x=doy, width=288)]

#+ rollingAvg-chla-3-fig, fig.width=3.5, fig.height=6.5, fig.cap="**Figure.** Rolling averages of chlorophyll in Peter Lake. The window size for the rolling average varies among panels.", results='hide'
par(mfrow=c(3,1), mar=c(2.5, 2.0, 1, 0.25), mgp=c(1, 0.25, 0), tcl=-0.15, ps=8, cex=1)
sos_12[,plot(x, y, xlab="DoY", ylab="Peter Chla (1 hr avg)", type='l')]
sos_144[,plot(x, y, xlab="DoY", ylab="Peter Chla (12 hr avg)", type='l')]
sos_288[,plot(x, y, xlab="DoY", ylab="Peter Chla (24 hr avg)", type='l')]
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
sos_ac1_fS_12 <- sos_12[,j={roll_ts(y=y, x=x, FUN=ac1, width=28)}]
sos_ac1_fS_144 <- sos_144[,j={roll_ts(y=y, x=x, FUN=ac1, width=28)}]
sos_ac1_fS_288 <- sos_288[,j={roll_ts(y=y, x=x, FUN=ac1, width=28)}]

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
#' Thus far, these results are very interesting. There are two important messages so far. 1) if you can collect high-frequency (which is becoming more and more common due to advances in sensor technologies) do it, because you can just increase the window size when calculating your AC EWS statistic and get a have more power for that calculation (more data points); 2) You don't have to worry too much about the effects of data aggregation when calculating your statistic --- you can use the points separately or in aggregate, and you apparently get very similar results. This is good news if you are worried about computation time, or simply don't want to worry about having to fine-tune another user-defined parameter when computing EWS. This result earns EWS a point for user-friendliness!  
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   


#' #A General Method for Determining Time Scale

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
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
library("timeScales")

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

#+ rollingAvg-chla-3-fig, fig.width=3.5, fig.height=6.5, fig.cap="**Figure.** Rolling averages of chlorophyll in Peter Lake. The window size for the rolling average varies among panels.", results='hide'
plot_agg_sub <- function(X, steps, ln="Peter", vn="chla", stat_name="", agg_tag="avg", type=c("agg", "samp")){
	# key use of this function is that it works for multiple time scales (amounts of aggregation or subsetting)
	# will not be useful for hanging many combinations of lakes or variables
	# however, takes arguments for 
	type <- match.arg(type)
	ns <- length(steps)
	if(type=="agg"){
		par(mfrow=c(ns,1), mar=c(1.25, 1.75, 0.25, 0.25), oma=c(1,0.1,1.25,0.1), mgp=c(1, 0.15, 0), tcl=-0.15, ps=8, cex=1)
	}else{
		# par(mfcol=c(3,4), mar=c(2.5, 1.75, 1, 0.25), mgp=c(1, 0.25, 0), tcl=-0.15, ps=8, cex=1)
	}
	for(s in 1:ns){
		iname <- interval_name(steps[s])
		tdat <- X[[paste0(type,steps[s])]][lake==ln & variable==vn]
		ylab <- paste0(stat_name, ln, " ", vn, " (", iname, " ", agg_tag, ")")
		plot(tdat[,x], tdat[,y], xlab="", ylab=ylab, type='l')
	}
	mtext("Day of Year", side=1, outer=TRUE, line=0.1)
}
plot_agg_sub(X=sos_agg, steps=agg_steps)

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
#' ##Fixed Window Size (number observations)
#+ rollingAC-fixedSize
sos_agg_ac_fS <- roll_ac.sos(X=sos_agg, window_elapsed=win_days, vars=vars, lakes=lakes)


#+ rollingAC-fixedSize-figure, fig.width=3.5, fig.height=6.5, fig.cap="**Figure.** Rolling window AC(1). Panels differ in the window size of the rolling average (applied before calcluating AC statistic). AC statistic is calculated from the same number of points in each panel.", results='hide'
plot_agg_sub(X=sos_agg_ac_fS, steps=agg_steps, stat_name="AC of ")
mtext(paste0("Windows have Fixed # of Steps (", win_days, " Obs.)\nbut Vary in Time Period Covered"), side=3, line=-0.25, outer=TRUE)
#'   
#'   
#' There is minimal signal in AR(1) for the hourly data, and a very clear increase in AR(1) for the daily data; 12-hr data is somewhere in between. Each point in all three panels is autocorrelation calculated from the same number of data points, but due to differences in sampling frequency, the amount of time covered by those data points differs.  
#'   
#' Using this approach, there is a strong influence of sampling frequency on not only the AR(1) value (note y-axis scaling), but also the change in the statisic over time. Namely, the bottom panel shows a clear increase leading up to day ~195, whereas the top panel shows no trend at all.  
#'   

#' ##Fixed Window Time (amount of time elapsed)
#'   
#'   
#+ rollingAC-fixedTime
sos_agg_ac_fT <- roll_ac.sos(X=sos_agg, window_elapsed=win_days*24*60/5/agg_steps, vars=vars, lakes=lakes)

#+ rollingAC-fixedTime-fig, fig.width=3.5, fig.height=6.5, fig.cap="**Figure.** Rolling window AC(1). Panels differ in the window size of the rolling average (applied before calcluating AC statistic). AC statistic is calculated over the same duration of time in each panel.", results='hide'
plot_agg_sub(X=sos_agg_ac_fT, steps=agg_steps, stat_name="AC of ")
mtext(paste0("Windows Cover Fixed Time Period (", win_days, " days)\nbut Vary in # Observations"), side=3, line=-0.25, outer=TRUE)
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

#' ##Fixed Window Time & Size (subsampling)
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
samp_sos <- function(aggsteps){
	ss <- function(x){sub_samp(x=x, n=aggsteps, phase=1)}
	out <- sosm[,j={list(x=ss(x=doy), y=ss(x=value))},by=c("lake","variable")]
	out
}
sos_samp <- lapply(agg_steps, samp_sos)
names(sos_samp) <- paste0("samp", agg_steps)

#+ rollingAC-subData-fullStat-fullWindow
sos_ac_fSfW <- roll_ac.sos(X=sos_samp, window_elapsed=win_days*24*60/5/agg_steps, vars=vars, lakes=lakes)

#+ rollingAC-subData-subStat-fullWindow
sos_ac_sSfW <- roll_ac.sos(X=sos_samp, window_elapsed=win_days*24*60/5/agg_steps, vars=vars, lakes=lakes, subStat=TRUE)

#+ rollingAC-subData-fullStat-subWindow
sos_ac_fSsW <- roll_ac.sos(X=sos_samp, window_elapsed=win_days*24*60/5/agg_steps, vars=vars, lakes=lakes, subWindow=TRUE)

#+ rollingAC-subSample-subStat-subWindow
sos_ac_sSsW <- roll_ac.sos(X=sos_samp, window_elapsed=win_days*24*60/5/agg_steps, vars=vars, lakes=lakes, subStat=TRUE, subWindow=TRUE)

#+ rollingAC-subSample-fig, fig.width=6.5, fig.height=6, fig.cap="**Figure.** Rolling window AC(1). Panels differ in the amount of subsampling applied before calculating statistic. AC statistic is calculated over the same duration of time in each panel.", results='hide'
par(mfcol=c(4,4), mar=c(1.5, 1.75, 1, 0.25), mgp=c(1, 0.25, 0), tcl=-0.15, ps=8, cex=1)
plot_agg_sub(sos_ac_fSfW, agg_steps, agg_tag="sample", type="samp")
mtext("Full Stat, Full Window", outer=TRUE, line=-1, side=3, adj=0.065)
plot_agg_sub(sos_ac_sSfW, agg_steps, agg_tag="sample", type="samp")
mtext("Sub-Stat, Full Window", outer=TRUE, line=-1, side=3, adj=0.375)
plot_agg_sub(sos_ac_fSsW, agg_steps, agg_tag="sample", type="samp")
mtext("Full Stat, Sub Window", outer=TRUE, line=-1, side=3, adj=0.675)
plot_agg_sub(sos_ac_sSsW, agg_steps, agg_tag="sample", type="samp")
mtext("Sub-Stat, Sub Window", outer=TRUE, line=-1, side=3, adj=0.985)

#' Of all the changes I made in this section relative to previous sections, I think the alteration with the biggest impact was not aggregating the original time series (e.g., daily average), and instead subsampling the original time series (e.g., one sample per day). In this round of results, I am most surprised that, when calculating the windowed AC statistic, using the full windowed series (much larger sample size for high-res time series) was not very different from using the subsampled series. Previously, my best guess as to why the high-res time series performed so well -- arguably *better* than the daily average -- in the fixed-time-elapsed + hourly-average scenario was because each calculation of AC had so many more observations available. This does not appear to be the case any longer. In the "Full Stat" columns of the above figure, the top row (hourly) has 24-times as many data points per AC calculation relative to the bottom row (daily). But apparently sample size supplied to the AC statistic is **not** the what was giving the high-frequency data such good performance in the early analyses.  
#'   
#' For next steps I need to 1) check my code; 2) try higher- and lower- resolutions than what I'm using here in order to "break" the AC indicator (get it to not show a warning); 3) If I can't "break" the indicator, I should apply & read about DFA (I'm thinking that if I can't break AC, then it is scale-free, though I'm pretty sure that DFA as a warning just means that the scale-free property changes [decays?] as the tipping point is approach, and in that case, DFA wouldn't tell me what I want to know); 4) after 2&3, plot a heat map of the ACF (`acf()`), because that will explicitly show me how AC is changing over time and across time scales. My guess is that this is essentially going to be identical to a heat map based on `spec.ar()`.  
#'   
#' **Note:** (25-Sept-2017) I've checked my code an made some changes. I found an error in my code to subset the statistic (simple didn't do the subset, oops). Once I fixed this error, weird looking oscillations appeared in some of the output (see above figures, and below for explanation); this caught me off guard and I spent a while understanding it (and check my code many more times). Furthermore, I realized that the presence of a trend is probably having a massive effect on autocorrelation (see further below for explanation). Of the above 'next steps', I think the most relevant is the suggestion to try `spec.ar()`. I've adapted code to play around with time series resolutions easily, but before I try that, I need to remove trends from the time series (that'll probably 'break' the high-frequency statistic).  
#'   


#' ###Explanation for why sub-stat sub-window oscillates
#+ exp-oscillations-largeSim, results="hide"
y <- arima.sim(model=list(ar=0.5), n=2*28*5)
testX <- list(samp144=data.table(lake="Peter",variable="chla",x=roundGrid(seq(0.25,length(y)/2, by=0.5), 0.5), y=y))
testOut <- roll_ac.sos(X=testX, window_elapsed=win_days*24*60/5/agg_steps[3], vars=vars, lakes=lakes, subStat=TRUE)
testOut[[1]][,plot(x,y, type='l')]
#' Reproduces oscillations in simple simulated data  
#'   
#' Next use even smaller simulation, and manually code the calculations for greater transparency
#+ exp-oscillations-smallSim-Manual, results="markup"
# manually reproduce procedure for calculating sub-stat rolling window
y <- arima.sim(model=list(ar=0.5), n=80) # simulate smaller time series for easier visual inspection
y_embed_window <- embed(y, 20) # each row is a different 'window', each column different time step in that window
y_embed_window_small <- y_embed_window[1:5,] # subset to the first 5 windows for even smaller example
y_embed_window_small_embed <- lapply(unlist(apply(y_embed_window_small, 1, list), rec=F), embed, 2)
t_ac1 <- function(x){cor(x, use="na.or.complete")[2,1]}
sapply(y_embed_window_small_embed, t_ac1) # don't subset the statistic; doesn't oscillate
sub_fun <- function(x){x[c(TRUE,FALSE),]} # subset to every-other row
sub_ac1 <- function(x){t_ac1(sub_fun(x))} # calculate autocorrelation for a matrix subset to every other row
plot(sapply(y_embed_window_small_embed, sub_ac1), type='l') # statistic oscillates

print(y) # follow a value from this time series throughout the next windowed groups
print(y_embed_window_small_embed[[1]][c(TRUE,FALSE),]) # first window, column 1 is time t, column 2 is time t-1
print(y_embed_window_small_embed[[2]][c(TRUE,FALSE),]) # second window
print(y_embed_window_small_embed[[3]][c(TRUE,FALSE),]) # third window
print(y_embed_window_small_embed[[4]][c(TRUE,FALSE),]) # fourth window
print(y_embed_window_small_embed[[5]][c(TRUE,FALSE),]) # fifth window

#' Examining the time series and then the double-embeded matrices (each embedded matrix above is a 'window' where column 1 is time t and column 2 is time t-1), it is clear that the high and low values of autocorrelation appear with alternative sets of adjacent observations. If you have pairs of observations $AB, BC, CD, DE, EF, FG, GH, HI$, one window (of size 3) might have ${AB, CD, EF}$. Then the next would have ${BC, DE, FG}$. Then you'd have ${CD, EF, GH}$ for the third window, and ${DE, FG, HI}$ for the fourth. Note that the 1st and 3rd windows share ${CD, EF}$, and the 2nd and 4th share ${DE, FG}$. By contrast, adjacent windows don't have any observation pairs in common. In these examples, every-other window shares 2/3 of the previous window's pairs. If the window sizes were larger, the commonality would be even higher, yet the adjacent windows would still have 0 pairs in common. Just by chance, calculating the AC statistic for two subsets of a time series will produce 2 non-identical values, one of which must be larger than the other. For large-sized windows (subsets), removing one observation pair and adding in a new one creates a relatively small difference in the information content of the new window. Thus, the initial arbitrary large-small difference between adjacent windows is maintained as the windows are incremented, because the inremented windows have a lot of the same information content (this is why a rolling window tends to produce a smoothe time series of the applied statistic). But because the statistic is also applying an every-other subset of its own, adjacent windows actually have very different information content whereas every-other window has very similar content. As a result, you get repeated up-down cycle. Later, when I increment the window by a larger interval (instead of only sliding it by 1, slide it by the same subsetting used in the AC calculation), then the AC's subset interval and the window's subset interval are synchronized, resulting in an output that is effectively either the "peaks" or "troughs" of the jagged time series produced when the AC is subset but the window is incremented by just 1. A time series of all "peaks" or all "troughs" looks as smooth as a time series produced by a typical rolling window statistic.  



#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Achieving Stationarity
#' Here I'm going to try several approaches to making a sonde chlorophyll time series stationary while calculating AC. Initially I'll be testing these approaches on the full time series. However, what's truly relevant is that the time series is stationary within a window being analyzed. I'm using the full time series first because it's easier to just visualize one set of residuals and fit 1 model. I figure that if I can get the full tiem series to become stationary with a reasonably simple approach, then that same approach should be pretty straightforward to adapt to a 28-day window (full time series is 119 days).  
#'   
#+ test-make-stationary-data, results="hide", fig.width=4, fig.height=6
dat <- sos_samp$samp12[lake==lakes[1] & variable==vars[1]] # take the hourly data for starters
d <- dat[,y]
dts <- ts(d, deltat=1/24)

par(mfrow=c(2,1)) # prepare to visualize
dat[,plot(x, y, type='l')] # full time series
dat[1:(24*5),plot(x, y, type='l')] # 5 days worth of data
#' This time series is clearly non-stationary. Throughout the series there are 24-hour cycles. Over the full time series there isn't a single strong trend. However, if the series were to be split into rolling windows, several of those windows would show a strong trend, or perhaps a increase then decrease or vice versa.

#+ decompose-stl-show-shortTS, fig.cap="**Figure.** 5 days of the decomposed hourly time series. Blue lines are at 6am (near sunrise), red lines are at noon, at orange lines are at 7pm (near sunset)"
plot(d[1:(24*5)], type='o')
abline(v=24*(0:5)+6, col='blue') # 6am
abline(v=24*(0:5)+12, col='red') # noon
abline(v=24*(0:5)+19, col='orange') # 9pm
#' Just as a reference, it's useful to note that the fluorescence is highest at night, lowest durign the day, and that there's more volatility in the time series at night when fluorescence is high.  
#'   


#' ##Achieving Stationarity: Time Series Decomposition
#+ decompose-stl, fig.cap="**Figure.** Top panel shows the full time series. The second panel shows the contribution of seasonal component (diel cycle; a regularly repeating pattern). The third panel is the trend component, which includes local slopes as well as 'cycles' (oscillations w/o a fixed period). The fourth panel is the residuals.", fig.width=4, fig.height=6
dts_interp <- ts(approx(x=seq_along(d)[!is.na(d)], y=d[!is.na(d)], xout=seq_along(d))$y, freq=24) # need to interpolate for stl
stl_fit <- stl(dts_interp, s.window="periodic", t.window=24*3)
plot(stl_fit)

#' Plot shows that diel cycle (second panel) is relatively small in magnitude relative to local trends (third panel). In this analysis, I set the trend window to 3 days; I also tried a version when I took the log, but it didn't affect things too much. It looks like the residuals (fourth panel) might still have a bit of periodicity. Oddly the, the seasonal panel shows that the seasonality might be a bit more complex than a simple wave of period = 1 day --- i.e., every day there seem to be some small fluctuations at a higher frequency than the diel cycle. At the peaks of the time series there appears to be some fluctuation that isn't present in the troughs, and I think that this is what the decomposition is picking up.  
#'   
#' Ecological note: these peaks occur at night (fluoresence appears to be highest in dark, probably as a result of recovery from photodmg and reduction in NPQ, b/c ineffeciency should be lower at night [and thus both fluoresence and NPQ]).  
#'   
#' Back to the data analysis side of thigns, it's also worth lookingat the trend panel to understand what order of a Fourier series we might want to use if `fourier` is to be used to produce a set of covariates that could be used in a timer series model that would help us make this time series stationary. Every regular wave probably requires a sine-cosine pair. So if there's an important bump in the time series that isn't repeated, it requires an extra pair. Looking at this time series, I see about 3 big important bumps that aren't part of any regular series: there's a peak around 15, one from 25-40, and another big one around 85. This is super approximate (I don't have to know the timing, I'm just trying to come up with an informed guess as to what order of Fourier series I should try using first). Furthermore, there's a lot of little increases and decreases throughout --- I'm not sure if these are regularly spaced, but it's plausible that they could be. So maybe those could be estimated by a single sine-cosine pair. So there's the daily cycle (seasonal panel), the 3 big bumps in the trend panel, and a possible repeating set of bumps in the trend panel. So maybe ~5 sine-cosine pairs would be a reasonable place to start if trying to generate a Fourier series. If I try this and it look like I'm missing some important local trends that would invalidate stationarity for a particular window, I could increase the order of the Fourier series.  
#'   
#' Furthermore, I could try fitting the fourier series separately for each window; in that case, I would want at least 1 pair for the diel cycle (for time series with a sub-daily observation frequency), and then maybe 1 more for some sort of trend or hump in the data (shrug).  
#'   


#' ##Achieving Stationarity: Seasonal Dummy & `auto.arima()`
#+ test-make-stationary-seasonalDummy, fig.width=4, fig.height=6
aa_basline <- auto.arima(dts, max.q=0, max.p=1, max.d=0, seasonal=FALSE)

seas_reg <- forecast::seasonaldummy(dts, h=length(dts))
aa_seasDummy <- auto.arima(dts, max.q=0, max.p=1, max.d=0, seasonal=FALSE, xreg=seas_reg)

dts_1window <- ts(dts[1:(24*28)], freq=24)
seas_reg_1window <- forecast::seasonaldummy(dts_1window, h=length(dts_1window))
aa_seasDummy_1window <- auto.arima(dts_1window, max.q=0, max.p=1, max.d=0, seasonal=FALSE, xreg=seas_reg_1window)


#+ test-make-stationary-seasonalDummy-acfFig, fig.width=4, fig.height=7
par(mfrow=c(2,1), mar=c(2,2,1,0.5), ps=8, cex=1, tcl=-0.15, mgp=c(1,0.15,0))
Acf(residuals(aa_basline), na.action=na.pass, lag.max=1000)
mtext("baseline ar1 model", side=3, line=0)
Acf(residuals(aa_seasDummy), na.action=na.pass, lag.max=1000)
mtext("ar1 model with seasonal dummy in xreg", side=3, line=0)
#' This plot shows that the baseline ar1 model clearly still has a lot of autocorrelation in the residuals, and still has a daily cycle. The model with the dummy xreg also has autocorrelation in the residuals, but the diel cycle is much less apparent. This doesn't even get into the issue of local trends (or overall trend), but these results alone indicate that the introducing a seasonal xreg to an ar1 model doesn't satisfy model assumptions.

#+ test-make-stationary-seasonalDummy-residualsFig, fig.width=4, fig.height=7
par(mfrow=c(4,1), mar=c(2,2,1,0.5), ps=8, cex=1, tcl=-0.15, mgp=c(1,0.15,0))
plot(residuals(aa_basline), main="residuals from baseline ar1 model")
plot(residuals(aa_seasDummy), main="residuals full model")
plot(ts(residuals(aa_seasDummy)[1:length(dts_1window)], freq=24), main="residuals full model, only first 28 days")
plot(residuals(aa_seasDummy_1window), main="residuals from model fit to first 28 days")
#' This plot is of residuals, and helps me understand how well I did at handling those local trends ... basically I wanted to eyeball whether the residuals looked like noise, or if there was a lot of structure leftover. Unfortuneately, I think there's still a ton of structure. Probably because there's a lot more going on in here than simply the seasonality, and that's all that the xreg really addressed. Looking at the time series decomposition above suggested that this would be the case, because the diel cycle was a small portion of the decomposition relative to the local trend (or residuals).  
#'   


#' ##Achieving Stationarity: Fourier Series & `auto.arima()`
#' First, a handy function to fill-in NA's  
#'   
#+ fill-na-function
fill_na <- function(x){
	xvec <- seq_along(x)
	navec <- is.na(x)
	filled_x <- approx(xvec[!navec], y=x[!navec], xout=xvec)$y
	if(is.ts(x)){
		filled_x <- ts(filled_x, freq=frequency(x))
	}
	return(filled_x)
}
#' ###Closer Inspection, & Approaches That don't Work
#'   
#'   
#' First let's look at the spectrum for the time series, and see how much variability we can get rid of by first removing linear trend and seasonality  
#+ spec-tslmResiduals-plot, fig.height=6, fig.width=4, fig.cap="**Figure.** Top panel is spectrum. Bottom panel is residuals after removing linear trend and seasonality using tslm(), which is a really handy function."
dts_interp <- fill_na(dts) # need to interpolate for stl
d_spec <- spec.ar(dts_interp, plot=FALSE) # was messing with this to find frequency, led me to discover tslm()
par(mfrow=c(2,1))
plot(d_spec)
plot(residuals(tslm(dts ~ trend + season)))
#' After removing trend a seasonality, there is still a lot of variability. The time series still looks pretty similar.  
#'   
#'   

#' I then started playing with high-order polynomials. The above time series is obviously much longer and more complex than what we'll ever have to deal with in a single 28-day window. However, it's clear that removing a linear trend wouldn't handle the complex patterns within a 28-day window, either (namely, a quadratic, or higher polynomial ... a window with multiple "humps"). So now I'm going to play around with adding polynomial predictors to `auto.arima`.  
#+ auto.arima-polyPredict-test1, results="markup"
seas_reg <- forecast::seasonaldummy(dts, h=length(dts))
trend_var <- seq_along(dts)/(1*10^(nchar(length(dts))-1))
trend_var2 <- trend_var^2
trend_var3 <- trend_var^3
mm <- cbind(seas_reg, trend_var, trend_var2, trend_var3)
# (aa_test <- auto.arima(dts_interp, xreg=mm)) # didn't select an AR term, slected differencing
(a_test <- arima(dts_interp, xreg=mm, order=c(1, 0, 0))) # forcing model to select at AR term, and no differencing
#'   
#'   

#' It's tricky to balance the ARIMA components with constant predictor variables. In particular, the differencing. I learned a thing or two about how differencing interacts with constants in an ARIMA model to effectively remove a trend: https://www.otexts.org/fpp2/arima-modelling-in-r.html#eq:mu I knew some of this, but I hadn't realized how the constant and differencing degree interacted. There's the problem that `Arima` won't allow you to slect d > 2 and still have a constant --- therefore, you can't capture a quadratic or higher pattern with this approach. Furthermore, I'm not sure how this works with the quadratic terms having different coefficients. But I decided to play around anyway.
#+ auto.arima-differencing, fig.height=6, fig.width=4, fig.cap="**Figure** A model that does not allow differencing but includes a linear and quadratic temporal predictor. This plot is the time series decomposition of that model's residuals. "
plot(stl(residuals(auto.arima(dts_interp, max.q=0, max.p=1, max.d=0, seasonal=FALSE, allowdrift=F, xreg=seas_reg)), s.window='periodic'))
#+ auto.arima-differencing2, fig.height=6, fig.width=4, fig.cap="**Figure** Trying to use differencing and constants to capture polynomial patterns. The figure is of the decomposition of the residuals of a model with ARIMA(1,1,0) errors and seasonal dummy predictors. The trend component still shows lots of pattern", results="markup"
# https://www.otexts.org/fpp2/arima-modelling-in-r.html#eq:mu
(aa_test2 <- auto.arima(dts_interp, max.q=0, max.p=1, max.d=7, seasonal=FALSE, allowdrift=TRUE, xreg=seas_reg))
plot(stl(residuals(aa_test2), s.window='periodic'))

# # try similar to aa_test2, but try forcing a quad trend?
# # auto.arima excludes this possibility b/c of interest in forecasting, but not as dangerous when not forecasting
# aa_test3 <- auto.arima(dts_interp, max.q=0, max.p=1, max.d=7, seasonal=FALSE, allowdrift=TRUE, xreg=cbind(seas_reg, trend_var, trend_var2))
# plot(stl(residuals(aa_test3), s.window='periodic'))
#' These plots are time series decompositions of residuals of two models: the first doesn't have differencing/ drift, the second does. Allowing for differencing does a pretty good job of absorbing a lot of the low-frequency pattern in the "trend" component of the decomposition. I also tried adding linear and quadratic temporal trends to the `xreg`, but this didn't change the decomposition output much.  
#'   
#' The problem with this approach becomes apparent when you look at the AR(1) coefficient of the models. For example, `aa_test2` has an AR1 coefficient of **`r unname(round(aa_test2$coef['ar1'], 3))`**. The coefficient is negative. This seems suspicious to me. *A priori*, I know that a first order polynomial isn't going to capture much of the local nonstationarity (patterns that would appear as nonstationary in a windowed time series), and my hunch is that the AR1 term is taking on a weird value to account for this limited description of low-frequency pattern.  
#'   
#' To really demonstrate why this matters, I need to examine a window-sized (28 days) subset of the time series, and analyze that. This will provide a much clearer and transferable/ to-the-point test of how to handle nonstationarity when calculating AC as an early warning statistic.  
#'   


#' ###Test the Fourier Approach
#' Before I get into using a subset of the time series, let's test the Fourier series as a way of removing the seasonality and other wobbles. The main challenge with this approach arises when you don't know the frequency of the oscillations *a priori*, or if there aren't oscillations so much as local trends. Ultimately, I decided against this approach, but I'll make note of my experience with it.  
#'   
#' First thing to note is that I have to define a `msts` --- multi-seasonality time series. To use the Fourier function from the `forecast` package, I need to define all possible frequencies ahead of time. Unfortunately, `forecast::findfrequency` only tries to estimate a dominant frequency; moreover, the source code for `forecast::findfrequency` reveals that this approach to finding the frequency depends on first making the time series stationary (the function uses `forecast::tslm`) ... the very problem we're trying to address in the first place. You can begin to see the problem I'm facing here. Anyway, I took a whack at it using daily, weekly, fortnightly, and monthly frequencies. I have no idea how meaningful these are, but ti was a starting place.
#+ auto.arima-fourier, fig.width=4, fig.height=6, fig.cap="**Figure.** Decomposition of an auto.arima model where some of the xreg are from a Fourier series."
dmsts <- msts(fill_na(d), seasonal.periods=c(24, 24*7, 24*14, 24*28), ts.freq=24) # define multiple periods
fourier_reg <- forecast::fourier(dmsts, K=c(6, 4, 4, 4)) # create Fourier series, defining the # of sine-cosine pairs per frequency
# aa_four <- auto.arima(dmsts, seasonal=FALSE, xreg=fourier_reg)
aa_four <- auto.arima(dmsts, seasonal=FALSE, xreg=fourier_reg, max.q=0, max.p=1, max.d=0)
aa_four_resid <- residuals(aa_four)
stl_aa_fit <- stl(aa_four_resid, s.window="periodic", t.window=24*3)
plot(stl_aa_fit)
#' The result of this model fit, as you can see, is really disappointing: the residuals have tons of structure and remove far less of the seasonal component than using the seasonal dummy variable. To be fair, there were 23 columns to the matrix of seasonal dummy variables, whereas I only had 12 columns dedicated to sine or cosine waves with a daily period. So this is with ~half the parameters. But still.  
#'   
#' Furthermore, the trend component of the decomposition reveals a lot of pattern. This is not much better than some of the less successful approaches that I tried previously. The patterns are similar, but just of smaller overall magntidue. So a bit better, but not really doing a great job. If the Fourier series can be made to work, it would likely require paying a lot of careful attention to the frequencies and the number of Fourier pairs. Fitting that, or trying all combinations, would be intensive.  
#'   
#'   


#' ##Achieving Stationarity: 1 Month of Data fit with Polynomials + Arima()
#' It'll be helpful to fit the models using 1 month of data. Based on previous experiments, the seasonal dummy term and a polynomial might do the trick. Fitting a polynomial to the full time series is not a very transferrable idea --- it would have to be a pretty high order polynomial, much higher than I would want for a 28-day window (most likely). Nonetheless, I know ahead of time that I wouldn't want the same order polynomial for all windows; so in this sense, the approach I take here needs to be both somewhat conservative (fitting high-order polynomials is sketchy) and generic (don't know what order is appropriate ahead of time).  
#'   
#+ stationarity-subsetTimeSeries, fig.width=4, fig.height=6, fig.cap="**Figure.** Decomposition of a month of data that has a hump in the middle of it."
# example of a subset of d that would be best-fit by a quadratic ... what happens?
dts_quad <- fill_na(ts(d[540:1211], freq=24)) # 28 days with a peak in the middle
# dts_quad
# plot(dts_quad)
plot(stl(dts_quad, s.window="per"))
#' This chunk of data is a good example of a potential challenge. It has a big hump in the middle. So I'd suspect a 2nd order polynomial. Of course, I also need to take care of the daily cycle.  
#'   

#+ trend-xreg-function
# This function creates dummy variables for polynomials ... to be used in regression
trend_xreg <- function(exp.order, y){
	x <- seq_along(y)
	hnames <- paste0("trend",1:exp.order)
	
	raise_order <- function(o){
		matrix(x^o, ncol=1)
	}
	omat <- sapply(1:exp.order, raise_order)
	dimnames(omat) <- list(NULL, hnames)
	return(scale(omat))
}


#+ select-sar-poly-function
# AR Model with X-reg (e.g., Seasonal Dummy) and Selected Degree of Polynomial Trend
select_ar_xreg_poly <- function(series, poly_max, xreg){
	
	if(missing(xreg)){
		xreg <- 1
	}else{
		xreg <- cbind(1, xreg)
	}
	
	fit_each <- function(poly_order){
		if(poly_order > 0){
			forecast::Arima(series, order=c(1,0,0), xreg=cbind(xreg, trend_xreg(poly_order, series)), include.mean=FALSE)
		}else{
			forecast::Arima(series, order=c(1,0,0), xreg=cbind(xreg), include.mean=FALSE)
		}
	}
	
	mod_list <- lapply(1:poly_max, fit_each)
	aiccs <- sapply(mod_list, function(x)x$aicc)
	
	# get_ar <- function(x){
# 		round(unname(x$coef['ar1']), 3)
# 	}
# 	ars <- sapply(mod_list, get_ar)
	
	aiccs_adjusted <- aiccs + seq_along(aiccs)*10 # to make sure a simpler model is chosen unless a more complicated one is more than 10 aicc points better (per extra parameter)
	mod_out <- mod_list[[which.min(aiccs_adjusted)]]
	return(mod_out)
}

#+ getARCoef-function
get_ar <- function(x){
	round(unname(x$coef['ar1']), 3)
}
get_ar_name <- function(x){
	paste("AR =", get_ar(x))
}

#+ ar-seasonDummy-lineTrend, fig.width=4, fig.height=4, fig.cap="**Figure** 28-day time series of Peter chlorophyll (black line) and model fitted values (red). Model is an AR(1) with fixed regressors of a seasonal dummy and a linear trend."
seas_reg_quad <- forecast::seasonaldummy(dts_quad, h=length(dts_quad))
xreg_linear <- cbind(1, seas_reg_quad, trend_xreg(1, dts_quad))
a_test_linear <- forecast::Arima(dts_quad, order=c(1,0,0), xreg=xreg_linear, include.mean=FALSE)
plot(dts_quad)
lines(fitted(a_test_linear), col='red')
mtext(text=get_ar_name(a_test_linear), side=3, line=-1, adj=0.95)
# plot(forecast(a_test_linear, xreg=cbind(1, seas_reg_quad, trend1=(1:length(dts_quad)))))
#' This figure shows the fitted values and AC coefficient when a linear trend is estimated alongside the AR(1) and a seasonal dummy. It's important to recognize that these are 1-step-ahead predictions (for the ar1 component). Import to note is how the predicted values are overestimates at the beginning and the end of the time series. Furthermore, note that the AR1 coefficient is `r get_ar(a_test_linear)`. These features contrast between this model witha  linear trend and an upcoming model with an AICc-selected quadratic model.  
#'   
#'   
#+ ar-seasonDummy-lineTrend-residualSTL, fig.width=4, fig.height=6, fig.cap="**Figure.** Decomposition of residuals of a AR1 model with seasonal dummies and linear trend as predictors. Point is that the linear trend is insufficient."
plot(stl(residuals(a_test_linear), s.window="per"))
#' This figure furthers the point that the linear trend is insufficient as a xreg in the ARIMA model --- a decomposition of the residuals of this model still shows quadratic nonstationarity in the "trend" panel.  
#'   


#+ ar-seasnDummy-polyTrend-select, fig.width=4, fig.height=4, fig.cap="**Figure.**  28-day time series of Peter chlorophyll (black line) and model fitted values (red). Model is an AR(1) with fixed regressors of a seasonal dummy and an AICc-selected polynomial trend (quadratic selected)."
seas_reg_quad <- forecast::seasonaldummy(dts_quad, h=length(dts_quad))
polyAR_fit <- select_ar_xreg_poly(dts_quad, poly_max=6, xreg=seas_reg_quad)

plot(dts_quad)
lines(fitted(polyAR_fit), col='red')
mtext(text=get_ar_name(polyAR_fit), side=3, line=-1, adj=0.95)
#' This figure shows that the quadratic model does a better job of estimating chlorophyll at the start and finish of the time series (the model with a linear trend overestimated chlorophyll at these points). Furthermore, choosing a linear trend vs a quadratic trend has a noticable impact on the AR estimate: the linear trend model estimated the AC as `r get_ar(a_test_linear)`, whereas the quadratic trend model is estimating the AC as `r get_ar(polyAR_fit)`.  
#'   
#'   

#+ ar-seasonDummy-polyTrend-resids, fig.width=4, fig.height=6, fig.cap="**Figure.**  Decomposition of residuals of a model that has AR(1) with fixed regressors of a seasonal dummy and an AICc-selected polynomial trend (quadratic selected)."
plot(stl(residuals(polyAR_fit), s.window="per"))
#' Here we see a decomposition of the residuals of an AR model with a term for a quadratic trend. The "trend" portion of the decomposition shows that, while there is a substantial amount of low-frequency of variability, the pattern appears to be stationary (note: I don't mean this literally. The residuals still contain pattern --- ar, ma, etc. I would estimate a full-blown ARIMA model for this, but I don't want to complicate the interpretation of the lag-1 autocorrelation). The "seasonal" panel shows some pattern, but the magnitude is tiny, indicating that most of the strictly periodic pattern has been accounted for.  
#'   
#' One issue I may face with a higher-frequency time series is that the 24-hour seasonality will be computationally intensive to capture with dummy variables (~288 parameters, vs currently I only need ~24). In those cases, I might have to adapt the Fourier approach to capture the seasonality.  
#'   

#' ##Achieving Stationarity: Remove 24-hr Seasonality w/ Fourier
#' With a higher sampling frequency, removing seasonality with seasonal dummy variables becomes computationally intensive b/c each observation per day gets its own parameter. So 24 observations per day (hourly) is managable, but 288 observations per day (5-minutes) is not. To solve this problem, I'll invest more in figuring out how to use the Fourier series to deal with daily seasonality.  
#'   
#' Previously, my primary disappointment in the Fourier was my hopes that it could remove lower-frequency variation. This problem has been primarilty solved through the windowing of the time series to 28 days, and the addition of polynomial detrending (with the degree of the polynomial selected by AICc, which protects against overfitting, the primary problem of polynomials).  
#'   
#+ ar-fourier-poly, results="markup", fig.width=4, fig.height=6, fig.cap="**Figure** Using Fourier to account for seasonality instead of seasonal dummy variables. Decomposition of residuals of model with AR term, quadratic trend, and Fourier series for seasonal regressors."
select_ar_fourier_poly <- function(series, poly_max, K_max, freq=NULL){
	if(!is.ts(series)){
		stopifnot(!is.null(freq))
		series <- ts(series, freq=freq)
	}
	# stopifnot(is.ts(series))
	
	# Determine the best polynomial order
	# **Just for the purpose of fitting the best fourier order**
	# Will re-fit ideal poly order along with AR term once fourier order determined
	poly_list <- rep(NA, poly_max)
	for(p in 1:poly_max){
		treg <- trend_xreg(p, series)
		t_mod <- lm(series ~ treg)
		poly_list[p] <- AIC(t_mod, k=20) # introduce a very strong penalty for additional parameters
	}
	temp_poly <- which.min(poly_list)
	
	# Functions for Fourier
	getF <- function(tf){forecast::fourier(series, tf)}
	fit_each_noPoly <- function(fourier_order){
		forecast::Arima(series, order=c(1,0,0), xreg=cbind(1, getF(fourier_order)), include.mean=FALSE)
	}
	fit_each_poly <- function(fourier_order){
		forecast::Arima(series, order=c(1,0,0), xreg=cbind(1, getF(fourier_order), trend_xreg(temp_poly, series)), include.mean=FALSE)
	}
	
	# Get AIC for fourier
	fourier_orders <- 1:K_max
	if(temp_poly>0){
		mod_list <- lapply(fourier_orders, fit_each_noPoly)
	}else{
		mod_list <- lapply(fourier_orders, fit_each_poly)
	}
	aiccs <- sapply(mod_list, function(x)x$aicc)
	# aiccs_adjusted <- aiccs + seq_along(aiccs) #*10
	best_fo <- fourier_orders[which.min(aiccs)]
	t_reg <- getF(best_fo)
	
	# Reselect poly order now that we have fourier selected
	mod_out <- select_ar_xreg_poly(series, poly_max, xreg=t_reg)
	return(mod_out)
}

(o <- select_ar_fourier_poly(dts_quad, 4, 6)) # a test run
plot(stl(residuals(o), s.window='per'))


#' #Stationary Rolling AC
#+ full-stationary-ac-fourier-poly, fig.width=4, fig.height=4, fig.cap="**Figure.** Time seies of rolling autocorrelation. No subsetting of the statistic or the window, but the time series is hourly. The autocorrelation was calculated using `Arima()`, with `xreg` set to a Fourier series and a polynomial trend. The size the Fourier series and order of the polynomial was selected by AIC for each window. Maximum size and order were each set to 6. This took a pretty long time to run.", results='markup'
Sys.time()
ac1.stationary <- function(series, poly_max=6, K_max=6, freq=NULL){
	get_ar(select_ar_fourier_poly(series, poly_max=poly_max, K_max=K_max, freq=freq))
}
ac1.stationary.model <- function(series, poly_max=6, K_max=6, freq=NULL){
	select_ar_fourier_poly(series, poly_max=poly_max, K_max=K_max, freq=freq)
}

y <- ts(sos_samp$samp12[lake==lakes[1] & variable==vars[1], y], freq=24*60/5/agg_steps[2])
x <- sos_samp$samp12[lake==lakes[1] & variable==vars[1], x]

blah <- roll_ts(y=y, width=win_days*24*60/5/agg_steps[2], FUN=ac1.stationary, x=x, freq=24*60/5/agg_steps[2]) # agg_steps[2] corresponds to samp12
plot(blah, type='l')
Sys.time()

#' This look a lot like the normal calculation of autocorrelation, except that 1) the AC doesn't increase as early in the time series (here the increase doesn't start until around 190, previously increase started around 165); 2) after about day 210 the AC declines down to original levels, whereas previously the decline to original levels didn't occur until after day 240.  
#'   
#'   

#' ##Stationary Rolling AC: Inspect Windows
#+ rolling-window-inspect-for-stationarity, results='markup'
Sys.time()
y <- ts(sos_samp$samp12[lake==lakes[1] & variable==vars[1], y], freq=24*60/5/agg_steps[2])
x <- sos_samp$samp12[lake==lakes[1] & variable==vars[1], x]
rollAC_12_64 <- roll_ts(y=y, width=win_days*24*60/5/agg_steps[2], by=64, FUN=ac1.stationary.model, x=x, freq=24*60/5/agg_steps[2])
rollAC_12_64[!is.na(x),plot(x, sapply(y, get_ar), type='l')]

plot_roll_pieces <- function(X, plot_interactive=FALSE){
	X2 <- X[!is.na(x)]
	
	pfun <- function(x, day_plot){
		if(missing(day_plot)){day_plot <- x[1, x]}
		x2 <- x[x==day_plot]
		mod <- x2[,y][[1]]
		day <- x2[,x]
		
		xreg_names <- colnames(mod$xreg)
		ts_freq <- frequency(mod$x)
		xreg_pred <- ts(mod$xreg %*% mod$coef[names(mod$coef)%in%xreg_names], freq=ts_freq, start=day)
		
		trend_label <- paste(xreg_names[grepl("trend[0-9]{1}", xreg_names)], collapse=" + ")
		ar_label <- get_ar_name(mod)
		
		par(mfrow=c(3,1), mar=c(2.5,2.5,1,0.5), ps=8, cex=1, mgp=c(1.25,0.25,0), tcl=-0.25)
		plot(ts(mod$x, freq=ts_freq, start=day), lwd=1) # observed time series
		lines(ts(fitted(mod), freq=ts_freq, start=day), col='red') # fitted time series
		lines(xreg_pred, col='blue') # pattern from fourier + polynomial trend
		lines(ts(mod$x, freq=ts_freq, start=day), lwd=4, col=adjustcolor('black', 0.2))
		mtext(ar_label, side=3, line=0.05, adj=0.95, font=2)
		mtext(trend_label, side=3, line=0.05, adj=0.05)
		
		plot(ts(residuals(mod), freq=ts_freq, start=day), type='l')
		
		Acf(residuals(mod), na.action=na.pass, ylab="ACF of Residuals")
	}
	
	nmax <- nrow(X2)
	if(interactive() & plot_interactive){
		day_ind <- 1
		input <- -1
		while(!input%in%c("q","Q")){
			pfun(X2, X2[day_ind, x])
			input <- readline("Enter day number, hit enter, or type 'q' or 'Q'\n\n")
			if(input=="" & day_ind+1>nmax){
				input <- 'q'
			}else{
				if(input==""){
					day_ind <- day_ind + 1L
				}
			}
		
		}
	}else{
		day_inds <- seq(1, nmax, length.out=6)
		holder <- lapply(day_inds, function(ind)pfun(X2, X2[ind, x]))
	}

}
Sys.time()

plot_roll_pieces(rollAC_12_64)

#' #Detrender
#+ seriousDetrending
detrendR <- function(x, max_poly=6, max_fourier=6, max_interaction=3){
	
	full_poly <- trend_xreg(max_poly, x)
	full_fourier <- forecast::fourier(x, max_fourier)
	
	get_interaction <- function(p, f){
		int_list <- structure(vector("list", ncol(p)), .Names=colnames(p))
		for(pp in 1:ncol(p)){
			PP <- p[,pp] #pick this trend
			pfmat <- matrix(NA, ncol=ncol(f), nrow=nrow(f))
			for(ff in 1:ncol(f)){	
				pfmat[,ff] <- f[,ff]*PP
			}
			colnames(pfmat) <- colnames(f)
			int_list[[pp]] <- pfmat
		}
		return(int_list)
	}
	full_interaction_list <- get_interaction(p=full_poly, f=full_fourier)
	
	get_mm <- function(p, f, i){
		pmm <- full_poly[,1:p, drop=FALSE]
		get_fmm <- function(x, f){
			x[,1:min(2*f, ncol(x))]
		}
		fmm <- get_fmm(full_fourier, f)
		
		if(i == 0){
			mm <- cbind(intercept=1, pmm, fmm)
		}else{
			imm_list <- lapply(full_interaction_list[1:i], get_fmm, f=f)
			for(l in 1:length(imm_list)){
				dimnames(imm_list[[l]])[[2]] <- paste(dimnames(imm_list[[l]])[[2]], names(imm_list)[l], sep='-')
			}
			imm <- do.call(cbind, imm_list)
		
			mm <- cbind(intercept=1, pmm, fmm, imm)
		}
		return(mm)
	}
	
	pfi_combos0 <- expand.grid(p=1:max_poly, f=1:max_fourier, i=0:max_interaction)
	limit_interactionOrder <- pfi_combos0[,3] <= pfi_combos0[,1] #& pfi_combos0[,3] <= pfi_combos0[,2]
	pfi_combos <- pfi_combos0[limit_interactionOrder,]
	
	
	mod_aicc <- function(m){
		X <- do.call(get_mm, as.list(pfi_combos[m,]))
		tX <- t(X)
		K <- ncol(X)
		bhat <- solve(tX%*%X)%*%tX%*%Y
		Yhat <- X%*%bhat
		resid <- Y - Yhat
		s2 <- sd(resid)
		NLL <- -sum(dnorm(resid, mean=0, sd=s2, log=TRUE))
		correction <- (2*K*(K+1))/(nrow(X) - K - 1)
		AICc <- 2*K + 2*NLL + correction
		return(AICc)
	}
	Y <- fill_na(as.numeric(x)) # use notation that's easier to recognize in the following
	AICcs <- sapply(1:nrow(pfi_combos), mod_aicc)
	# for(m in 1:nrow(pfi_combos)){
# 		mod_aicc(m)
# 	}
	X <- do.call(get_mm, as.list(pfi_combos[which.min(AICcs),]))
	resid <- c(Y - X%*%(solve(t(X)%*%X)%*%t(X)%*%Y))
	resid <- ts(resid, freq=frequency(x))
	return(resid)
}
detrendR(dts, 10, 6, 4)

# Rprof(tmp <- tempfile(), memory.profiling=TRUE)
# detrendR(dts_quad)
# Rprof()
# summaryRprof(tmp, memory="both")


#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Effect of Trend on AC over Different Time Horizons
#' *Update on Previous Thinking:* I think the high-frequency time series show increases in autocorrelation for a long time horizon but not a short time horizon because of a trend in the data, and that the ac1 estimator doesn't detrend.
#+ trend-impact-AC, results="markup"
set.seed(42)
x <- 1:100
y0 <- as.numeric(arima.sim(model=list("ar"=0.25), n=length(x)))
y <- y0 + 0.075*x

# correlation of full time series w/o trend, as baseline
# should be ~0.25
emat0 <- embed(y0, 2)
cor(emat0[,1], emat0[,2]) # ~0.18 in one run, 0.29 in another; a bit off, but close

# correlation of full time series w/ trend
emat <- embed(y, 2)
cor(emat[,1], emat[,2]) # is ~0.76 in one run, 0.86 in another; much higher than true 0.25

# correlation over full series, but with fewer pairs of observations
randInd <- sort(sample(x-1, 28))
cor(emat[randInd,1], emat[randInd,2]) # is ~0.90, almost same as correlation over full series

# correlation over a shorter period of time
cor(emat[1:28,1], emat[1:28,2]) # is ~0.37, much lower than over full series

# detrended time series, then calculate correlations again
detrend <- function(x, time=1:length(x)){
	stats::residuals(stats::lm(x~time))
}
emat2 <- embed(detrend(y),2)
cor(emat2[,1], emat2[,2]) # is ~0.29, nearly identical to the correlation of y0, which did not have trend
cor(emat2[randInd,1], emat2[randInd,2]) # this coefficient seems to vary a lot among random seeds :p
cor(emat2[1:28,1], emat2[1:28,2]) # ~0.33

# another approach to calculating detrended ar1:
arima(y, order=c(1,0,0), xreg=x) # estimates ar1 as ~0.29, and trend as ~0.08 --- pretty close on both
#' These results indicate that a linear trend can inflate the autocorrelation estimate, and that this inflation is greater when a longer portion of the time series (with the trend) is analyzed. As a result, estimates of first order autocorrelation will appear to be different over short time horizons (over which the influence of the trend is smaller, and the AC less inflated) and long time horizons (over which the trend is much more apparent, and the AC more inflated).  
#'   
#' 



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
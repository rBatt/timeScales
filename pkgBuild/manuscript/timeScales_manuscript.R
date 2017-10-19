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
# 	"~/Documents/School&Work/epaPost/timeScales/pkgBuild/manuscript/timeScales_manuscript.R",
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


#' #Options
#+ options
win_days <- 28 # window size in days covered
agg_steps <- c(1, 12, 288, 288*2) # step sizes for aggregation
lakes <- c("Peter","Paul") # can be vector; lakes to analyze (Paul, Peter)
vars <- "chla" # can be vector; variables to analyze (wtr, bga, chla)


#'   
#' \FloatBarrier  
#'   
#' ***  
#'   


#' #Data Prep
#' ##Subset, Restructure, Define as ts Object
#+ data-prep-basic
# ---- drop tuesday lake ----
sos <- sos_data[Lake!="Tuesday"]

# ---- shorter names ----
setnames(sos, 
	old=c("Year", "Lake", "DoY", "DateTime", "Temp_HYLB", "Chla_Conc_HYLB", "BGA_Conc_HYLB"),
	new=c("year","lake","doy","datetime","wtr","chla","bga")
)

# ---- ensure numeric, re-structure data set ----
sos[,bga:=as.numeric(bga)]
sosm <- melt(sos, id.vars=c("year","lake","doy","datetime"))[variable%in%vars & lake%in%lakes]

# ---- make measured values of class "ts" with frequency = 288 samples per day ----
set_ts <- function(y, x, freq=288){
	ts(y, freq=288, start=x)
}
sosm[, value:=set_ts(y=value, x=doy[1]), by=c("year","lake","variable")]

# ---- grab range limits (primarily for plotting) ----
doy_range <- sos[,range(doy, na.rm=TRUE)]
chla_range <- sos[,range(chla, na.rm=TRUE)]
# bga_range <- sos[,range(bga, na.rm=TRUE)]
# wtr_range <- sos[,range(wtr, na.rm=TRUE)]

#' ##Aggregate Data for Typical Rolling Window Calculation
#+ data-prep-aggregation
agg_sos <- function(aggsteps){
	out <- sosm[,j={agg_ts(y=value, x=doy, width=aggsteps)},by=c("lake","variable")]
	out
}
sos_agg <- lapply(agg_steps, agg_sos)
names(sos_agg) <- paste0("agg", agg_steps)
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   


#'  #Chlorophyll Time Series for each Lake
#+ chla-timeSeries-figure, fig.width=3.5, fig.height=5, fig.cap="**Figure.** High frequency chlorophyll (chla, micrograms per liter) time series in Peter (red) and Paul (blue) lakes in 2015.", results='hide'
par(mfcol=c(2,1), mar=c(2, 2.0, 1, 0.25), mgp=c(1, 0.25, 0), tcl=-0.15, ps=8, cex=1)
sos[lake=="Paul", plot(doy, chla, xlim=doy_range, col="black", type='l', xlab="", ylab="Chlorophyll")]
mtext("Paul Lake (Reference)", side=3, line=-0.1, adj=0.05, font=2)
sos[lake=="Peter", plot(doy, chla, xlim=doy_range, col="black", type='l', xlab="Day of year", ylab="Chlorophyll")]
mtext("Peter Lake (Manipulated)", side=3, line=-0.1, adj=0.05, font=2)
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


#' #ACF Plots for Chla in each Lake
#+ functions-plotting
plot_acf <- function(ln=c("Paul","Peter"), v=c("chla", "bga"), na.action=na.exclude, lag.max=288*12, ...){
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
#+ chlorophyll-acf-figure, fig.width=3.5, fig.height=5, fig.cap="**Figure.** Autocorrelation function (ACF) of chlorophyll a (indicator of algal biomass) from Peter Lake (manipulated) and Paul Lake (reference).", results='hide'
par(mfrow=c(2,1), mar=c(2, 2.0, 0.25, 0.25), mgp=c(1, 0.25, 0), tcl=-0.15, ps=8, cex=1)
plot_acf(ylab="Paul Lake Chlorophyll ACF", main="")
plot_acf(ln='Peter', ylab="Peter Lake Chlorophyll ACF", main="")
#' Autocorrelation is time scale dependent in both the manipulated and the reference lake. 

#' #Changes in Autocorrelation Over Time and Across Time Scales
#' ##Rolling Window Autocorrelation for Select Time Scales
#+ rollingWindowAC-calculation
steps_per_day <- 60*24/(5 * agg_steps) # observations per day = (60 minutes / 1 hour) * (24 hours / 1 day) * (1 observation / 5*n min)
steps_per_window <- steps_per_day*win_days # steps per window = (n steps / 1 day) * (n days / 1 window)
AC_list <- roll_ac.sos(sos_agg, window_elapsed=steps_per_window, vars=vars, lakes=lakes, DETREND=TRUE, by=c(24, 2, 1, 1))

#+ rollingWindowAC-PaulPeterDifference, fig.width=6, fig.height=6, fig.cap="**Figure.** Rolling windows of first-order autocorrelation from detrended chlorophyll time series. Blue lines are from Paul Lake (reference), red lines are Peter Lake (manipulated). In the second column, the black lines represent the difference (Peter - Paul) between the lines in the first column (positive values indicate that autocorrelation was higher in Peter than in Paul). Each row of the figure has a different sampling frequency."
plotac <- function(X, ...){
	X <- copy(X)
	ylim <- X[,range(y, na.rm=TRUE)]
	
	ydiff <- X[lake=="Peter", y] - X[lake=="Paul", y]
	zdiff <- data.table(lake="zdiff", variable=X[,variable[1]], x=X[lake==lake[1], x], y=ydiff)
	X <- rbind(X, zdiff)
	
	ul <- X[,unique(lake)]
	for(l in 1:length(ul)){
		# tX <- X[lake==ul[l]]
# 		tcol <- c("Paul"="blue","Peter"="red", "zdiff"="black")[ul[l]]
# 		if(ul[l]=="Paul"){
# 			plot(tX[,x],tX[,y], type='l', col=tcol, ylim=ylim, ...)
# 		}else if(ul[l]=="Peter"){
# 			lines(tX[,x],tX[,y], col=tcol)
# 		} else if(ul[l]=="zdiff"){
# 			plot(tX[,x],tX[,y], type='l', col=tcol, ...)
# 		}
		
		dud <- X[lake==ul[l],j={
			tcol <- c("Paul"="blue","Peter"="red", "zdiff"="black")[lake[1]]
			if(lake[1]=="Paul"){
				plot(x,y, type='l', col=tcol, ylim=ylim, ...)
			}else if(lake[1]=="Peter"){
				lines(x, y, col=tcol)
			} else if(lake[1]=="zdiff"){
				p2 <- function(..., ylab, ylab2="") plot(..., ylab=ylab2)
				p2(x,y, type='l', col=tcol, ...)
				# abline(h=0, lty="dashed")
			}
		
			NULL
		}]
	}
	
	invisible()
}

ylabs <- paste0("Chl-a AR(1) (", sapply(agg_steps, interval_name), ")")

xlabs <- rep("", length(agg_steps))
xlabs[length(agg_steps)] <- "Day of Year"

par(mfrow=c(length(agg_steps),2), mar=c(2,2,0.5,0.5), cex=1, tcl=-0.15, mgp=c(1,0.2,0), ps=8)
invisible(mapply(plotac, X=AC_list, ylab=ylabs, xlab=xlabs))







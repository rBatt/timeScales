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
#'     fig_caption: true
#' geometry: margin=1.0in
#' lineno: true
#' lineSpacing: false
#' titlesec: true
#' documentclass: article
#' placeins: true
#' ---

#+ deleted-pandoc-headers, include=FALSE, echo=FALSE
# #'      pandoc_args: [
# #'      "--chapters"
# #'      ]


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
# 	"~/Documents/School&Work/pinskyPost/spatialDiversity/pkgBuild/manuscript/timeScales_report.R",
# 	output_format=o_f,
# 	output_dir='~/Documents/School&Work/pinskyPost/spatialDiversity/pkgBuild/manuscript',
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

#' #Fluoresence Time Series

#'   
#' \FloatBarrier  
#'   
#' ***  
#'   


#' #Influence of Sampling Frequency on AC

#'   
#' \FloatBarrier  
#'   
#' ***  
#'   


#' #Influence of Sampling Frequency on AC as EWS

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


#' #Functions in Developement

#'   
#' \FloatBarrier  
#'   
#' ***  
#'   


#' #Info
#+ systemSettings, results='markup'
Sys.time()
sessionInfo()
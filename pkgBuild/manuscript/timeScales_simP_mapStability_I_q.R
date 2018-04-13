#' ---
#' title: "Map the stability properties of a lake P model for different values of P input (I) and an exponent (q)"
#' author: "Ryan Batt"
#' date: "2018-04-10"
#' abstract: I have been varying the parameter I, P input rate, as a bifurcation parameter. When doing all past simulations, I set q, and exponent in the model, to a fairly large value (10). An ecologically reasonable value is between 2 and 10, so the value I chose is a bit on the high end. Larger values of q cause the system to shift more abruptly, whereas small values of q cause the state variables to trend more strongly before a somewhat smaller jump. So I chose a large value of q to make the "regime shift" aspect of the simulation visually apparent. However, there may have been further consequences for the dynamical stability of the simulated time series. The objective of this report is to precisely identify the equilibria and critical points (critical values of I) for different values of q. This analysis will help produce a function that can suggest a reasonable range of I to simulate across when using this model to test early warning statistics, and should also help me decide what value of q to use in the paper.
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
# 	"~/Documents/School&Work/epaPost/timeScales/pkgBuild/manuscript/timeScales_simP_mapStability_I_q.R",
# 	output_format=o_f,
# 	output_dir='~/Documents/School&Work/epaPost/timeScales/pkgBuild/manuscript',
# 	clean = TRUE
# )

Sys.setenv(PATH=paste0("/Library/TeX/texbin:",Sys.getenv("PATH")))
opts_chunk$set(
	fig.path = 'timeScales_simP_mapStability_I_q/', 
	cache.path='timeScales_simP_mapStability_I_q/',
	echo=TRUE, 
	include=TRUE, 
	cache=FALSE,
	autodep=TRUE,
	results='asis',
	warning=FALSE,
	fig.show="hold",
	fig.lp = if(o_f=="html_document"){"**Figure.**"}else{NULL}
)

#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Solve for dX/dt as a Function of X
#' 
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
#' \end{array}
#'   
#' It is Eq 7 above that I will be working from for the rest of this report. It gives me the derivative of X as a function of X (and parameters).  When this equation equal zero, there is a (un)stable point. When the derivative of Eq 7 is 0 AND Eq 7 itself is 0, that's a critical point.  
#'   


#' #Functions for dX/dt=f(X)
#' Note that for these functions, I will be including "commented-out" code prior to them that will serve as a template for generating the documentation for the function if I choose to include them in the timeScales R package.
#' ##Function to calculate dX/dt
#+ dXdt_function
# Put whole model in terms of X (no M)
# gives us dX/dt as a function of X when dM/dt is 0

# Change in X with respect to time as a function of X
# 
# Calculates the change in water P per unit time as a function of water P and parameters
# 
# @param X numeric, water P
# @param I numeric, P input to lake
# @param q numeric, exponent affecting sharpness of transition
# 
# @return a numeric value indicating dX/dt
# @seealso \code{\link{modelDeterministicXM}}
# @export
dX_dt_ofXIq <- function(X, I, q){
	pars <- unlist(formals(modelDeterministicXM)[c("s", "m", "r", "h", "b")])
	for(i in 1:length(pars)){assign(names(pars)[i], unname(pars)[i])}

	# h <- 0.15
	# s <- 0.7
	# m <- 2.4
	# b <- 0.001
	# r <- 0.019
	
	R <- X^q/(m^q + X^q)
	dXdt <- I - X*(s+h) + r*R*((s*X)/(b+r*R))
	return(dXdt)
}
dX_dt_ofXIq(1, 1, 10)

#' ##Function for Plotting dX/dt vs X
#+ plot_dXdt_function
# Plot dXdt vs X
# 
# Plot the change in water P vs water P
# 
# @param Xvals numeric vector of X (water P) values
# @param Ivals numeric vector of I (P input; bifurcation parameter) values
# @param q numeric scalar of exponents
# @param ... arguments to be passed to \code{\link{modelDeterministicXM}}
# 
# @return no value returned, but as a byproduct q plots, each with length(Ivals) lines, is produced. 
plot_dXdt_X <- function(Xvals, Ivals, q, ...){
	dXdt_q10 <- outer(Xvals, Ivals, FUN=dX_dt_ofXIq, q=q)
	cols <- viridis::viridis(n=ncol(dXdt_q10))
	ylim <- range(dXdt_q10)
	ylim[1] <- max(-10, min(-0.5, ylim[1]))
	ylim[2] <- min(10, max(0.5, ylim[2]))
	plot(NA, ylim=ylim, xlim=range(Xvals), xlab="X", ylab="dX/dt", ...)
	abline(h=0, lty=2, col='gray')
	
	for(j in 1:ncol(dXdt_q10)){
		if(any(abs(dXdt_q10[,j])>20)){
			bp <- which.min(diff(dXdt_q10[,j]))
			inds <- list(1:bp, (bp+1):nrow(dXdt_q10))
			for(i in 1:2){
				lines(Xvals[inds[[i]]], dXdt_q10[inds[[i]],j], col=cols[j], pch=20, cex=0.25)
			}
		}else{
			lines(Xvals, dXdt_q10[,j], col=cols[j], pch=20, cex=0.25)
		}
	}
	legend("topright", pch=20, col=cols, legend=paste("I =", round(Ivals,4)), bg='white', x.intersp=0.5, y.intersp=0.5, inset=c(-0.1, -0.1), xpd=TRUE)
}

#' #Map dX/dt vs X, varying I and q
#+ figure1-dXdt-vs-X-qI, fig.width=3.5, fig.height=7, fig.cap="**Figure 1.** Changes in water phosphorus (X) per unit time as a function of X and I (input of P to lake) and q (an exponent in the model which indicates the sharpness of the transition). Different colors of lines indicate different values of I. Different panels show different values of q. When dX/dt=0, there is a stable or unstable point; the point is stable if 0 is crossed from above (a negative slope in dX/dt), and unstable if crossed from below (a positive slope in dX/dt). If the slope (derivative) of dX/dt is 0 AND dX/dt is 0, then there is a critical point at that value of X (and parameter combination). For example, when q=8, there is a critical point near X=7 and I= ~0.75; if I is decreasing, and usntable and stable point collide, if I is increasing, a saddle-node point emerges (I think)."
Xvals <- seq(0,7, length.out=250)
Ivals <- seq(0, 3, length.out=5)
# dev.new(width=3.5, height=7)
par(mfrow=c(4, 1), mar=c(1.25,2,1,1.5), cex=1, ps=8, mgp=c(1, 0.25, 0), tcl=-0.25, oma=c(0.75, 0.1, 0.1, 0.1))
# qs <- seq(2, 10, by=2)
qs <- c(2, 5, 8, 10)
for(j in 1:length(qs)){
	plot_dXdt_X(Xvals, Ivals, q=qs[j], main=paste("q =", qs[j]))
}
mtext("X", side=1, line=1, outer=FALSE)

#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Phase Portrait
#+ phasePortrait, fig.width=6, fig.height=6, fig.cap="**Figure 2.** A phase portrait of the system for varying values of P input when q=8. The vector field (indicating the direction and speed that the system moves through phase space at that point) is represented by gray arrows. Nullclines are represented red and blue lines, indicating where dX/dt and dM/dt are equal to zero, respectively.  Trajectories starting at arbitrary initial points (open diamonds) and continuing the along the accompanying solid black line indicate how the system moves from the initial point through phase space for 50 years. Equilibria are indicated by points: solid filled circle is a stable node, an 'X' is a saddle point. An equilibrium occurs whereever the nullclines cross. The different panels correspond to different values of P loading (I). "
par(mfrow=c(2,2), mar=c(2,2,1,0.5), mgp=c(1,0.25,0), tcl=-0.15, cex=1, ps=8)
Is <- c(0.75, 1, 1.25, 1.5)
for(i in 1:length(Is)){
	timeScales::phasePortrait(pars=c(I=Is[i]), nFlow=10, addLeg=(i==2))
	mtext(paste0("I = ",Is[i]), side=3, line=0, adj=0, font=2)
}


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

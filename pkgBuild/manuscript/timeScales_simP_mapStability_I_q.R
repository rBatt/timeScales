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
#' \text{differentiating Eq. 7 with respect to } X\\[0.5em]
#' f''(X) = \frac{rsX^q(b(m^q(q+1)+X^q)+rX^q)}{(bm^q+bX^q+rX^q)^2} -s -h   & (8)\\
#' \end{array}
#'   
#' It is Eqs 7 and 8 above that I will be working from for the rest of this report. Eq 7 gives me the derivative of X as a function of X (and parameters).  When this equation equal zero, there is a equilibrium. When Eq 8 is 0 AND Eq 7 itself is 0, that's a critical point.  
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
dX_dt_ofXI <- function(X, I, pars){
	parsF <- unlist(formals(modelDeterministicXM)[c("s", "m", "r", "h", "b", "q")])
	if(missing(pars)){
		pars <- parsF
	}else{
		pars <- c(pars, parsF[!names(parsF)%in%names(pars)])
	}
	# pars <- unlist(formals(modelDeterministicXM)[c("s", "m", "r", "h", "b")])
	for(i in 1:length(pars)){assign(names(pars)[i], unname(pars)[i])}

	# h <- 0.15
	# s <- 0.7
	# m <- 2.4
	# b <- 0.001
	# r <- 0.019
	
	R <- X^q/(m^q + X^q)
	dXdt <- I - X*(s+h) + r*R*((s*X)/(b+r*R))
	# I - X*(s+h) + r*(X^q/(m^q + X^q))*((s*X)/(b+r*(X^q/(m^q + X^q))))
	return(dXdt)
}
dX_dt_ofXI(1, 1, pars=c(q=10))

#' ##Function to Calculate d2X/dt2
#+ secondDeriv, results='markup'
d2Xdt <- function(X, pars){
	parsF <- unlist(formals(modelDeterministicXM)[c("s", "m", "r", "h", "b", "q")])
	if(missing(pars)){
		pars <- parsF
	}else{
		pars <- c(pars, parsF[!names(parsF)%in%names(pars)])
	}
	with(as.list(pars), {
		# (((b*r*s*q*X^q*m^q)/(X^q+m^q)^2)+((r*s*b*X^q)/(m^q+X^q))+((r^2*s*X^(2*q))/(m^q+X^q)^2))/(b^2 + ((2*b*r*X^q)/(m^q+X^q)) + ((r*X^q)/(m^q+X^q))^2) - s - h
		(r*s*X^q*(b*(m^q*(q+1)+X^q)+r*X^q))/(b*m^q+b*X^q+r*X^q)^2 - s - h
	})
}

uniroot.all(d2Xdt, c(0,10), n=1E4) # X values when I hits a critical value (where f'(dX/dt)==0)
uniroot.all(dX_dt_ofXI, c(0,10), I=1, pars=c(q=8), n=1E4) # equilibrium val of X when I = 1
uniroot.all(dX_dt_ofXI, c(0,10), I=0.997593522, pars=c(q=8), n=1E6) # at this val of I, two equilibria are very close to colliding; collision would happen as I is increasing
uniroot.all(dX_dt_ofXI, c(0,10), I=0.519496874, pars=c(q=8), n=1E6) # different value of I, but again, two equilibria very close; collision would happen as I is decreasing



#' ##Function to Find Critical Values
#+ function-findCriticalValues, results='markup'
findCrit <- function(pars, critRange=c(0.01,10), tol=.Machine$double.eps^0.5, nGrid=1E6, xRange=c(0,100)){
	parsF <- unlist(formals(modelDeterministicXM)[c("s", "m", "r", "h", "b", "q")])
	if(missing(pars)){
		pars <- parsF
	}else{
		pars <- c(pars, parsF[!names(parsF)%in%names(pars)])
	}
	
	# xRange <- c(dX_dt_ofXI(), dX_dt_ofXI())
	# xRange <- c(0, 100)
	
	x_targets <- uniroot.all(d2Xdt, xRange, n=nGrid)
	target_dist <- function(I, target){
		with(as.list(pars),{
			# uniroot.all(dX_dt_ofXI, c(0,10), I=I, q=q, n=1E4)
			# outer(uniroot.all(dX_dt_ofXI, c(0,10), I=I, q=q, n=1E4), x_targets, FUN="-")
			# sum(outer(uniroot.all(dX_dt_ofXI, c(0,10), I=I, q=q, n=1E4), target, FUN="-")^2)
			diffs <- outer(uniroot.all(dX_dt_ofXI, xRange, I=I, pars=c(q=8), n=nGrid), target, FUN="-")
			diffs[which.min(abs(diffs))]
		})
	}
	
	
	# testI <- seq(0.1,2,by=0.01)
	# tdist <- vector('numeric', length(testI))
	# for(i in 1:length(testI)){
	# 	tI <- testI[i]
	# 	tdist[i] <- target_dist(tI, target=x_targets[1])
	# }
	# plot(testI, tdist) # I cannot use uniroot() b/c there are severe discontinuities; as uniroot.all warns, it is really bad at finding 0's that just barely touch the 0 line; I think the algorithm looks for a change in sign or something
	# uniroot(target_dist, interval=c(critRange[1],critRange[2]), target=xt, maxiter=1E4, tol=.Machine$double.eps/2)
	
	abs_target_dist <- function(I, target){abs(target_dist(I=I, target=target))}
	# optim(0.5, abs_target_dist, target=x_targets[1])
	# optimize(f=abs_target_dist, interval=critRange, target=x_targets[1], tol=.Machine$double.eps^0.5)
	
	criticalValue <- vector('numeric', length(x_targets))
	for(i in 1:length(x_targets)){
		xt <- x_targets[i]
		criticalValue[i] <- optimize(f=abs_target_dist, interval=critRange, target=x_targets[i], tol=tol)$minimum
	}
	
	return(criticalValue)
	# diff(sort(root))*0.5*c(-1,1)+sort(root) # a good range of critical values
}


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
	dXdt_q10 <- outer(Xvals, Ivals, FUN=dX_dt_ofXI, pars=c(q=q))
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

#' #Identify Critical Values and Plot f(X)=dX/dt and f'(X)
#+ findCriticalValues, result='markup'
critVals <- sort(findCrit())
critVals # these are the critical values of I --- at these values, the X equilibria emerge/ collid (which depends on whether I is increasing or decreasing)

#+ figure2-fX-fprimeX, fig.height=6, fig.width=3.5, fig.cap="**Figure 2.** Black line is the rate of change in water phosphorus (dX/dt) vs X when the rate of change in sediment phosphorus is 0. The blue line is the derivative of the black line with respect to X."
par(mfrow=c(2,1), mar=c(2,2,0.75,0.5), mgp=c(1,0.25,0), tcl=-0.25, ps=8, cex=1)
Xvals <- seq(0,8, by=0.01)
ddXdt <- d2Xdt(Xvals)
for(i in 1:length(critVals)){
	dXdt <- dX_dt_ofXI(Xvals, I=critVals[i], pars=c(q=unlist(formals(modelDeterministicXM)[c("q")])))
	ylim <- range(c(dXdt, ddXdt))
	plot(Xvals, dXdt, ylim=ylim, type='l', xlab='Water P', ylab="f(X) or f '(X)")
	lines(Xvals, ddXdt, col='blue', type='l')
	legend("topright", lty=1, col=c("black","blue"), legend=c("f(X)=dX/dt", "f '(X)"))
	abline(h=0, lty=2, col='gray')
	grid()
	mtext(paste0("I = ", round(critVals[i],6)), adj=0.1, font=2)
}



#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Phase Portrait
#+ figure3-phasePortrait, fig.width=6, fig.height=6, fig.cap="**Figure 3.** A phase portrait of the system for varying values of P input when q=8. The vector field (indicating the direction and speed that the system moves through phase space at that point) is represented by gray arrows. Nullclines are represented red and blue lines, indicating where dX/dt and dM/dt are equal to zero, respectively.  Trajectories starting at arbitrary initial points (open diamonds) and continuing the along the accompanying solid black line indicate how the system moves from the initial point through phase space for 50 years. Equilibria are indicated by points: solid filled circle is a stable node, an 'X' is a saddle point. An equilibrium occurs whereever the nullclines cross. The different panels correspond to different values of P loading (I). "
par(mfrow=c(2,2), mar=c(2,2,1,0.5), mgp=c(1,0.25,0), tcl=-0.15, cex=1, ps=8)
# Is <- c(0.75, 1, 1.25, 1.5)
Is0 <- diff(sort(critVals))*0.5*c(-1,1)+sort(critVals) # a good range of critical values
# Is <- round(c(Is0[1], critVals, Is0[2])-0.1,2)
Is <- round(sort(c(outer(critVals, c(-0.1, .1), FUN='-'))), 2)
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
Sys.time()
sessionInfo()

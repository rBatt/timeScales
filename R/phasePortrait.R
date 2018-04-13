#' Phase Portrait
#' 
#' Plot a phase portrait of a dynamical system (defaults to lake eutrophication model) that includes the vector field, nullclines, example trajectories, and the identification and classification of equilibria.
#' 
#' @param I numeric scalar for P input
#' @param pars optional named numeric vector for additional parameters to be passed to \code{func}
#' @param func, function that returns local derivatives, defaults to \code{\link{mDXM}}, which is just a formatting wrapper for \code{\link{modelDeterministicXM}}; currently only works for default due to use of \code{\link{getRoot_df}} in \code{\link{stabClass}}
#' @param x.lim,y.lim minimum and maximum values for each of the two states variables to be used when plotting and searching for equilibria
#' @param add local, if FALSE (default) makes a new plot
#' @param nFlow square root of the number of arrows to draw in the vector field
#' @param nNull numeric scalar, number of points used in drawing the nullclines
#' @param t.end number of time steps (in model time units) to use when simulating the trajectories
#' @param inits matrix of initial values to use for trajectories
#' @param addLeg if TRUE, adds a legend
#' @param legPos character position of legend; see \code{graphics::legend}
#' 
#' @details This function is largely a wrapper for \code{phaseR::flowField}, \code{phaseR::trajectory}, \code{phaseR::nullclines}, and \code{phaseR::stability}. This function does not call \code{phaseR::stability} directly, and instead uses \code{\link{stabClass}}.
#' 
#' @return a data.frame with stability information; see \code{\link{stabClass}}
#' @examples 
#' shortInits <- eval(formals(phasePortrait)$inits)[1:3,]
#' pars <- c(q=8,s=0.7,h=0.15, b=0.001)
#' (phasePortrait(
#' 	I=0.75, pars=pars,
#' 	inits=shortInits, nFlow=3, nNull=10,t.end=10
#' ))
#' @export
phasePortrait <- function(I, pars, func=mDXM, x.lim=c(0,15), y.lim=c(0,1000), add=FALSE, nFlow=20, nNull=300, t.end=50, inits=matrix(c(0.5,800, 0.1,950, 3,700, 1.5,600, 0.5,100, 2,150, 5,75, 4.5,115, 1,300, 2,400, 4,500, 10,300, 7,200), byrow=TRUE, ncol=2), addLeg=FALSE, legPos="topright"){
	requireNamespace("phaseR", quietly=TRUE)
	if(missing(I)){
		if('I' %in% names(pars)){
			I <- unname(pars["I"])
		}else{
			stop("Provide a value for I, P input")
		}
	}else{
		if(!'I' %in% names(pars)){
			pars <- c(I=I, pars)
		}
	}
	if(missing(pars)){pars <- NULL}
		
		ff <- phaseR::flowField(func, x.lim=x.lim, y.lim=y.lim, parameters=c(pars), add=add, points=nFlow, xlab="Water P (X)", ylab="Sediment P (M)")
		nulls <- phaseR::nullclines(func, x.lim=x.lim, y.lim=y.lim, parameters=c(pars), add=TRUE, points=nNull)
		traj <- phaseR::trajectory(func, y0=inits, t.end=t.end, parameter=c(pars), colour=rep("black", nrow(inits)), pch=23)
		stab <- stabClass(I=I, pars=pars)
		stabPCH <- c(
			"Saddle" = 4,
			"Indeterminate" = 11,
			"Stable node" = 19,
			"Unstable node" = 21,
			"Stable focus" = 10,
			"Unstable focus" = 8,
			"Centre" = 14
		)
		points(stab[,c("X","M")], pch=stabPCH[stab[,"classification"]])
		
		if(addLeg){
			nd <- !duplicated(stabPCH[stab[,"classification"]])
			upch <- stabPCH[stab[,"classification"]][nd]
		
			legText <- c("dX/dt=0", "dM/dt=0", "trajectory", names(upch))
			cex <- par("cex")
			text.font <- NULL
			tw <- max(abs(strwidth(legText, units = "user", cex = cex, font = text.font)))*0.9
			
			legend(legPos, inset=-c(0.0, 0.0), lty=c(1,1,1, rep(-1, length(upch))), pch=c(NA,NA,23,upch), col=c("red","blue","black",rep('black',length(upch))), legend=legText, ncol=2, bg=adjustcolor('white',0.5), xpd=TRUE, merge=FALSE, x.intersp=0.8, y.intersp=0.5, text.width=tw)
		}
		
		
		invisible(stab)
}

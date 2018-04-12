# See simLakeP.R for the functions being tested/ played with in this script
# This script was written while writing those functions, 
# in order to test them and as a way of determining what functions would be needed/ useful.


# ---- Find roots and eigenvalues for a range of loading and water P ----
initialValues <- getInit(18, Mrange=c(200, 600), Mn=6)
rootGrid_I <- getRoot_df(initialValues)
eigGrid_I <- getEigs_df(rootGrid_I)


# ---- Plot the roots, using the eigenvalues to determine stability and reflect stable-unstable in point type ----
# set point type based on stability
princEigPos <- Re(eigGrid_I[,c("X1")]) > 0
secEigPos <- Re(eigGrid_I[,c("X2")]) > 0
eitherEigPos <- princEigPos | secEigPos
eigPCH <- rep(20, nrow(eigGrid_I))
eigPCH[eitherEigPos] <- 21

# make plot
par(mar=c(2,2,0.5,2), mgp=c(1,0.25,0), tcl=-0.25, ps=8)
plot(eigGrid_I[,"I"], eigGrid_I[,"X"], col="blue", type='p', xlab="P Loading (g/m^2)", ylab="Water P (g/ m^2)", pch=eigPCH)
par(new=TRUE)
plot(eigGrid_I[,"I"], eigGrid_I[,"M"], col="red", xaxt='n', yaxt='n', type='p', xlab="", ylab="", pch=eigPCH)
axis(side=4)


# ---- show point migration towards stability ----
nTime <- 250
dT <- 1/4
stateArray <- array(dim=c(nrow(eigGrid_I), 3, nTime), dimnames=list(NULL,variable=c("I","X","M"),time=c()))
dStateArray <- stateArray
stateArray[,,1] <- data.matrix(eigGrid_I[,c("I","init.X","init.M")])
stateArray[,"I",] <- stateArray[,"I",1]
stateArray[,"M",] <- stateArray[,"M",1] #eigGrid_I[,"M"] #eigGrid_I[,"M"]/3
dStateArray[,"I",] <- stateArray[,"I",1]
for(ti in 2:nTime){
	dState_dt <- t(apply(stateArray[,,ti-1], 1, function(x)modelDeterministicXM(state=x[c("X","M")], pars=c(x["I"]))))
	dStateArray[,c("X","M"),ti-1] <- dState_dt
	stateArray[,c("X","M"),ti] <- stateArray[,c("X","M"),ti-1] + dState_dt*dT
}

plotState <- function(stateArray, t){
	plot(stateArray[,"I",t], stateArray[,"X",t], pch=20, ylim=range(stateArray[,"X",], na.rm=TRUE))
}

# nGif <- 75
# nSlow <- 20
# gifSeq <- c(1:nSlow, trunc((seq(from=nSlow+1, to=nTime, length.out=nGif-nSlow))))
# gifSeq[1] <- 1
# animation::saveGIF({
# 	for(g in 1:nGif){
# 		if(!(complete.cases(stateArray)[gifSeq[g]])){next}
# 		plotState(stateArray, t=gifSeq[g])
# 		mtext(paste("t =", gifSeq[g]), side=3, adj=0.9, line=0)
# 	}
# }, movie.name="~/Desktop/test_approachEqui.gif")
nGif <- 75
nSlow <- 20
gifSeq <- c(1:nSlow, trunc((seq(from=nSlow+1, to=nTime, length.out=nGif-nSlow))))
gifSeq[1] <- 1
base_file <- "~/Desktop/test_approachEqui/"
if(!dir.exists(base_file)){dir.create(base_file)}
gFileList <- vector("character", nGif)
for(g in 1:nGif){
	if(!(complete.cases(stateArray)[gifSeq[g]])){next}
	g_file <- paste0(base_file,formatC(g, width=nchar(nGif), flag="0"), ".png")
	gFileList[g] <- g_file
	png(file=g_file, res=150, units='in', width=5, height=5)
	plotState(stateArray, t=gifSeq[g])
	mtext(paste("t =", gifSeq[g]), side=3, adj=0.9, line=0)
	dev.off()
}
system_comm <- paste0("convert -delay 50 -loop 0 ", base_file, "*.png ", "~/Desktop/equilibria.gif")
system(system_comm)
unlink(base_file, rec=TRUE)

dev.new()
par(mfrow=c(3,3), mar=c(2,2,0.5,0.5), ps=8, cex=1, mgp=c(1, 0.25, 0), tcl=-0.2)
time_seq <- round(seq(from=1, to=nTime, length.out=9))
for(ts in 1:length(time_seq)){
	plot(stateArray[,"X",time_seq[ts]], stateArray[,"M",time_seq[ts]], xlab="X (water P)", ylab="M (mud P)")
}




uniqueI <- unique(dStateArray[,"I",1])
i1 <- uniqueI[which.min(abs(uniqueI-1.5))]
state_i1_0 <- stateArray[dStateArray[,"I",1]==i1,"X",1]
state_i1 <- state_i1_0[order(state_i1_0)]
dState_i1 <- dStateArray[dStateArray[,"I",1]==i1,"X",1][order(state_i1_0)]
dev.new()
plot(state_i1, dState_i1, type='l')
abline(h=0)
abline(v=c(0.5, 2.15), lty='dashed')



# below will eventually be stochastic model
# function(N=200, dt=0.1, F=14.6, c=0.00115, b=0.001, h=0.15, m=2.4, r=0.019, q=8, s=0.7, sigma=0.01, lambda=0.35){
#
# 	R <- function(X){(X^q/(m^q + X^q))}
#
# 	dZ <- rnorm(N)
# 	Z <- cumsum(dZ)
# 	H <- function(t){
# 		Zt <- Z[t]
# 		exp(lambda*Zt-((t*lambda^2)/2))
# 	}
#
# 	dW_dt <- rnorm(N)
#
# 	# dU_dt <- F - c*U*H
# 	# dX_dt <- c*U*H - (s+h)*X + r*M*R(X) + sigma*M*R(X)*dW_dt
# 	# dM_dt <- s*X - b*M - r*M*R(X) - sigma*M*R(X)*dW_dt
#
# 	MRX <- M*R(X)
# 	sigmaMRX <- sigma*MRX
# 	rMRX <- r*MRX
# 	dU_dt <- F - c*U*H
# 	dX_dt <- c*U*H - (s+h)*X + rMRX + sigmaMRX*dW_dt # need to subset noise to just at time t
# 	dM_dt <- s*X - b*M - rMRX - sigmaMRX*dW_dt
#
# }

# ======================
# = Try Phase Portrait =
# ======================
library(deSolve)
initialValues <- getInit(2, Mrange=c(200, 600), In=2)
mDXM_jac <- function(t, state, pars, ...){
	modelDeterministicXM(state, pars, ...)
}


# times <- seq(0, 100, by=1/24)
# state <- c(X=initialValues[1,"X"], M=initialValues[1,"M"])
# o <- vode(state, times, mDXM, parms=c(I=initialValues[1,"I"]), q=8)
# out <- as.data.frame(vode(state, times, mDXM, parms=c(I=initialValues[1,"I"]), q=8))
#
#
# phasePortrait <- function(states, pars=c(I=1), i1=5, Times, ...){
# 	stopifnot(is.matrix(states))
# 	if(missing(Times)){
# 		# Times <- c(seq(0,3,by=0.01), 3^(seq(1,12,length.out=400))) #seq(0, 20E3, length.out=1E2) #seq(0, nt, by=dt)
# 		Times <- seq(0, 100, by=0.1)
# 	}
# 	mDXM <- function(t, state, pars, ...){
# 		list(modelDeterministicXM(state, pars, ...))
# 	}
# 	out <- vector("list", nrow(states))
# 	for(i in 1:nrow(states)){
# 		st <- states[i,]
# 		out[[i]] <-  as.data.frame(ode(st, Times, mDXM, parms=pars, ...))
# 		# out[[i]] <-  as.data.frame(ode(s, Times, mDXM, parms=pars))
# 	}
#
# 	allDists <- c()
# 	for(i in 1:nrow(states)){
# 		x <- out[[i]][[2]]
# 		y <- out[[i]][[3]]
# 		distTrav <- sum(sqrt(diff(x)^2 + diff(y)^2))
# 		allDists[i] <- distTrav
# 	}
# 	nArrows <- (allDists>median(allDists))*2+1
#
# 	xlim <- range(sapply(out, function(x)range(x[[2]])))
# 	ylim <- range(sapply(out, function(x)range(x[[3]])))
# 	plot(NA, xlim=xlim, ylim=ylim, xlab=names(out[[1]])[2], ylab=names(out[[1]])[3])
#
# 	for(i in 1:nrow(states)){
# 		x <- out[[i]][[2]]
# 		y <- out[[i]][[3]]
# 		lines(x, y)
#
# 		cumDistTrav <- c(0,cumsum(sqrt(diff(x)^2 + diff(y)^2)))
# 		distTrav <- cumDistTrav[length(cumDistTrav)]
#
# 		# probs <- c(0.01, 0.5, 0.99)
# 		probs <- c(0, 0.25, 0.75)[1:nArrows[i]]
# 		qx <- cumDistTrav/distTrav
# 		inds <- sapply(probs, function(x)which.min(abs(qx-x)))
# 		for(p in 1:length(probs)){
# 			ind <- inds[p]
# 			shape::Arrows(x[ind], y[ind], x[ind+1], y[ind+1], arr.length=0.2, segment=FALSE, code=2, arr.adj=1)
# 		}
#
# 		# rx <- rev(x)
# # 		ry <- rev(y)
# 		# shape::Arrows(rx[i1+1], ry[i1+1], rx[i1], ry[i1], col='red', arr.length=0.2)
# 	}
# 	invisible(NULL)
# }
# # states <- rbind(
# # 	as.matrix(expand.grid(X=seq(0.05, 1, length.out=3), M=seq(600, 800, length.out=3)))
# # 	,as.matrix(expand.grid(X=seq(1, 5, length.out=3), M=seq(200, 400, length.out=3)))
# # )
# # states <- as.matrix(expand.grid(X=seq(0.05, 10, length.out=3), M=seq(100, 1000, length.out=5)))
# # states <- as.matrix(expand.grid(X=c(seq(0.05,1,length.out=2), seq(1.5, 10, length.out=3)), M=c(seq(100, 500, length.out=3), seq(600, 1000, length.out=3))))
#
# states <- as.matrix(expand.grid(X=c(0.05, 1, 2, seq(3, 15, length.out=3)), M=seq(100, 1000, length.out=10)))
# # phasePortrait(states, pars=c(I=0.5), q=5, i1=5, nt=2E2)
#
#
#
# qs <- c(5, 8, 10)
# Is <- c(0.5, 0.75, 1, 1.2, 1.3, 1.5)
#
# dev.new()
# par(mfrow=c(length(qs),length(Is)), mar=c(2,2,0.5,0.5), mgp=c(1,0.25,0), tcl=-0.25, cex=1, ps=8, oma=c(0.1, 0.5, 0.5, 0.1))
#
# for(j in 1:length(qs)){
# 	tq <- qs[j]
#
# 	for(i in 1:length(Is)){
# 		ti <- Is[i]
# 		phasePortrait(states, pars=c(I=ti), q=tq, i1=2, nt=2E2)
#
# 		if(j==1){mtext(paste0("I = ", ti), side=3, line=0.25, font=2, cex=1.3)}
# 		if(i==1){mtext(paste0("q = ", tq), side=2, line=1.75, font=2, cex=1.3)}
#
# 	}
# }

# getRoot_df(cbind(states,I=1.2), q=10)
# points()



nearby <- function(x, adj=c(-0.1, 0.1)){
	a <- adj # + 1
	n <- x[0,]
	na <- length(a)
	origs <- do.call("rbind", replicate(na, x, simplify=FALSE))#[order(rep(1:nrow(x), na)),]
	
	for(j in 1:ncol(x)){
		rx <- diff(range(x[,j]))
		if(!rx>0){
			a <- adj+1
			v <- pmax(0,c(outer(x[,j],a, FUN="*")))
		}else{
			v <- pmax(0,c(outer(x[,j],a*rx, FUN="+")))
		}
		
		o <- origs
		o[,j] <- v
		n <- rbind(n, o)
	}
	return(n)
}

# rs <- getRoot_df(cbind(states,I=1.2), q=10)#[,c("X","M")]
# urs <- rs[!duplicated(paste0(round(rs[,"X"],6),round(rs[,"M"],6))), c("I","X","M"), drop=FALSE]
# urs
# e0 <- getEigs_df(urs, q=10)
# eigs <- cbind(e0[,c("I","X","M")], X1_r=Re(e0[,"X1"]), X1_i=Im(e0[,"X1"]), X2_r=Re(e0[,"X2"]), X2_i=Im(e0[,"X2"]))
#
# class1 <- c()
# class2 <- c()
# for(i in 1:nrow(eigs)){
# 	r1 <- eigs[i,"X1_r"]
# 	r2 <- eigs[i, "X2_r"]
# 	i1 <- eigs[i, "X1_i"]
# 	i2 <- eigs[i, "X2_i"]
# 	if(r1<0 & r2<0){
# 		# is a node or focal point
# 		class1[i] <- "stable"
# 	}else if(sum(c(r1,r2)<0)==1){
# 		# is a saddle
# 		class1[i] <- "saddle"
# 	}else{
# 		# is unstable
# 		class1[i] <- "unstable"
# 	}
#
# 	if(i1!=0 | i2!=0){
# 		class2[i] <- "spiral"
# 		if(Conj(i1)==i2){
# 			class2[i] <- "Hopf"
# 		}
# 	}else{
# 		class2[i] <- "node"
# 	}
# }
#
# eig_bg <- c("stable"="black", "unstable"="gray", "saddle"="red")
# eig_col <- c("node"="black", "spiral"="blue", "Hopf"="green")
#
#
# # near_urs <- rbind(cbind(I=urs[,"I"], nearby(urs[,-1], adj=c(-0.15,0.15))), cbind(I=urs[,"I"], nearby(urs[,-1], adj=c(-0.05,0.05))), cbind(I=urs[,"I"], nearby(urs[,-1])))
# adj_list <- list(c(-0.05,0.05),  c(-0.01, 0.01))
# near_urs <- do.call("rbind", lapply(adj_list, nearby, x=urs[complete.cases(urs),-1, drop=FALSE]))
# dev.new()
# phasePortrait(near_urs, pars=c(I=1.2), q=10, i1=1, Times=c(seq(0,100,by=0.1), 2^seq(log(100+1)/log(2),15,length.out=1E3)))
#
# points(urs[,c("X","M")], bg=eig_bg[class1], pch=21, col=eig_col[class2])
# legend("topright", pch=21, col=c(rep(NA,3),eig_col), pt.bg=c(eig_bg, rep(NA,3)), legend=c(names(eig_bg), names(eig_col)), ncol=2)

# ==============
# = try phaseR =
# ==============
library(phaseR)
library(timeScales)

mDXM <- function(t, y, parameters, ...){
	names(y) <- c("X","M")
	list(modelDeterministicXM(y, pars=parameters, ...))
}

# I <- 1.0
# q <- 10
# b <- 0.001 #* 22
# s <- 0.7 #0.748 * 0.5
# h <- 0.15


rs <- getRoot_df(cbind(expand.grid(X=seq(0,15,length.out=15), M=seq(0,1E3,length.out=15)),I=1), pars=c(q=10, b=0.001))#[,c("X","M")]
rs <- rs[complete.cases(rs),,drop=FALSE]
urs <- rs[!duplicated(paste0(round(rs[,"X"],5),round(rs[,"M"],5))), c("I","X","M"), drop=FALSE]

# stabClass <- function(I, pars, Xvals=seq(0,15,length.out=15), Mvals=seq(0,1E3,length.out=15)){
# 	if(missing(pars)){pars <- NULL}
# 	if("I"%in%names(pars)){pars <- pars[names(pars)!="I"]}
# 	gridVals <- cbind(expand.grid(X=Xvals, M=Mvals), I=I)
# 	rs <- getRoot_df(gridVals, pars=pars)
# 	rs <- rs[complete.cases(rs),,drop=FALSE]
# 	urs <- rs[!duplicated(paste0(round(rs[,"X"],4),round(rs[,"M"],4))), c("I","X","M"), drop=FALSE]
#
# 	stab <- function(x){
# 		st <- stability(mDXM, y.star=x[-1,drop=FALSE], parameters=c(I=I,pars), summary=FALSE)
# 		# if(is.null(names(st$parameters))){names(st$parameters) <- "I"}
# 		o <- cbind(data.frame(X=st$y.star[1], M=st$y.star[2], classification=st$classification, Delta=st$Delta, discriminant=st$discriminant, tr=st$tr),as.list(st$parameters))
# 		rownames(o) <- NULL
# 		o$classification <- as.character(o$classification)
# 		return(o)
# 	}
# 	do.call('rbind', apply(urs, 1, stab))
# }

stabClass(I=1)


# stab <- stability(mDXM, y.star=urs[1,-1, drop=FALSE], parameters=c(I=1, q=q, b=b))




# phasePortrait <- function(I, pars, func=mDXM, x.lim=c(0,15), y.lim=c(0,1000), add=FALSE, nFlow=20, nNull=300, t.end=50, inits=matrix(c(0.5,800, 0.1,950, 3,700, 1.5,600, 0.5,100, 2,150, 5,75, 4.5,115, 1,300, 2,400, 4,500, 10,300, 7,200), byrow=TRUE, ncol=2), addLeg=FALSE, legPos="topright"){
# 	if(missing(I)){
# 		if('I' %in% names(pars)){
# 			I <- unname(pars["I"])
# 		}else{
# 			stop("Provide a value for I, P input")
# 		}
# 	}
# 	if(missing(pars)){pars <- NULL}
#
# 		ff <- flowField(mDXM, x.lim=x.lim, y.lim=y.lim, parameters=c(pars), add=add, points=nFlow, xlab="Water P (X)", ylab="Sediment P (M)")
# 		nulls <- nullclines(mDXM, x.lim=x.lim, y.lim=y.lim, parameters=c(pars), add=TRUE, points=nNull)
# 		traj <- trajectory(mDXM, y0=inits, t.end=t.end, parameter=c(pars), colour=rep("black", nrow(inits)), pch=23)
# 		stab <- stabClass(I=I, pars=pars)
# 		stabPCH <- c(
# 			"Saddle" = 4,
# 			"Indeterminate" = 11,
# 			"Stable node" = 19,
# 			"Unstable node" = 21,
# 			"Stable focus" = 10,
# 			"Unstable focus" = 8,
# 			"Centre" = 14
# 		)
# 		points(stab[,c("X","M")], pch=stabPCH[stab[,"classification"]])
#
# 		nd <- !duplicated(stabPCH[stab[,"classification"]])
# 		upch <- stabPCH[stab[,"classification"]][nd]
#
# 		legText <- c("dX/dt=0", "dM/dt=0", "trajectory", names(upch))
# 		cex <- par("cex")
# 		text.font <- NULL
# 		tw <- max(abs(strwidth(legText, units = "user", cex = cex, font = text.font)))*0.9
#
# 		legend(legPos, inset=-c(0.0, 0.0), lty=c(1,1,1, rep(-1, length(upch))), pch=c(NA,NA,23,upch), col=c("red","blue","black",rep('black',length(upch))), legend=legText, ncol=2, bg=adjustcolor('white',0.5), xpd=TRUE, merge=FALSE, x.intersp=0.8, y.intersp=0.5, text.width=tw)
#
# 		return(stab)
# }

dev.new()
par(mfrow=c(2,2), mar=c(2,2,0.5,0.5), mgp=c(1,0.25,0), tcl=-0.15, cex=1, ps=8)
pars <- c(I=0.25, q=8, b=0.001, s=0.7, h=0.15)
phasePortrait(pars=pars, addLeg=TRUE)
pars <- c(I=0.75, q=8, b=0.001, s=0.7, h=0.15)
phasePortrait(pars=pars, addLeg=TRUE)
pars <- c(I=0.98, q=8, b=0.001, s=0.7, h=0.15) # by I=0.9981, the saddle and node have disappeared
phasePortrait(pars=pars, addLeg=TRUE)
pars <- c(I=1, q=8, b=0.001, s=0.7, h=0.15)
phasePortrait(pars=pars, addLeg=TRUE)










# See simLakeP.R for the functions being tested/ played with in this script
# This script was written while writing those functions, 
# in order to test them and as a way of determining what functions would be needed/ useful.


# ---- Find roots and eigenvalues for a range of loading and water P ----
initialValues <- getInit(30, Mrange=c(200, 600), Mn=6)
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


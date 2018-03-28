library(data.table)
library(timeScales)
#phi[j-1] + 0.95*cos(2*pi*dt*j)*(dt*2*pi)
# dt <- 1/160 # period of the sine wave [so not really a 'dt', don't be confused by this]
 

y <- c(0) # observed state
z <- c(0) # true state
n <- 300 # number of time steps
sigma_w <- 0.5 # sd of the process error
sigma_v <- 0.01 # sd of the observation error
sigma_xi <- 0.0001 # sd of the mean parameter's random walk
sigma_eps <- 0.025 # sd of the AR(1) coefficient's random walk
phi0 <- 0.05 # starting value for the AR(1) coefficient
phi <- c(phi0) # will be a vector of AR(1) coefficients
C <- rep(0, n) # vector of the parameter governing the time series mean; will be filled in with non-0 later

# simulate a time-varying AR(1) model
# both the AR(1) coefficient and the mean parameter vary.
# based on Dakos & Ives 2012 Ecosphere
for(j in 2:n){
	phi[j] <- 1*phi[j-1] + rnorm(1, mean=0, sd=sigma_eps) #AR(1) coefficient
	C[j] <- C[j-1] + rnorm(1, mean=0, sd=sigma_xi) # mean parameter
	z[j] <- (z[j-1]-C[j])*phi[j] + C[j] + rnorm(1, mean=0, sd=sigma_w) # true state
	y[j] <- z[j] + rnorm(1, mean=0, sd=sigma_v) # observed state
}
# set vectors of interest to time series objects for convenient plotting
z <- ts(z)
y <- ts(y)
C <- ts(C)
phi <- ts(phi)


# ---- 'dlm' ----
# initial attempt at 'dlm' package
# # this package confused the hell out of me
# # steve said he had a hard time getting it to work too
# library(dlm)
# buildFun <- function(x){
# 	dlmModARMA(ar=x[1], sigma2=exp(x[2]))
# }
# dlmMLE(y, parm=c(0,0), build=buildFun)

# ---- 'tvReg' ----
# initial attempt at 'tvReg' package
# this package basically just does rolling window ac(), but instead of doing a regression of
# x_t vs x_{t-1}, it puts 'weights' on x_{t-1} (weighted least squares) where the weights are determined by the kernel
# library(tvReg) # https://cran.r-project.org/web/packages/tvReg/vignettes/tvReg-vignette.html
# tvmod <- tvAR(z, p=1, type='none', bw=0.1)
# plot(tvmod$tvcoef[,1], phi[-1])
# plot(tvmod)


# ===============
# = Ryan + JAGS =
# ===============
library(R2jags)
Y <- y # will supply JAGS with Y, in case I want to modify it from y (y simulated above)
# Y[sample(10:n, 0.1*(length(Y)))] <- NA # note that JAGS v4.0 and later can now handle missing observations in this type of model; i tested, works great

# ---- Fit Model w/o Mean ----
out_noMean <- tvarss(Y, niter=2E3)
summarize.tvarss(as.tvarss(out_noMean))[]
dev.new()
par(mfrow=c(2,2), mar=c(2.5,2.5,1,1), mgp=c(1.5,0.5,0), tcl=-0.2)
plot(phi, type='l')
plot(C, type='l')
plotPost.tvarss(out_noMean, varName=c('Phi'))

# ---- Fit Model w/ Time-Varying Mean ----
out_mean <- tvarss(Y, niter=2E3, tvMean=TRUE)
summarize.tvarss(as.tvarss(out_noMean))[]
dev.new()
par(mfrow=c(2,2), mar=c(2.5,2.5,1,1), mgp=c(1.5,0.5,0), tcl=-0.2)
plot(phi, type='l')
plot(C, type='l')
plotPost.tvarss(out_mean, varName=c('Phi','C'))




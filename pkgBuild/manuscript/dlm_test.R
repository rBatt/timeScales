# Test fitting a dlm


y <- c(0)
n <- 1000
w <- 0.5
wphi <- 0.025
phi0 <- 0.5
phi <- c(phi0)


for(j in 2:n){
	phi[j] <- phi[j-1] + rnorm(1, mean=0, sd=wphi)
	y[j] <- y[j-1]*phi[j] + rnorm(1, mean=0, sd=w)
}
y <- ts(y)
phi <- ts(phi)

par(mfrow=c(2,2))
plot(y)
plot(phi)

# this package confused the hell out of me
# library(dlm)
# buildFun <- function(x){
# 	dlmModARMA(ar=x[1], sigma2=exp(x[2]))
# }
#
# dlmMLE(y, parm=c(0,0), build=buildFun)


library(tvReg) # https://cran.r-project.org/web/packages/tvReg/vignettes/tvReg-vignette.html
tvmod <- tvAR(y, p=1, type='none', bw=0.1)
plot(tvmod$tvcoef[,1], phi[-1])
plot(tvmod)





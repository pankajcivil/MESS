#===============================================================================
# Preliminary testing with the Naveau (i) model.
# AKA, Papastathopoulos and Tawn (2013) model (iii)
# Note: this is antiquated by the "fit_many_gev_and_naveau_model.R" script.
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================


#===============================================================================
# set up parameters - parametric family (i) from Naveau et al (2016)
lnonstat <- vector('list',3); names(lnonstat) <- c('kappa','sigma','xi')
lnonstat$kappa <- FALSE
lnonstat$sigma <- FALSE
lnonstat$xi    <- FALSE

# set parameter names, with an error check, that the nonstationary parameter
# combinations are supported
if(length(which(unlist(lnonstat)==TRUE))==0) {
  parnames <- c('kappa','sigma','xi')
} else if(length(which(unlist(lnonstat)==TRUE))==1) {
  if(lnonstat$kappa) {parnames <- c('kappa0','kappa1','sigma','xi')}
  else {print('ERROR - if only one nonstationary parameter, it must be kappa')}
} else if(length(which(unlist(lnonstat)==TRUE))==2) {
  if(lnonstat$kappa & lnonstat$sigma) {parnames <- c('kappa0','kappa1','sigma0','sigma1','xi')}
  else {print('ERROR - if only two nonstationary parameters, they must be kappa and scale')}
} else if(length(which(unlist(lnonstat)==TRUE))==3) {
  parnames <- c('kapp0','mu1','sigma0','sigma1','xi0','xi1')
} else {print('ERROR - unknown number of non-stationary parameters')}
#===============================================================================


#===============================================================================
# read temperature hindcast (1880-2016) and forecast (2016-2100) data
# yields time_forc, temperature_forc (temperature relative to 1901-2000, deg C)

library(ncdf4)
source('read_data_temperature.R')
#===============================================================================


#===============================================================================
# read and process the tide gauge data from Delfzijl

source('read_data_tidegauge.R')

data_calib <- lsl.max
auxiliary <- temperature_forc[1:length(data_calib)]
#===============================================================================


#===============================================================================
# set up some priors

# uniform priors
priors.unif <- vector('list', length(parnames)); names(priors.unif) <- parnames
for (p in 1:length(parnames)) {
  priors.unif[[p]] <- vector('list', 2); names(priors.unif[[p]]) <- c('lb','ub')
}
#===============================================================================


#===============================================================================

# TODO
# TODO -- check these
# TODO

naveau_level <- function(p, kappa, sigma, xi){
  level <- NULL
  level <- (sigma/xi)*( (1-(p^(1/kappa)))^(-xi) -1)
  return(level)
}

naveau_pdf <- function(x, kappa, sigma, xi){
  prob <- NULL
  prob <- kappa*((1-(1+xi*x/sigma)^(-1/xi))^(kappa-1))*((1+xi*x/sigma)^(-(1+1/xi)))/sigma
  return(prob)
}

naveau_cdf <- function(x, kappa, sigma, xi){
  q <- NULL
  q <- (1-(1+xi*x/sigma)^(-1/xi))^kappa
  return(q)
}

# inverse cdf for generalized pareto distribution:
gpinv <- function(x, mu, sigma, xi){output <- mu + (sigma/xi)*(((1-x)^(-xi))-1); return(output)}
# verified (9 June 2017) that 'extRemes' package function qevd(..., type='GP')
# indeed reproduces this. Need for Naveau likelihood function (their Appendix)

# good testing, "calibrate" by eyeball metric to see how the parameters interact
tmp <- naveau_pdf(x, 20000, 5, 0.3); hist(lsl.max, xlim=c(0,1000)); plot(x, tmp, type='l', xlim=c(0,1000))

# verify that naveau_pdf(...) is coded correctly - against Fig 1 from Naveau et al:
xtmp <- seq(from=0, to=7, by=0.02)
ptmp1 <- naveau_pdf(x=xtmp, 1, 1, 0.5)
ptmp2 <- naveau_pdf(x=xtmp, 2, 1, 0.5)
ptmp3 <- naveau_pdf(x=xtmp, 5, 1, 0.5)
plot(xtmp, ptmp1, type='l', xlim=c(0,7), ylim=c(0,1)); lines(xtmp, ptmp2, type='l', lty=2); lines(xtmp, ptmp3, type='l', lty=3)
#===============================================================================


#===============================================================================
# fit MLE Naveau-(i) forms

library(extRemes)
library(adaptMCMC)
library(lhs)
library(DEoptim)

source('likelihood_naveau.R')

# preliminary latin hypercube

# ALL STATIONARY
parnames <- c('kappa','sigma','xi')
bound.lower <- c(0, 0, 0)
bound.upper <- c(50000, 100, 10)
auxiliary <- NULL

# KAPPA NONSTATIONARY
parnames <- c('kappa0','kappa1','sigma','xi')
bound.lower <- c(0, -10, 0, 0)
bound.upper <- c(50000, 10, 40, 5)
auxiliary <- trimmed_forcing(year.unique, time_forc, temperature_forc)$temperature

niter.lhs <- 1e5
nparam.lhs <- length(bound.lower)




parameters.lhs0 <- randomLHS(n=niter.lhs, k=nparam.lhs)
parameters.lhs  <- parameters.lhs0
for (p in 1:nparam.lhs) {
  parameters.lhs[,p] <- parameters.lhs0[,p]*(bound.upper[p]-bound.lower[p]) + bound.lower[p]
}

llik.lhs <- rep(NA, niter.lhs)
pb <- txtProgressBar(min=0,max=niter.lhs,initial=0,style=3)
for (i in 1:niter.lhs) {
  llik.lhs[i] <- log_like_naveau(parameters=parameters.lhs[i,], parnames=parnames, data_calib=data_calib, auxiliary=auxiliary)
  setTxtProgressBar(pb, i)
}
close(pb)
lhs.sample <- data.frame(cbind(parameters.lhs, llik.lhs))
colnames(lhs.sample) <- c(parnames,'llike')

# sort lhs results from highest likelihood to lowest
lhs.sample <- lhs.sample[rev(order(llik.lhs)),]

# Differential evolution optimization
# (can only minimize, so use negative log-likelihood function)

# calculate BIC with each of the maximum likelihood sets (not that 'bestval' is
# from the negative log-likelihood, cancelling out the - in BIC)

NP.deoptim <- 200
niter.deoptim <- 200
F.deoptim <- 0.8
CR.deoptim <- 0.9


#========================
# all stationary
parnames <- c('kappa','sigma','xi')
bound.lower <- c(0, 0, -10)
bound.upper <- c(2e7, 1000, 10)
deoptim.nav3 <- DEoptim(neg_log_like_naveau, lower=bound.lower, upper=bound.upper,
                        DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                        parnames=parnames, data_calib=data_calib)
names(deoptim.nav3$optim$bestmem) <- parnames
bic.nav3 <- 2*deoptim.nav3$optim$bestval + length(parnames)*log(length(data_calib))


# plot this survival function against the empirical one
parameters <- deoptim.nav3$optim$bestmem
x <- seq(0, 10000, 10)
nav.cdf <- naveau_cdf(x=x, kappa=parameters[1], sigma=parameters[2], xi=parameters[3])
plot(esf.levels, log10(esf.values), ylim=c(-2.5,0), xlim=c(0,6000), xlab='Level [cm]', ylab='Survival function [1-cdf]', yaxt='n')
axis(2, at=seq(-3, 0, by=1), label=parse(text=paste("10^", seq(-3,0), sep="")))
lines(x, log10(1-nav.cdf), col='red')

# plot with the fitted GEV
gev.cdf <- pevd(q=x, loc=deoptim.gev3$optim$bestmem[1], scale=deoptim.gev3$optim$bestmem[2], shape=deoptim.gev3$optim$bestmem[3])
plot(x, log10(1-gev.cdf), type='l', ylim=c(-2.5,0), xlim=c(0,6000), xlab='Level [cm]', ylab='Survival function [1-cdf]', yaxt='n')
axis(2, at=seq(-3, 0, by=1), label=parse(text=paste("10^", seq(-3,0), sep="")))
lines(x, log10(1-nav.cdf), col='red')
points(esf.levels, log10(esf.values))
#========================


#========================
# lower tail parameter nonstationary
parnames <- c('kappa0','kappa1','sigma','xi')
bound.lower <- c(0, -100, 0, -10)
bound.upper <- c(2e7, 100, 1000, 10)
auxiliary <- trimmed_forcing(year.unique, time_forc, temperature_forc)$temperature
deoptim.nav4 <- DEoptim(neg_log_like_naveau, lower=bound.lower, upper=bound.upper,
                        DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                        parnames=parnames, data_calib=data_calib, auxiliary=auxiliary)
names(deoptim.nav4$optim$bestmem) <- parnames
bic.nav4 <- 2*deoptim.nav4$optim$bestval + length(parnames)*log(length(data_calib))

# plot this survival function against the empirical one
kappa0 <- deoptim.nav4$optim$bestmem[1]; kappa1 <- deoptim.nav4$optim$bestmem[2]
sigma <- deoptim.nav4$optim$bestmem[3]; xi <- deoptim.nav4$optim$bestmem[4]
kappa <- kappa0 + kappa1*auxiliary; sigma <- rep(sigma, length(auxiliary)); xi <- rep(xi, length(auxiliary))
cdf.tmp <- naveau_cdf(x=lsl.max, kappa=kappa, sigma=sigma, xi=xi)
plot(esf.levels, log10(esf.values))
#lines(x, log10(1-cdf.tmp), col='red')
points(lsl.max, log10(1-cdf.tmp), pch=16, col='blue')
#========================


#========================
# lower tail and scale nonstationary
parnames <- c('kappa0','kappa1','sigma0','sigma1','xi')
bound.lower <- c(0, -100, 0, -100, -10)
bound.upper <- c(2e7, 100, 100, 100, 10)
auxiliary <- trimmed_forcing(year.unique, time_forc, temperature_forc)$temperature
deoptim.nav5 <- DEoptim(neg_log_like_naveau, lower=bound.lower, upper=bound.upper,
                        DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                        parnames=parnames, data_calib=data_calib, auxiliary=auxiliary)
names(deoptim.nav5$optim$bestmem) <- parnames
bic.nav5 <- 2*deoptim.nav5$optim$bestval + length(parnames)*log(length(data_calib))

# plot this survival function against the empirical one
kappa0 <- deoptim.nav5$optim$bestmem[1]; kappa1 <- deoptim.nav5$optim$bestmem[2]
sigma0 <- deoptim.nav5$optim$bestmem[3]; sigma1 <- deoptim.nav5$optim$bestmem[4]
xi <- deoptim.nav5$optim$bestmem[5]
kappa <- kappa0 + kappa1*auxiliary; sigma <- exp(sigma0 + sigma1*auxiliary)
xi <- rep(xi, length(auxiliary))
cdf.tmp <- naveau_cdf(x=lsl.max, kappa=kappa, sigma=sigma, xi=xi)
plot(esf.levels, log10(esf.values))
#lines(x, log10(1-cdf.tmp), col='red')
points(lsl.max, log10(1-cdf.tmp), pch=16, col='blue')
#========================


#========================
# location, scale and shape all nonstationary
parnames <- c('kappa0','kappa1','sigma0','sigma1','xi0','xi1')
bound.lower <- c(0, -100, 0, -100, -10, -10)
bound.upper <- c(2e7, 100, 100, 100, 10, 10)
auxiliary <- trimmed_forcing(year.unique, time_forc, temperature_forc)$temperature
deoptim.nav6 <- DEoptim(neg_log_like_naveau, lower=bound.lower, upper=bound.upper,
                        DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                        parnames=parnames, data_calib=data_calib, auxiliary=auxiliary)
names(deoptim.nav6$optim$bestmem) <- parnames
bic.nav6 <- 2*deoptim.nav6$optim$bestval + length(parnames)*log(length(data_calib))

# plot this survival function against the empirical one
kappa0 <- deoptim.nav6$optim$bestmem[1]; kappa1 <- deoptim.nav6$optim$bestmem[2]
sigma0 <- deoptim.nav6$optim$bestmem[3]; sigma1 <- deoptim.nav6$optim$bestmem[4]
xi0 <- deoptim.nav6$optim$bestmem[5]; xi1 <- deoptim.nav6$optim$bestmem[6]
kappa <- kappa0 + kappa1*auxiliary; sigma <- exp(sigma0 + sigma1*auxiliary)
xi <- xi0 + xi1*auxiliary;
cdf.tmp <- naveau_cdf(x=lsl.max, kappa=kappa, sigma=sigma, xi=xi)
plot(esf.levels, log10(esf.values))
#lines(x, log10(1-cdf.tmp), col='red')
points(lsl.max, log10(1-cdf.tmp), pch=16, col='blue')
#========================





#===============================================================================
# End
#===============================================================================

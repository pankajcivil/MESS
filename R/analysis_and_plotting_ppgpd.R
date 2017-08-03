#===============================================================================
# analysis_and_plotting.R
#
# Analysis and plotting for wong, klufas and keller pp/gpd analysis.
# This assumes you've downloaded the Github repo and are running from the
# 'R' directory.
#
# If any this is vague/confusing and you have questions, email me. Seriously.
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================


#===============================================================================
# Some preliminaries

setwd('~/codes/EVT/R')

library(ncdf4)
library(extRemes)

# get some colorblind-friendly colors to plot with
source('colorblindPalette.R')

# set useful directories -- assumes you are in the 'R' directory within the repo
# directory structure
plot.dir <- '../figures/'
output.dir <- '../output/'

# calibrated parameter sets (samples; all should be same size)
# these are the results from 'calibration_dayPOT-experiments_driver.R'
filename.norfolk.normalgamma <-  paste(output.dir,'calibratedParameters_ppgpd-experiments_norfolk_normalgamma_27Jul2017.nc', sep='')
filename.norfolk.uniform <-      paste(output.dir,'calibratedParameters_ppgpd-experiments_norfolk_uniform_27Jul2017.nc',     sep='')
filename.balboa.normalgamma <-   paste(output.dir,'calibratedParameters_ppgpd-experiments_balboa_normalgamma_28Jul2017.nc',  sep='')
filename.balboa.uniform <-       paste(output.dir,'calibratedParameters_ppgpd-experiments_balboa_uniform_28Jul2017.nc',      sep='')
filename.delfzijl.normalgamma <- paste(output.dir,'calibratedParameters_ppgpd-experiments_delfzijl_normalgamma_29Jul2017.nc',sep='')
filename.delfzijl.uniform <-     paste(output.dir,'calibratedParameters_ppgpd-experiments_delfzijl_uniform_29Jul2017.nc',    sep='')

filename.norfolk.data <- '../data/tidegauge_processed_norfolk_26Jul2017.rds'
filename.balboa.data <- '../data/tidegauge_processed_balboa_26Jul2017.rds'
filename.delfzijl.data <- '../data/tidegauge_processed_delfzijl_26Jul2017.rds'

filename.priors <- '../output/surge_priors_normalgamma_ppgpd_26Jul2017.rds'

# file to save progress as you run
filename.saveprogress <- '../output/analysis_inprogress.RData'

#===============================================================================



#===============================================================================
#===============================================================================
# ANALYSIS
#===============================================================================
#===============================================================================




#===============================================================================
# read temperature data for projecting surge levels covarying with global mean
# surface temperature
#===============================================================================

# yields 'temperature_forc' and 'time_forc', between 1850 and 2100
# temperatures are relative to 1901-2000 mean temperature.

source('read_data_temperature.R')

# check the normalization
temperature.check <- mean(temperature_forc[which(time_forc==1901):which(time_forc==2000)])
print(paste('mean(temperature_forc between 1901 and 2000) is: ', temperature.check, sep=''))
print('(should be normalized to 0 during that period)')

#===============================================================================
# read parameter sets
#===============================================================================

# create list objects to store the parameters
site.names <- c('Balboa','Delfzijl','Norfolk'); nsites <- length(site.names)
data.lengths <- c(30,50,70,90,110,130); ndata.exp <- length(data.lengths)
data.experiment.names <- rep(NA, length(data.lengths)); for (dd in 1:length(data.lengths)) {data.experiment.names[dd] <- paste('y',data.lengths[dd],sep='')}
types.of.gpd <- c('gpd3','gpd4','gpd5','gpd6'); nmodel <- length(types.of.gpd)

# define a generic list object that will have values for (i) each site,
# (ii) each data length experiment and (iii) each GPD model structure. (in that
# order)
list.init <- vector('list', nsites); names(list.init) <- site.names
for (site in site.names) {
  list.init[[site]] <- vector('list', length(data.lengths)); names(list.init[[site]]) <- data.experiment.names
  for (data.exp in data.experiment.names) {
    list.init[[site]][[data.exp]] <- vector('list', nmodel); names(list.init[[site]][[data.exp]]) <- types.of.gpd
  }
}

# initialize anything you would like to be defined for each site, for each data
# length and for each GPD model structure
gpd.parameters <- list.init
rl100 <- list.init

# index to store within the list level [[data.experiments]] which of the
# experiments corresponds to using all of the data. so you would use:
#    gpd.parameters$delfzijl[[all.data]]$gpd3
# for example to get at the calibrated gpd3 (stationary model) parameters for
# Delfzijl.
all.data <- rep(NA, nsites); names(all.data) <- site.names
data.lengths <- vector('list', nsites); names(data.lengths) <- site.names
data.sites <- vector('list', nsites); names(data.sites) <- site.names

# hard-coding here, even though it is bad practice. the different sites ahve
# different lengths of record, making things a bit complicated.

# parameter nmaes
parnames <- vector('list', nmodel); names(parnames) <- types.of.gpd

for (site in site.names) {
  if(site=='Norfolk') {
    ncdata <- nc_open(filename.norfolk.normalgamma)
    for (model in types.of.gpd) {
      gpd.parameters[[site]]$y30[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd30.',model,sep='')))
      gpd.parameters[[site]]$y50[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd50.',model,sep='')))
      gpd.parameters[[site]]$y70[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd70.',model,sep='')))
      gpd.parameters[[site]]$y90[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd89.',model,sep='')))
      parnames[[model]] <- ncvar_get(ncdata, paste('parnames.',model,sep=''))
    }
    n.ensemble <- nrow(gpd.parameters[[site]]$y30$gpd3)
    data.sites[[site]] <- readRDS(filename.norfolk.data)
    data.lengths[[site]] <- names(data.sites[[site]])[intersect(which(nchar(names(data.sites[[site]]))>3) , grep('gpd', names(data.sites[[site]])))]
    all.data[[site]] <- length(data.lengths[[site]])
  } else if(site=='Balboa') {ncdata <- nc_open(filename.balboa.normalgamma)
    for (model in types.of.gpd) {
      gpd.parameters[[site]]$y30[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd30.',model,sep='')))
      gpd.parameters[[site]]$y50[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd50.',model,sep='')))
      gpd.parameters[[site]]$y70[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd70.',model,sep='')))
      gpd.parameters[[site]]$y90[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd90.',model,sep='')))
      gpd.parameters[[site]]$y110[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd107.',model,sep='')))
    }
    data.sites[[site]] <- readRDS(filename.balboa.data)
    data.lengths[[site]] <- names(data.sites[[site]])[intersect(which(nchar(names(data.sites[[site]]))>3) , grep('gpd', names(data.sites[[site]])))]
    all.data[[site]] <- length(data.lengths[[site]])
    if(nrow(gpd.parameters[[site]]$y30$gpd3) != n.ensemble) {print('ERROR - all sites ensembles must be the same size')}
  } else if(site=='Delfzijl') {ncdata <- nc_open(filename.delfzijl.normalgamma)
    for (model in types.of.gpd) {
      gpd.parameters[[site]]$y30[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd30.',model,sep='')))
      gpd.parameters[[site]]$y50[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd50.',model,sep='')))
      gpd.parameters[[site]]$y70[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd70.',model,sep='')))
      gpd.parameters[[site]]$y90[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd90.',model,sep='')))
      gpd.parameters[[site]]$y110[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd110.',model,sep='')))
      gpd.parameters[[site]]$y130[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd137.',model,sep='')))
    }
    data.sites[[site]] <- readRDS(filename.delfzijl.data)
    data.lengths[[site]] <- names(data.sites[[site]])[intersect(which(nchar(names(data.sites[[site]]))>3) , grep('gpd', names(data.sites[[site]])))]
    all.data[[site]] <- length(data.lengths[[site]])
    if(nrow(gpd.parameters[[site]]$y30$gpd3) != n.ensemble) {print('ERROR - all sites ensembles must be the same size')}
  } else {print('ERROR - unrecognized site name')}
}

#===============================================================================
# calculate return levels
#===============================================================================

# years to grab the return levels
rl.years <- c(2000, 2016, 2065); nyears <- length(rl.years)
year.names <- rep(NA, length(rl.years)); for (year in 1:length(rl.years)) {year.names[year] <- paste('y',rl.years[year],sep='')}
temperature.years <- temperature_forc[which(time_forc == rl.years)]; names(temperature.years) <- year.names
print('you probably just got an error message - dont worry about it, probably')

# this is definitely coded like a neanderthal. can come back and code it more
# efficiently..?
for (site in site.names) {
  for (data.len in data.experiment.names[1:all.data[[site]]]) {
    for (model in types.of.gpd) {
      rl100[[site]][[data.len]][[model]] <- vector('list', nyears); names(rl100[[site]][[data.len]][[model]]) <- year.names
      for (year in year.names) {
        rl100[[site]][[data.len]][[model]][[year]] <- rep(NA, n.ensemble)
        print(paste('calculating return levels...',site,data.len,model,year,sep=' - '))
        tbeg <- proc.time()
        pb <- txtProgressBar(min=0,max=n.ensemble, initial=0,style=3)
        for (sow in 1:n.ensemble) {
          if (length(parnames[[model]])==3) {
            lambda <- gpd.parameters[[site]][[data.len]][[model]][sow,match('lambda', parnames[[model]])]
            sigma <- gpd.parameters[[site]][[data.len]][[model]][sow,match('sigma', parnames[[model]])]
            xi <- gpd.parameters[[site]][[data.len]][[model]][sow,match('xi', parnames[[model]])]
          } else if(length(parnames[[model]])==4) {
            lambda <- gpd.parameters[[site]][[data.len]][[model]][sow,match('lambda0', parnames[[model]])] +
                      gpd.parameters[[site]][[data.len]][[model]][sow,match('lambda1', parnames[[model]])]*temperature.years[[year]]
            sigma <- gpd.parameters[[site]][[data.len]][[model]][sow,match('sigma', parnames[[model]])]
            xi <- gpd.parameters[[site]][[data.len]][[model]][sow,match('xi', parnames[[model]])]
          } else if(length(parnames[[model]])==5) {
            lambda <- gpd.parameters[[site]][[data.len]][[model]][sow,match('lambda0', parnames[[model]])] +
                      gpd.parameters[[site]][[data.len]][[model]][sow,match('lambda1', parnames[[model]])]*temperature.years[[year]]
            sigma <- exp(gpd.parameters[[site]][[data.len]][[model]][sow,match('sigma0', parnames[[model]])] +
                         gpd.parameters[[site]][[data.len]][[model]][sow,match('sigma1', parnames[[model]])]*temperature.years[[year]])
            xi <- gpd.parameters[[site]][[data.len]][[model]][sow,match('xi', parnames[[model]])]
          } else if(length(parnames[[model]])==6) {
            lambda <- gpd.parameters[[site]][[data.len]][[model]][sow,match('lambda0', parnames[[model]])] +
                      gpd.parameters[[site]][[data.len]][[model]][sow,match('lambda1', parnames[[model]])]*temperature.years[[year]]
            sigma <- exp(gpd.parameters[[site]][[data.len]][[model]][sow,match('sigma0', parnames[[model]])] +
                         gpd.parameters[[site]][[data.len]][[model]][sow,match('sigma1', parnames[[model]])]*temperature.years[[year]])
            xi <- gpd.parameters[[site]][[data.len]][[model]][sow,match('xi0', parnames[[model]])] +
                  gpd.parameters[[site]][[data.len]][[model]][sow,match('xi1', parnames[[model]])]*temperature.years[[year]]
          }
          rl100[[site]][[data.len]][[model]][[year]][sow] <- rlevd(100, scale=sigma, shape=xi,
                                                                   threshold=data.sites[[site]]$gpd$threshold,
                                                                   type='GP',
                                                                   npy=365.25,
                                                                   rate=lambda)
          setTxtProgressBar(pb, sow)
        }
        close(pb)
      }
    }
  }
  save.image(file=filename.saveprogress)
}

#===============================================================================
# Bayesian model averaging ensemble
#===============================================================================

# need the likelihood functions and the priors
source('likelihood_ppgpd.R')
priors <- readRDS(filename.priors)

llik <- list.init
lpri <- list.init
lpost <- list.init
llik.mod <- list.init
lpost.mod <- list.init
llik.ref <- vector('list', nsites); names(llik.ref) <- site.names
lpost.ref <- vector('list', nsites); names(lpost.ref) <- site.names
for (site in site.names) {
  llik.ref[[site]] <- vector('list', ndata.exp); names(llik.ref[[site]]) <- data.experiment.names
  lpost.ref[[site]] <- vector('list', ndata.exp); names(lpost.ref[[site]]) <- data.experiment.names
}

for (site in site.names) {
  for (ind.data in 1:all.data[[site]]) {
    for (model in types.of.gpd) {
      print(paste('calculating log-posts/-likelihoods/-priors ',site,data.experiment.names[ind.data],model, sep=' - '))
      if(model=='gpd3') {auxiliary <- NULL
      } else {auxiliary <- trimmed_forcing(data.sites[[site]][[data.lengths[[site]][ind.data]]]$year, time_forc, temperature_forc)$temperature}
      lpri[[site]][[ind.data]][[model]] <- sapply(1:n.ensemble, function(sow) {log_prior_ppgpd(parameters=gpd.parameters[[site]][[ind.data]][[model]][sow,], parnames=parnames[[model]], priors=priors, model=model)})
      llik[[site]][[ind.data]][[model]] <- sapply(1:n.ensemble, function(sow) {log_like_ppgpd(parameters=gpd.parameters[[site]][[ind.data]][[model]][sow,], parnames=parnames[[model]], data_calib=data.sites[[site]][[data.lengths[[site]][ind.data]]], auxiliary=auxiliary)})
      lpost[[site]][[ind.data]][[model]] <- lpri[[site]][[ind.data]][[model]] + llik[[site]][[ind.data]][[model]]
      llik.mod[[site]][[ind.data]][[model]] <- mean(llik[[site]][[ind.data]][[model]][is.finite(llik[[site]][[ind.data]][[model]])])
      lpost.mod[[site]][[ind.data]][[model]] <- mean(lpost[[site]][[ind.data]][[model]][is.finite(lpost[[site]][[ind.data]][[model]])])
    }
    # reference likelihood and posterior values for BMA weight calculation
    llik.ref[[site]][[ind.data]] <- median(unlist(llik.mod[[site]][[ind.data]]))
    lpost.ref[[site]][[ind.data]] <- median(unlist(lpost.mod[[site]][[ind.data]]))
  }
}
save.image(file=filename.saveprogress)

# create the BMA-weighted ensemble

# throw this to HPC; need the RData files and the raw parameter sets

TODO <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< HERE NOW!!
TODO <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< HERE NOW!!


#===============================================================================
# Calculate AIC, BIC, DIC for each ensemble
# (of course, this could be wrapped up in the loop structure above, but it is
# nice to compartmentalize the code)
#===============================================================================

# new initialization object for the model comparison metrics, so we have a table
# comparing each of the models instead of a list for each of them
metric.init <- vector('list', nsites); names(metric.init) <- site.names
for (site in site.names) {
  metric.init[[site]] <- mat.or.vec(all.data[[site]], nmodel)
  rownames(metric.init[[site]]) <- data.experiment.names[1:all.data[[site]]]
  colnames(metric.init[[site]]) <- types.of.gpd
}

aic <- metric.init
bic <- metric.init
dic <- metric.init
  dev <- list.init      # deviance
  par.mean <- list.init # expected value of the model parameters
  dev.mean <- list.init # deviance at expected value of the parameters
  np.eff <- list.init   # effective number of parameters (p_D)
bma.weight <- metric.init

for (site in site.names) {
  for (ind.data in 1:all.data[[site]]) {
    ndata <- sum(data.sites[[site]][[data.lengths[[site]][ind.data]]]$counts)
    for (model in types.of.gpd) {
      imax <- which.max(llik[[site]][[ind.data]][[model]])
      aic[[site]][ind.data, model] <- -2*llik[[site]][[ind.data]][[model]][imax] + length(parnames[[model]])*2
      bic[[site]][ind.data, model] <- -2*llik[[site]][[ind.data]][[model]][imax] + length(parnames[[model]])*log(ndata)
      dev[[site]][[ind.data]][[model]] <- -2*llik[[site]][[ind.data]][[model]]
      if(model=='gpd3') {auxiliary <- NULL
      } else {auxiliary <- trimmed_forcing(data.sites[[site]][[data.lengths[[site]][ind.data]]]$year, time_forc, temperature_forc)$temperature}
      par.mean[[site]][[ind.data]][[model]] <- apply(gpd.parameters[[site]][[ind.data]][[model]], MARGIN=2, FUN=mean)
      dev.mean[[site]][[ind.data]][[model]] <- log_like_ppgpd(parameters=par.mean[[site]][[ind.data]][[model]], parnames=parnames[[model]], data_calib=data.sites[[site]][[data.lengths[[site]][ind.data]]], auxiliary=auxiliary)
      np.eff[[site]][[ind.data]][[model]] <- mean(dev[[site]][[ind.data]][[model]]) - dev.mean[[site]][[ind.data]][[model]]
      dic[[site]][ind.data, model] <- np.eff[[site]][[ind.data]][[model]] + mean(dev[[site]][[ind.data]][[model]])
      bma.weight[[site]][ind.data, model] <- exp(lpost.mod[[site]][[ind.data]][[model]] - lpost.ref[[site]][[ind.data]])/sum(exp(unlist(lpost.mod[[site]][[ind.data]]) - lpost.ref[[site]][[ind.data]]))
    }
  }
}

TODO <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< HERE NOW!!
TODO <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< HERE NOW!!



#===============================================================================




#===============================================================================
#===============================================================================
# PLOTTING
#===============================================================================
#===============================================================================


#
#===============================================================================
# SOM FIGURE S1 -- show the histograms of the MLE parameters and superimpose the
#                  prior distribution (normal or gamma)

load('../output/fitting_priors.RData')

nbins <- 12 # note that there are only 30 sites...
frac.ran <- 0.35 # extend axes above/below max/min range

# get standard limits for each of the 6 parameters
lims <- mat.or.vec(6,2)
pp <- 1 # lambda0
  parameters.pooled <- c(deoptim.all$gpd3[,1], deoptim.all$gpd4[,1], deoptim.all$gpd5[,1], deoptim.all$gpd6[,1])
  ran <- diff(quantile(parameters.pooled, c(0,1)))
  lims[pp,] <- c(min(parameters.pooled) - frac.ran*ran , max(parameters.pooled) + frac.ran*ran)
pp <- 2 # lambda1
  parameters.pooled <- c(deoptim.all$gpd4[,2], deoptim.all$gpd5[,2], deoptim.all$gpd6[,2])
  ran <- diff(quantile(parameters.pooled, c(0,1)))
  lims[pp,] <- c(min(parameters.pooled) - frac.ran*ran , max(parameters.pooled) + frac.ran*ran)
pp <- 3 # sigma0
  parameters.pooled <- c(log(deoptim.all$gpd3[,2]), log(deoptim.all$gpd4[,3]), deoptim.all$gpd5[,3], deoptim.all$gpd6[,3])
  ran <- diff(quantile(parameters.pooled, c(0,1)))
  lims[pp,] <- c(min(parameters.pooled) - frac.ran*ran , max(parameters.pooled) + frac.ran*ran)
pp <- 4 # sigma1
  parameters.pooled <- c(deoptim.all$gpd5[,4], deoptim.all$gpd6[,4])
  ran <- diff(quantile(parameters.pooled, c(0,1)))
  lims[pp,] <- c(min(parameters.pooled) - frac.ran*ran , max(parameters.pooled) + frac.ran*ran)
pp <- 5 # xi0
  parameters.pooled <- c(deoptim.all$gpd5[,5], deoptim.all$gpd6[,5])
  ran <- diff(quantile(parameters.pooled, c(0,1)))
  lims[pp,] <- c(min(parameters.pooled) - frac.ran*ran , max(parameters.pooled) + frac.ran*ran)
pp <- 6 # xi1
  parameters.pooled <- c(deoptim.all$gpd6[,6])
  ran <- diff(quantile(parameters.pooled, c(0,1)))
  lims[pp,] <- c(min(parameters.pooled) - frac.ran*ran , max(parameters.pooled) + frac.ran*ran)



pdf(paste(plot.dir,'priors_normalgamma.pdf',sep=''), height=7, width=10, colormodel='cmyk')
par(mfrow=c(4,6))
##=============================
model <- 'gpd3'
pp <- 1; par(mai=c(.25,.59,.25,.01))
box.width <- diff(lims[1,])/nbins; box.edges <- seq(from=lims[1,1], to=lims[1,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[1,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[1,1], to=lims[1,2], length.out=1000); lines(x.tmp, dgamma(x=x.tmp, shape=priors_normalgamma[[model]]$lambda$shape, rate=priors_normalgamma[[model]]$lambda$rate), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
mtext('Probability density', side=2, line=1.2, cex=1)
mtext('ST', side=2, line=3, cex=1)
plot.new()
pp <- 2; par(mai=c(.25,.35,.25,.25))
box.width <- diff(lims[3,])/nbins; box.edges <- seq(from=lims[3,1], to=lims[3,2], by=box.width)
hist(log(deoptim.all[[model]][,pp]), xlim=lims[3,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
rate.tmp <- median(log(deoptim.all[[model]][,2])) / (0.5*(max(log(deoptim.all[[model]][,2]))-min(log(deoptim.all[[model]][,2]))))^2
shape.tmp <- median(log(deoptim.all[[model]][,2])) * rate.tmp
x.tmp <- seq(from=lims[3,1], to=lims[3,2], length.out=1000); lines(x.tmp, dgamma(x=x.tmp, shape=shape.tmp, rate=rate.tmp), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
plot.new()
pp <- 3; par(mai=c(.25,.3,.25,.3))
box.width <- diff(lims[5,])/nbins; box.edges <- seq(from=lims[5,1], to=lims[5,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[5,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[5,1], to=lims[5,2], length.out=1000); lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[model]]$xi$mean, sd=priors_normalgamma[[model]]$xi$sd), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
plot.new()
##=============================
model <- 'gpd4'
pp <- 1; par(mai=c(.25,.59,.25,.01))
box.width <- diff(lims[1,])/nbins; box.edges <- seq(from=lims[1,1], to=lims[1,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[1,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[1,1], to=lims[1,2], length.out=1000); lines(x.tmp, dgamma(x=x.tmp, shape=priors_normalgamma[[model]]$lambda0$shape, rate=priors_normalgamma[[model]]$lambda0$rate), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
mtext('Probability density', side=2, line=1.2, cex=1)
mtext('NS1', side=2, line=3, cex=1)
pp <- 2; par(mai=c(.25,.35,.25,.25))
box.width <- diff(lims[2,])/nbins; box.edges <- seq(from=lims[2,1], to=lims[2,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[2,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[2,1], to=lims[2,2], length.out=1000); lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[model]]$lambda1$mean, sd=priors_normalgamma[[model]]$lambda1$sd), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
pp <- 3; par(mai=c(.25,.3,.25,.3))
box.width <- diff(lims[3,])/nbins; box.edges <- seq(from=lims[3,1], to=lims[3,2], by=box.width)
hist(log(deoptim.all[[model]][,pp]), xlim=lims[3,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
rate.tmp <- median(log(deoptim.all[[model]][,3])) / (0.5*(max(log(deoptim.all[[model]][,3]))-min(log(deoptim.all[[model]][,3]))))^2
shape.tmp <- median(log(deoptim.all[[model]][,3])) * rate.tmp
x.tmp <- seq(from=lims[3,1], to=lims[3,2], length.out=1000); lines(x.tmp, dgamma(x=x.tmp, shape=shape.tmp, rate=rate.tmp), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
plot.new()
pp <- 4; par(mai=c(.25,.3,.25,.3))
box.width <- diff(lims[5,])/nbins; box.edges <- seq(from=lims[5,1], to=lims[5,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[5,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[5,1], to=lims[5,2], length.out=1000); lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[model]]$xi$mean, sd=priors_normalgamma[[model]]$xi$sd), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
plot.new()
##=============================
model <- 'gpd5'
pp <- 1; par(mai=c(.25,.59,.25,.01))
box.width <- diff(lims[1,])/nbins; box.edges <- seq(from=lims[1,1], to=lims[1,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[1,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[1,1], to=lims[1,2], length.out=1000); lines(x.tmp, dgamma(x=x.tmp, shape=priors_normalgamma[[model]]$lambda0$shape, rate=priors_normalgamma[[model]]$lambda0$rate), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
mtext('Probability density', side=2, line=1.2, cex=1)
mtext('NS2', side=2, line=3, cex=1)
pp <- 2; par(mai=c(.25,.35,.25,.25))
box.width <- diff(lims[2,])/nbins; box.edges <- seq(from=lims[2,1], to=lims[2,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[2,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[2,1], to=lims[2,2], length.out=1000); lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[model]]$lambda1$mean, sd=priors_normalgamma[[model]]$lambda1$sd), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
pp <- 3; par(mai=c(.25,.3,.25,.3))
box.width <- diff(lims[3,])/nbins; box.edges <- seq(from=lims[3,1], to=lims[3,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[3,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[3,1], to=lims[3,2], length.out=1000); lines(x.tmp, dgamma(x=x.tmp, shape=priors_normalgamma[[model]]$sigma0$shape, rate=priors_normalgamma[[model]]$sigma0$rate), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
pp <- 4; par(mai=c(.25,.3,.25,.3))
box.width <- diff(lims[4,])/nbins; box.edges <- seq(from=lims[4,1], to=lims[4,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[4,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[4,1], to=lims[4,2], length.out=1000); lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[model]]$sigma1$mean, sd=priors_normalgamma[[model]]$sigma1$sd), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
pp <- 5; par(mai=c(.25,.3,.25,.3))
box.width <- diff(lims[5,])/nbins; box.edges <- seq(from=lims[5,1], to=lims[5,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[5,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[5,1], to=lims[5,2], length.out=1000); lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[model]]$xi$mean, sd=priors_normalgamma[[model]]$xi$sd), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
plot.new()
##=============================
model <- 'gpd6'
pp <- 1; par(mai=c(.5,.59,.01,.01))
box.width <- diff(lims[1,])/nbins; box.edges <- seq(from=lims[1,1], to=lims[1,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[1,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[1,1], to=lims[1,2], length.out=1000); lines(x.tmp, dgamma(x=x.tmp, shape=priors_normalgamma[[model]]$lambda0$shape, rate=priors_normalgamma[[model]]$lambda0$rate), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
mtext('Probability density', side=2, line=1.2, cex=1)
mtext(expression(mu[0]), side=1, line=2.7, cex=1)
mtext('NS3', side=2, line=3, cex=1)
pp <- 2; par(mai=c(.5,.35,.01,.25))
box.width <- diff(lims[2,])/nbins; box.edges <- seq(from=lims[2,1], to=lims[2,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[2,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[2,1], to=lims[2,2], length.out=1000); lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[model]]$lambda1$mean, sd=priors_normalgamma[[model]]$lambda1$sd), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
mtext(expression(mu[1]), side=1, line=2.7, cex=1)
pp <- 3; par(mai=c(.5,.3,.01,.3))
box.width <- diff(lims[3,])/nbins; box.edges <- seq(from=lims[3,1], to=lims[3,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[3,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[3,1], to=lims[3,2], length.out=1000); lines(x.tmp, dgamma(x=x.tmp, shape=priors_normalgamma[[model]]$sigma0$shape, rate=priors_normalgamma[[model]]$sigma0$rate), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
mtext(expression(sigma[0]), side=1, line=2.7, cex=1)
pp <- 4; par(mai=c(.5,.3,.01,.3))
box.width <- diff(lims[4,])/nbins; box.edges <- seq(from=lims[4,1], to=lims[4,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[4,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[4,1], to=lims[4,2], length.out=1000); lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[model]]$sigma1$mean, sd=priors_normalgamma[[model]]$sigma1$sd), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
mtext(expression(sigma[1]), side=1, line=2.7, cex=1)
pp <- 5; par(mai=c(.5,.3,.01,.3))
box.width <- diff(lims[5,])/nbins; box.edges <- seq(from=lims[5,1], to=lims[5,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[5,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[5,1], to=lims[5,2], length.out=1000); lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[model]]$xi0$mean, sd=priors_normalgamma[[model]]$xi0$sd), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
mtext(expression(xi[0]), side=1, line=2.7, cex=1)
pp <- 6; par(mai=c(.5,.3,.01,.3))
box.width <- diff(lims[6,])/nbins; box.edges <- seq(from=lims[6,1], to=lims[6,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[6,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[6,1], to=lims[6,2], length.out=1000); lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[model]]$xi1$mean, sd=priors_normalgamma[[model]]$xi1$sd), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
mtext(expression(xi[1]), side=1, line=2.7, cex=1)
##=============================
dev.off()

#===============================================================================
#


#
#===============================================================================
# SOM FIGURE S2 -- show the calibrated distributions of parameters using the
#                  normal and gamma priors


#===============================================================================
#


#
#===============================================================================
# SOM FIGURE S3 -- show the calibrated distributions of parameters using the
#                  uniform priors


#===============================================================================
#


#
#===============================================================================
# FIGURE 1 - Motivation figure, showing 30-year blocks and the estimated
#            distributions (box-whisker/box-lighter-box) for each block, for
#            each site.


#===============================================================================
#


#
#===============================================================================
# FIGURE 2 – Comparison of the empirical survival function calculated from the
#            observed tide gauge data at each site (red points, different columns)
#            against the modeled survival function in the ensemble median (black
#            points) and 5-95% credible range (error bars) for models (a) ST:
#            all parameters stationary, (b) NS1: lambda non-stationary, (c) NS3:
#            lambda and sigma non-stationary, (d) NS3: all non-stationary and
#            (e) the BMA-weighted ensemble.


#===============================================================================
#


#
#===============================================================================
# FIGURE 3 - Top row: current surge levels; bottom row: projected 2065 surge
#            levels relative to present. Columns: different sites.
#            Projected distributions of 100-year surge level by BMA, relative to
#            each of the individual model structures.


#===============================================================================
#


#
#===============================================================================
# FIGURE 4 - Box-whisker (or box-lighter-box) distributions (horizontally) of
#            100-year return level for varying lengths of data employed for the
#            BMA ensemble. Different sites are 3 horizontally oriented panels.


#===============================================================================
#


#
#===============================================================================
# FIGURE 5 - Bayesian model averaging weights (equation (XX)) for the four
#            candidate models, using (a) 30 years of tide gauge data from
#            Delfzijl, (b) 50 years of data, (c) 70 years of data, (d) 90 years
#            of data, (e) 110 years of data and (f) 137 years of data. Higher
#            values imply better model-data match. Different sites are different
#            columns.


#===============================================================================
#




#===============================================================================
# End
#===============================================================================

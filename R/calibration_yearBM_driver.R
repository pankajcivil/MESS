#===============================================================================
# Calibration of GEV, Naveau-i (and other...?) model(s) for storm surge
# at Delfzijl, The Netherlands.
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================

rm(list=ls())

niter_mcmc_prelim000 <- 1e3      # number of MCMC iterations (PRELIMINARY chains)
nnode_mcmc_prelim000 <- 1        # number of CPUs to use (PRELIMINARY chains)
niter_mcmc_prod000 <- 2e3        # number of MCMC iterations (PRODUCTION chains)
nnode_mcmc_prod000 <- 2          # number of CPUs to use (PRODUCTION chains)
gamma_mcmc000 <- 0.5             # speed of adaptation (0.5=faster, 1=slowest)

filename.priors   <- 'surge_priors_gev_nav_19Jun2017.rds'  # file holding the 'priors' object
filename.temperature <- 'temperature_forcing_19Jun2017.csv'  # temperature forcing used
filename.initvals <- 'surge_initialvalues_gev_nav_19Jun2017.rds'  # file holding the 'deoptim.delfzijl' object
filename.mles <- 'surge_MLEs_gev_nav_19Jun2017.rds' # file holding the 'mle.fits' object
filename.datacalib <- 'datacalib_05Jul2017.rds' # file holding the 'data_calib' object (calibration data)

output.dir <- '../output/'

#setwd('/storage/home/axw322/work/codes/EVT/R')
setwd('/home/scrim/axw322/codes/EVT/R')
#setwd('/Users/axw322/codes/EVT/R')
appen <- 'gev_nav'

# Name the calibrated parameters output file
today <- Sys.Date(); today=format(today,format="%d%b%Y")
appen <- 'gev-nav'
filename.parameters <- paste(output.dir,'evt_models_calibratedParameters_',appen,'_',today,'.nc',sep='')

#
#===============================================================================
# relevant libraries - do 'install.pacakges([library name])' if you do not have
# one yet
#===============================================================================
#

library(extRemes)
library(zoo)
library(adaptMCMC)
library(lhs)
library(DEoptim)
#library(ncdf4)

#
#===============================================================================
# read and process data for temperature (auxiliary covariate for nonstationary
# parameters)
# yields: temperature_forc, time_forc
#===============================================================================
#

print('reading temperature data...')

source('read_data_temperature.R')

print('...done.')

#
#===============================================================================
# read and process data for Dutch station (Delfzijl)
#===============================================================================
#

print('reading Delfzijl, Netherlands, tide gauge data...')

#source('read_data_tidegauge_delfzijl.R')
data_calib <- readRDS(paste(output.dir,filename.datacalib,sep=''))

print('...done.')

#
#===============================================================================
# set up GEV and Naveau (i) model parameters
#===============================================================================
#

print('setting up GEV and Naveau-i model parameters for DE optimization...')

priors <- readRDS(paste(output.dir,filename.priors,sep=''))
deoptim.delfzijl <- readRDS(paste(output.dir,filename.initvals,sep=''))
mle.fits <- readRDS(paste(output.dir,filename.mles,sep=''))

source('parameter_setup_yearBM.R')

print('...done.')

source('likelihood_gev.R')
source('likelihood_naveau.R')

#
#===============================================================================
# set up and run PRELIMINARY Markov chain Monte Carlo (MCMC) calibration
#===============================================================================
#

# first, do a set of single-chain preliminary calibrations to get estimates of
# the jump covariance matrix

nnode_mcmc <- nnode_mcmc_prelim000
niter_mcmc <- niter_mcmc_prelim000
gamma_mcmc <- gamma_mcmc000
startadapt_mcmc <- max(500,round(0.05*niter_mcmc))
stopadapt_mcmc <- round(niter_mcmc*1.0)
accept_mcmc_few <- 0.44         # optimal for only one parameter
accept_mcmc_many <- 0.234       # optimal for many parameters
amcmc_prelim <- vector('list', nmodel); names(amcmc_prelim) <- types.of.model

for (model in types.of.gev) {

  print(paste('Starting preliminary calibration for model ',model,' (',nnode_mcmc,' cores x ',niter_mcmc,' iterations)...', sep=''))

  initial_parameters <- as.numeric(deoptim.delfzijl[[model]])
  if(model=='gev3') {auxiliary <- NULL
  } else {auxiliary <- trimmed_forcing(data_calib$gev_year$year, time_forc, temperature_forc)$temperature}
  accept_mcmc <- accept_mcmc_many + (accept_mcmc_few - accept_mcmc_many)/length(parnames_all[[model]])
  step_mcmc <- as.numeric(0.05*apply(X=mle.fits[[model]], MARGIN=2, FUN=sd))
  tbeg=proc.time()
  amcmc_prelim[[model]] = MCMC(log_post_gev, niter_mcmc, initial_parameters,
                            adapt=TRUE, acc.rate=accept_mcmc, scale=step_mcmc,
                            gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                            parnames=parnames_all[[model]], data_calib=data_calib$gev_year$lsl_max,
                            priors=priors, auxiliary=auxiliary, model=model)
  tend=proc.time()

  print(paste('... done. Took ',round(as.numeric(tend-tbeg)[3]/60,2),' minutes', sep=''))
}

for (model in types.of.nav) {

  print(paste('Starting preliminary calibration for model ',model,' (',nnode_mcmc,' cores x ',niter_mcmc,' iterations)...', sep=''))

  initial_parameters <- as.numeric(deoptim.delfzijl[[model]])
  if('lkappa'  %in% parnames_all[[model]]) {initial_parameters[match('lkappa' ,parnames_all[[model]])] <- log(initial_parameters[match('lkappa' ,parnames_all[[model]])])}
  if('lkappa0' %in% parnames_all[[model]]) {initial_parameters[match('lkappa0',parnames_all[[model]])] <- log(initial_parameters[match('lkappa0',parnames_all[[model]])])}
  if(model=='nav3') {auxiliary <- NULL
  } else {auxiliary <- trimmed_forcing(data_calib$gev_year$year, time_forc, temperature_forc)$temperature}
  accept_mcmc <- accept_mcmc_many + (accept_mcmc_few - accept_mcmc_many)/length(parnames_all[[model]])
  step_mcmc <- as.numeric(0.05*apply(X=mle.fits[[model]], MARGIN=2, FUN=sd))
  tbeg=proc.time()
  amcmc_prelim[[model]] = MCMC(log_post_naveau, niter_mcmc, initial_parameters,
                            adapt=TRUE, acc.rate=accept_mcmc, scale=step_mcmc,
                            gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                            parnames=parnames_all[[model]], data_calib=data_calib$gev_year$lsl_max,
                            priors=priors, auxiliary=auxiliary, model=model, Tmax=Tmax)
  tend=proc.time()

  print(paste('... done. Took ',round(as.numeric(tend-tbeg)[3]/60,2),' minutes', sep=''))
}

#
#===============================================================================
# set up and run PRODUCTION MCMC calibration
#===============================================================================
#

# then use these initial estimates of step_mcmc to launch the parallel chains
# (from amcmc_prelim[[model]]$cov.jump)

nnode_mcmc <- nnode_mcmc_prod000
niter_mcmc <- niter_mcmc_prod000
gamma_mcmc <- gamma_mcmc000
startadapt_mcmc <- max(500,round(0.05*niter_mcmc))
stopadapt_mcmc <- round(niter_mcmc*1.0)
accept_mcmc_few <- 0.44         # optimal for only one parameter
accept_mcmc_many <- 0.234       # optimal for many parameters
amcmc_out <- vector('list', nmodel); names(amcmc_out) <- types.of.model

for (model in types.of.model) {

  print(paste('Starting production calibration for model ',model,' (',nnode_mcmc,' cores x ',niter_mcmc,' iterations)...', sep=''))

  if (model %in% types.of.gev) {
    initial_parameters <- amcmc_prelim[[model]]$samples[amcmc_prelim[[model]]$n.sample,]
    if(model=='gev3') {auxiliary <- NULL
    } else {auxiliary <- trimmed_forcing(data_calib$gev_year$year, time_forc, temperature_forc)$temperature}
    accept_mcmc <- accept_mcmc_many + (accept_mcmc_few - accept_mcmc_many)/length(parnames_all[[model]])
    step_mcmc <- amcmc_prelim[[model]]$cov.jump
    if(nnode_mcmc==1) {
      # do single chain
      tbeg=proc.time()
      amcmc_out[[model]] = MCMC(log_post_gev, niter_mcmc, initial_parameters,
                            adapt=TRUE, acc.rate=accept_mcmc, scale=step_mcmc,
                            gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                            parnames=parnames_all[[model]], data_calib=data_calib$gev_year$lsl_max,
                            priors=priors, auxiliary=auxiliary, model=model)
      tend=proc.time()
    } else if(nnode_mcmc > 1) {
      # do parallel chains
      tbeg <- proc.time()
      amcmc_out[[model]] <- MCMC.parallel(log_post_gev, niter_mcmc, initial_parameters,
                            n.chain=4, n.cpu=4, packages='extRemes',
				            scale=step_mcmc, adapt=TRUE, acc.rate=accept_mcmc,
                            gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                            parnames=parnames_all[[model]], data_calib=data_calib$gev_year$lsl_max,
                            priors=priors, auxiliary=auxiliary, model=model)
      tend <- proc.time()
    }
  } else if (model %in% types.of.nav) {
    initial_parameters <- amcmc_prelim[[model]]$samples[amcmc_prelim[[model]]$n.sample,]
    # don't need to transform, since the preliminary chains sampled log(kappa)
#    if('lkappa'  %in% parnames_all[[model]]) {initial_parameters[match('lkappa' ,parnames_all[[model]])] <- log(initial_parameters[match('lkappa' ,parnames_all[[model]])])}
#    if('lkappa0' %in% parnames_all[[model]]) {initial_parameters[match('lkappa0',parnames_all[[model]])] <- log(initial_parameters[match('lkappa0',parnames_all[[model]])])}
    if(model=='nav3') {auxiliary <- NULL
    } else {auxiliary <- trimmed_forcing(data_calib$gev_year$year, time_forc, temperature_forc)$temperature}
    accept_mcmc <- accept_mcmc_many + (accept_mcmc_few - accept_mcmc_many)/length(parnames_all[[model]])
    step_mcmc <- amcmc_prelim[[model]]$cov.jump
    if(nnode_mcmc==1) {
      # do single chain
      tbeg=proc.time()
      amcmc_out[[model]] = MCMC(log_post_naveau, niter_mcmc, initial_parameters,
                            adapt=TRUE, acc.rate=accept_mcmc, scale=step_mcmc,
                            gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                            parnames=parnames_all[[model]], data_calib=data_calib$gev_year$lsl_max,
                            priors=priors, auxiliary=auxiliary, model=model, Tmax=Tmax)
      tend=proc.time()
    } else if(nnode_mcmc > 1) {
      # do parallel chains
      tbeg <- proc.time()
      amcmc_out[[model]] <- MCMC.parallel(log_post_naveau, niter_mcmc, initial_parameters, n.chain=4, n.cpu=4,
				             scale=step_mcmc, adapt=TRUE, acc.rate=accept_mcmc,
                             gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                             parnames=parnames_all[[model]], data_calib=data_calib$gev_year$lsl_max,
                             priors=priors, auxiliary=auxiliary, model=model, Tmax=Tmax)
      tend <- proc.time()
    }
  } else {print('error - unknown model type')}
  print(paste('... done. Took ',round(as.numeric(tend-tbeg)[3]/60,2),' minutes', sep=''))
}


# for now, just save RData file so you can play around with the GR stats,
# subtracting off burn-in, thinning (if needed), etc.

#
#===============================================================================
# save raw results
#===============================================================================
#

today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.everythingmcmc <- paste(output.dir,'kitchen_sink_mcmc_',appen,'_',today,'.RData', sep='')

print(paste('saving results as .RData file (',filename.everythingmcmc,') to read and use later...',sep=''))

save.image(file=filename.everythingmcmc)

print('...done.')

# test plot
if (FALSE) {
model <- 'nav4'
par(mfrow=c(6,2))
for (p in 1:length(parnames_all[[model]])) {
    plot(amcmc_out[[model]]$samples[,p], type='l', ylab=parnames_all[[model]][p])
    hist(amcmc_out[[model]]$samples[round(0.5*niter_mcmc):niter_mcmc,p], xlab=parnames_all[[model]][p], main='')
}
}

#
#===============================================================================
# convergence diagnostics
#===============================================================================
#

# Gelman and Rubin diagnostics - determine and chop off for burn-in
niter.test <- seq(from=round(0.1*niter_mcmc), to=niter_mcmc, by=round(0.05*niter_mcmc))
gr.test <- mat.or.vec(length(niter.test), nmodel)
gr.tmp <- rep(NA, length(niter.test))
colnames(gr.test) <- types.of.model


for (model in types.of.model) {
  if(nnode_mcmc == 1) {
    # don't do GR stats, just cut off first half of chains
    print('only one chain; will lop off first half for burn-in instead of doing GR diagnostics')
  } else if(nnode_mcmc > 1) {
    # this case is FAR more fun
    # accumulate the names of the soon-to-be mcmc objects
    string.mcmc.list <- 'mcmc1'
    for (m in 2:nnode_mcmc) {
      string.mcmc.list <- paste(string.mcmc.list, ', mcmc', m, sep='')
    }
    for (i in 1:length(niter.test)) {
      for (m in 1:nnode_mcmc) {
        # convert each of the chains into mcmc object
        eval(parse(text=paste('mcmc',m,' <- as.mcmc(amcmc_out[[model]][[m]]$samples[1:niter.test[i],])', sep='')))
      }
      eval(parse(text=paste('mcmc_chain_list = mcmc.list(list(', string.mcmc.list , '))', sep='')))

      gr.test[i,model] <- as.numeric(gelman.diag(mcmc_chain_list)[2])
    }
  } else {print('error - nnode_mcmc < 1 makes no sense')}
}

# Monitor posterior 5, 50 and 95% quantiles for drift
# Only checking for one of the chains
quant <- vector('list', nmodel); names(quant) <- types.of.model
names.monitor <- c('q05', 'q50', 'q95')
for (model in types.of.model) {
  quant[[model]] <- vector('list', 3); names(quant[[model]]) <- names.monitor
  for (q in names.monitor) {
    quant[[model]][[q]] <- mat.or.vec(length(niter.test)-1, length(parnames_all[[model]]))
    for (i in 1:(length(niter.test)-1)) {
      if(nnode_mcmc==1) {
        quant[[model]][[q]][i,] <- apply(X=amcmc_out[[model]]$samples[niter.test[i]:niter_mcmc,], MARGIN=2, FUN=quantile, probs=as.numeric(substr(q, 2,3))*0.01)
      } else {
        quant[[model]][[q]][i,] <- apply(X=amcmc_out[[model]][[1]]$samples[niter.test[i]:niter_mcmc,], MARGIN=2, FUN=quantile, probs=as.numeric(substr(q, 2,3))*0.01)
      }
    }
  }
}

if(FALSE) {
# examples monitoring of stability of quantiles:
model <- 'gev3'
par(mfrow=c(3,2))
for (p in 1:length(parnames_all[[model]])) {
  ran <- max(quant[[model]]$q95[,p])-min(quant[[model]]$q05[,p])
  lb <- min(quant[[model]]$q05[,p]) - 0.05*ran; ub <- max(quant[[model]]$q95[,p]) + 0.05*ran
  plot(niter.test[1:(length(niter.test)-1)], quant[[model]]$q50[,p], type='l',
    ylim=c(lb,ub), ylab=parnames_all[[model]][p], xlab='From HERE to end of chain')
  lines(niter.test[1:(length(niter.test)-1)], quant[[model]]$q05[,p], lty=2); lines(niter.test[1:(length(niter.test)-1)], quant[[model]]$q95[,p], lty=2);
}
model <- 'nav5'
par(mfrow=c(3,2))
for (p in 1:length(parnames_all[[model]])) {
  ran <- max(quant[[model]]$q95[,p])-min(quant[[model]]$q05[,p])
  lb <- min(quant[[model]]$q05[,p]) - 0.05*ran; ub <- max(quant[[model]]$q95[,p]) + 0.05*ran
  plot(niter.test[1:(length(niter.test)-1)], quant[[model]]$q50[,p], type='l',
    ylim=c(lb,ub), ylab=parnames_all[[model]][p], xlab='From HERE to end of chain')
  lines(niter.test[1:(length(niter.test)-1)], quant[[model]]$q05[,p], lty=2); lines(niter.test[1:(length(niter.test)-1)], quant[[model]]$q95[,p], lty=2);
}
# the thing to note is that these stabilize as you include more members (i.e.,
# as you move from right to left)
}


# Heidelberger and Welch diagnostics?
hw.diag <- vector('list', nmodel); names(hw.diag) <- types.of.model
for (model in types.of.model) {
  hw.diag[[model]] <- heidel.diag(as.mcmc(amcmc_out[[model]][[1]]$samples), eps=0.1, pvalue=0.05)
}

#
#===============================================================================
# Chop off burn-in
#===============================================================================
#

# Note: here, we are only using the Gelman and Rubin diagnostic. But this is
# only after looking at the quantile stability as iterations increase, as well
# as the Heidelberger and Welch diagnostics, which suggest the chains are okay.
# 'ifirst' is the first spot where the GR stat gets to and stays below gr.max
# for all of the models.
if(nnode_mcmc==1) {
  ifirst <- round(0.5*niter_mcmc)
} else {
  gr.max <- 1.1
  lgr <- rep(NA, length(niter.test))
  for (i in 1:length(niter.test)) {lgr[i] <- all(gr.test[i,] < gr.max)}
  ifirst <- NULL
  for (i in seq(from=length(niter.test), to=1, by=-1)) {
    if( all(lgr[i:length(lgr)]) ) {ifirst <- niter.test[i]}
  }
}

chains_burned <- vector('list', nmodel); names(chains_burned) <- types.of.model
for (model in types.of.model) {
  if(nnode_mcmc > 1) {
    chains_burned[[model]] <- vector('list', nnode_mcmc)
    for (m in 1:nnode_mcmc) {
      chains_burned[[model]][[m]] <- amcmc_out[[model]][[m]]$samples[(ifirst+1):niter_mcmc,]
    }
  } else {
    chains_burned[[model]] <- amcmc_out[[model]]$samples[(ifirst+1):niter_mcmc,]
  }
}


#
#===============================================================================
# possible thinning?
#===============================================================================
#

# If no thinning, then this initialization will remain
chains_burned_thinned <- chains_burned

if(FALSE) {#==========================

acf_cutoff <- 0.05
lag_max <- 0.01*niter_mcmc # if we end up with fewer than 100 samples, what are we even doing?
niter_thin <- rep(0, nmodel); names(niter_thin) <- types.of.model
for (model in types.of.model) {
  for (p in 1:length(parnames_all[[model]])) {
    if(nnode_mcmc > 1) {acf_tmp <- acf(x=chains_burned[[model]][[1]][,p], plot=FALSE, lag.max=lag_max)}
    else {acf_tmp <- acf(x=chains_burned[[model]][,p], plot=FALSE, lag.max=lag_max)}
    niter_thin[[model]] <- max(niter_thin[[model]], acf_tmp$lag[which(acf_tmp$acf < acf_cutoff)[1]])
  }
  nthin <- max(niter_thin, na.rm=TRUE)
  if(nnode_mcmc > 1) {
    for (m in 1:nnode_mcmc) {
      chains_burned_thinned[[model]][[m]] <- chains_burned[[model]][[m]][seq(from=1, to=nrow(chains_burned[[model]][[m]]), by=nthin),]
    }
  } else {
    chains_burned_thinned[[model]] <- chains_burned[[model]][seq(from=1, to=nrow(chains_burned[[model]]), by=nthin),]
  }
}

}#====================================


# thin to a target number of samples?
if(TRUE) {#===========================

n.sample <- 10000   # dike ring 15 safety standard is 1/2000, so at least 2000...

for (model in types.of.model) {
  if(nnode_mcmc == 1) {
    ind.sample <- sample(x=1:nrow(chains_burned[[model]]), size=n.sample, replace=FALSE)
    chains_burned_thinned[[model]] <- chains_burned[[model]][ind.sample,]
  } else {
    n.sample.sub <- rep(NA, nnode_mcmc)
    # the case where desired sample size is divisible by the number of chains
    if(round(n.sample/nnode_mcmc) == n.sample/nnode_mcmc) {
      n.sample.sub[1:nnode_mcmc] <- n.sample/nnode_mcmc
    } else {
    # the case where it is not
      n.sample.sub[2:nnode_mcmc] <- round(n.sample/nnode_mcmc)
      n.sample.sub[1] <- n.sample - sum(n.sample.sub[2:nnode_mcmc])
    }
    for (m in 1:nnode_mcmc) {
      ind.sample <- sample(x=1:nrow(chains_burned[[model]][[m]]), size=n.sample.sub[m], replace=FALSE)
      chains_burned_thinned[[model]][[m]] <- chains_burned[[model]][[m]][ind.sample,]
    }
  }
}

}#====================================


# Combine all of the chains from 'ifirst' to 'niter_mcmc' into a potpourri of
# [alleged] samples from the posterior. Only saving the transition covariance
# matrix for one of the chains (if in parallel).
parameters.posterior <- vector('list', nmodel); names(parameters.posterior) <- types.of.model
covjump.posterior <- vector('list', nmodel); names(covjump.posterior) <- types.of.model

for (model in types.of.model) {
  if(nnode_mcmc==1) {
    parameters.posterior[[model]] <- chains_burned_thinned[[model]]
    covjump.posterior[[model]] <- amcmc_out[[model]]$cov.jump
  } else {
    parameters.posterior[[model]] <- chains_burned_thinned[[model]][[1]]
    covjump.posterior[[model]] <- amcmc_out[[model]][[1]]$cov.jump
    for (m in 2:nnode_mcmc) {
      parameters.posterior[[model]] <- rbind(parameters.posterior[[model]], chains_burned_thinned[[model]][[m]])
    }
  }
}

#
#===============================================================================
# write output file
#===============================================================================
#

## Get maximum length of parameter name, for width of array to write to netcdf
## this code will write an n.parameters (rows) x n.ensemble (columns) netcdf file
## to get back into the shape we expect, just transpose it
lmax=0
for (model in types.of.model) {for (i in 1:length(parnames_all[[model]])){lmax=max(lmax,nchar(parnames_all[[model]][i]))}}

dim.parameters <- vector('list', nmodel); names(dim.parameters) <- types.of.model
dim.parnames   <- vector('list', nmodel); names(dim.parnames)   <- types.of.model
var.parameters <- vector('list', nmodel); names(var.parameters) <- types.of.model
var.parnames   <- vector('list', nmodel); names(var.parnames)   <- types.of.model
var.covjump    <- vector('list', nmodel); names(var.covjump)    <- types.of.model
dim.ensemble   <- vector('list', nmodel); names(dim.ensemble)   <- types.of.model
dim.name <- ncdim_def('name.len', '', 1:lmax, unlim=FALSE)
dim.time <- ncdim_def('ntime', '', (time_forc), unlim=FALSE)
var.time <- ncvar_def('time', '', dim.time, -999)
var.temperature <- ncvar_def('temperature', 'Global mean temperature relative to 1901-2000', dim.time, -999)
for (model in types.of.model) {
  dim.parameters[[model]] <- ncdim_def(paste('n.parameters.',model,sep=''), '', 1:length(parnames_all[[model]]), unlim=FALSE)
  dim.ensemble[[model]]   <- ncdim_def(paste('n.ensemble.',model,sep=''), 'ensemble member', 1:nrow(parameters.posterior[[model]]), unlim=FALSE)
  var.parameters[[model]] <- ncvar_def(paste('parameters.',model,sep=''), '', list(dim.parameters[[model]],dim.ensemble[[model]]), -999)
  var.parnames[[model]]   <- ncvar_def(paste('parnames.',model,sep=''), '', list(dim.name,dim.parameters[[model]]), prec='char')
  var.covjump[[model]]    <- ncvar_def(paste('covjump.',model,sep=''), '', list(dim.parameters[[model]],dim.parameters[[model]]), -999)
}
output.to.file <- vector('list', length(var.parameters) + length(var.parnames) + length(var.covjump) + 2)
output.to.file[[1]] <- var.time
output.to.file[[2]] <- var.temperature
cnt <- 3
for (model in types.of.model) {
  output.to.file[[cnt]] <- var.parameters[[model]]
  output.to.file[[cnt+1]] <- var.parnames[[model]]
  output.to.file[[cnt+2]] <- var.covjump[[model]]
  cnt <- cnt + 3
}
#outnc <- nc_create(filename.parameters, list(var.parameters, var.parnames, var.covjump))
outnc <- nc_create(filename.parameters, output.to.file)
ncvar_put(outnc, var.time, time_forc)
ncvar_put(outnc, var.temperature, temperature_forc)
for (model in types.of.model) {
  ncvar_put(outnc, var.parameters[[model]], t(parameters.posterior[[model]]))
  ncvar_put(outnc, var.parnames[[model]], parnames_all[[model]])
  ncvar_put(outnc, var.covjump[[model]], covjump.posterior[[model]])
}
nc_close(outnc)

# link straight into flood risk
parameters <- parameters.posterior

# ... or...

#
#===============================================================================
# read a previous calibration output file, to use for flood risk
#===============================================================================
#

types.of.gev <- c('gev3','gev4','gev5','gev6')
types.of.nav <- c('nav3','nav4','nav5','nav6')
types.of.model <- c(types.of.gev, types.of.nav)

parameters <- vector('list', nmodel); names(parameters) <- types.of.model
parnames_all <- vector('list', nmodel); names(parnames_all) <- types.of.model
covjump <- vector('list', nmodel); names(covjump) <- types.of.model
ncdata <- nc_open('evt_models_calibratedParameters__13Jun2017.nc')
  time_forc <- ncvar_get(ncdata, 'time')
  temperature_forc <- ncvar_get(ncdata, 'temperature')
  for (model in types.of.model) {
    parameters[[model]]   <- t(ncvar_get(ncdata, paste('parameters.',model,sep='')))
    parnames_all[[model]] <- ncvar_get(ncdata, paste('parnames.',model,sep=''))
    covjump[[model]]      <- ncvar_get(ncdata, paste('covjump.',model,sep=''))
  }
nc_close(ncdata)


}


#
#===============================================================================
# End
#===============================================================================
#

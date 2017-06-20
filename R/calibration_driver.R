#===============================================================================
# Calibration of GEV, Naveau-i (and other...?) model(s) for storm surge
# at Delfzijl, The Netherlands.
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================

rm(list=ls())

niter_mcmc_prelim000 <- 2e3      # number of MCMC iterations (PRELIMINARY chains)
nnode_mcmc_prelim000 <- 1        # number of CPUs to use (PRELIMINARY chains)
niter_mcmc_prod000 <- 2e3        # number of MCMC iterations (PRODUCTION chains)
nnode_mcmc_prod000 <- 1          # number of CPUs to use (PRODUCTION chains)
gamma_mcmc000 <- 0.5             # speed of adaptation (0.5=faster, 1=slowest)
filename.priors   <- 'surge_priors_gev_nav_19Jun2017.rds'  # file holding the 'priors' object
filename.temperature <- 'temperature_forcing_19Jun2017.csv'  # temperature forcing used
filename.initvals <- 'surge_initialvalues_gev_nav_19Jun2017.rds'  # file holding the 'deoptim.delfzijl' object
#setwd('/storage/home/axw322/work/codes/EVT/R')
setwd('/Users/axw322/codes/EVT/R')
appen <- 'gev_nav'

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

source('read_data_tidegauge_delfzijl.R')

print('...done.')

#
#===============================================================================
# set up GEV and Naveau (i) model parameters
#===============================================================================
#

print('setting up GEV and Naveau-i model parameters for DE optimization...')

priors <- readRDS(filename.priors)
deoptim.delfzijl <- readRDS(filename.initvals)

source('parameter_setup.R')

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
amcmc_prelim <- vector('list', length(types.of.model)); names(amcmc_prelim) <- types.of.model

for (model in types.of.gev) {

  print(paste('Starting preliminary calibration for model ',model,' (',nnode_mcmc,' cores x ',niter_mcmc,' iterations)...', sep=''))

  initial_parameters <- as.numeric(deoptim.delfzijl[[model]])
  if(model=='gev3') {auxiliary <- NULL
  } else {auxiliary <- trimmed_forcing(data_calib$year_unique, time_forc, temperature_forc)$temperature}
  accept_mcmc <- accept_mcmc_many + (accept_mcmc_few - accept_mcmc_many)/length(parnames_all[[model]])
  step_mcmc <- as.numeric(0.05*apply(X=mle.fits[[model]], MARGIN=2, FUN=sd))
  tbeg=proc.time()
  amcmc_prelim[[model]] = MCMC(log_post_gev, niter_mcmc, initial_parameters,
                            adapt=TRUE, acc.rate=accept_mcmc, scale=step_mcmc,
                            gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                            parnames=parnames_all[[model]], data_calib=data_calib$lsl_max,
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
  } else {auxiliary <- trimmed_forcing(data_calib$year_unique, time_forc, temperature_forc)$temperature}
  accept_mcmc <- accept_mcmc_many + (accept_mcmc_few - accept_mcmc_many)/length(parnames_all[[model]])
  step_mcmc <- as.numeric(0.05*apply(X=mle.fits[[model]], MARGIN=2, FUN=sd))
  tbeg=proc.time()
  amcmc_prelim[[model]] = MCMC(log_post_naveau, niter_mcmc, initial_parameters,
                            adapt=TRUE, acc.rate=accept_mcmc, scale=step_mcmc,
                            gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                            parnames=parnames_all[[model]], data_calib=data_calib$lsl_max,
                            priors=priors, auxiliary=auxiliary, model=model)
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
amcmc_out <- vector('list', length(types.of.model)); names(amcmc_out) <- types.of.model

for (model in types.of.model) {

  print(paste('Starting production calibration for model ',model,' (',nnode_mcmc,' cores x ',niter_mcmc,' iterations)...', sep=''))

  if (model %in% types.of.gev) {
    initial_parameters <- amcmc_prelim[[model]]$samples[amcmc_prelim[[model]]$n.sample,]
    if(model=='gev3') {auxiliary <- NULL
    } else {auxiliary <- trimmed_forcing(data_calib$year_unique, time_forc, temperature_forc)$temperature}
    accept_mcmc <- accept_mcmc_many + (accept_mcmc_few - accept_mcmc_many)/length(parnames_all[[model]])
    step_mcmc <- amcmc_prelim[[model]]$cov.jump
    if(nnode_mcmc==1) {
      # do single chain
      tbeg=proc.time()
      amcmc_out[[model]] = MCMC(log_post_gev, niter_mcmc, initial_parameters,
                            adapt=TRUE, acc.rate=accept_mcmc, scale=step_mcmc,
                            gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                            parnames=parnames_all[[model]], data_calib=data_calib$lsl_max,
                            priors=priors, auxiliary=auxiliary, model=model)
      tend=proc.time()
    } else if(nnode_mcmc > 1) {
      # do parallel chains
      tbeg <- proc.time()
      amcmc_out[[model]] <- MCMC.parallel(log_post_gev, niter_mcmc, initial_parameters,
                            n.chain=4, n.cpu=4, packages='extRemes',
				            scale=step_mcmc, adapt=TRUE, acc.rate=accept_mcmc,
                            gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                            parnames=parnames_all[[model]], data_calib=data_calib$lsl_max,
                            priors=priors, auxiliary=auxiliary, model=model)
      tend <- proc.time()
    }
  } else if (model %in% types.of.nav) {
    initial_parameters <- amcmc_prelim[[model]]$samples[amcmc_prelim[[model]]$n.sample,]
    # don't need to transform, since the preliminary chains sampled log(kappa)
#    if('lkappa'  %in% parnames_all[[model]]) {initial_parameters[match('lkappa' ,parnames_all[[model]])] <- log(initial_parameters[match('lkappa' ,parnames_all[[model]])])}
#    if('lkappa0' %in% parnames_all[[model]]) {initial_parameters[match('lkappa0',parnames_all[[model]])] <- log(initial_parameters[match('lkappa0',parnames_all[[model]])])}
    if(model=='nav3') {auxiliary <- NULL
    } else {auxiliary <- trimmed_forcing(data_calib$year_unique, time_forc, temperature_forc)$temperature}
    accept_mcmc <- accept_mcmc_many + (accept_mcmc_few - accept_mcmc_many)/length(parnames_all[[model]])
    step_mcmc <- amcmc_prelim[[model]]$cov.jump
    if(nnode_mcmc==1) {
      # do single chain
      tbeg=proc.time()
      amcmc_out[[model]] = MCMC(log_post_naveau, niter_mcmc, initial_parameters,
                            adapt=TRUE, acc.rate=accept_mcmc, scale=step_mcmc,
                            gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                            parnames=parnames_all[[model]], data_calib=data_calib$lsl_max,
                            priors=priors, auxiliary=auxiliary, model=model)
      tend=proc.time()
    } else if(nnode_mcmc > 1) {
      # do parallel chains
      tbeg <- proc.time()
      amcmc_out[[model]] <- MCMC.parallel(log_post_naveau, niter_mcmc, initial_parameters, n.chain=4, n.cpu=4,
				             scale=step_mcmc, adapt=TRUE, acc.rate=accept_mcmc,
                             gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                             parnames=parnames_all[[model]], data_calib=data_calib$lsl_max,
                             priors=priors, auxiliary=auxiliary, model=model)
      tend <- proc.time()
    }
  } else {print('error - unknown model type')}
  print(paste('... done. Took ',round(as.numeric(tend-tbeg)[3]/60,2),' minutes', sep=''))
}


# for now, just save RData file so you can play around with the GR stats,
# subtracting off burn-in, thinning (if needed), etc.


# later, can actually turn this into a simple pipeline that will spit out a
# calibrated parameters netcdf file


#
#===============================================================================
# save results (testing?)
#===============================================================================
#

today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.everythingmcmc <- paste('kitchen_sink_mcmc_',appen,'_',today,'.RData', sep='')

print(paste('saving results as .RData file (',filename.everythingmcmc,') to read and use later...',sep=''))

save.image(file=filename.everything)

print('...done.')





#===============================================================================
# DANGER: UNDER CONSTRUCTION
#===============================================================================


if (FALSE) {



# test plot
if (FALSE) {
model <- 'nav4'
par(mfrow=c(6,2))
for (p in 1:length(parnames_all[[model]])) {
    plot(amcmc_out[[model]]$samples[,p], type='l', ylab=parnames_all[[model]][p])
    hist(amcmc_out[[model]]$samples[round(0.5*niter_mcmc):niter_mcmc,p], xlab=parnames_all[[model]][p], main='')
}
}


# Gelman and Rubin diagnostics - determine and chop off for burn-in
if(FALSE) {

# TODO

    # TODO - tidy up and figure out how to store, burn-in, etc.

    chain1=amcmc.par1[[1]]$samples
    chain2=amcmc.par1[[2]]$samples
    chain3=amcmc.par1[[3]]$samples
    chain4=amcmc.par1[[4]]$samples

  niter.test = seq(from=10000, to=nrow(chain1), by=5000)
  gr.test = rep(NA,length(niter.test))

  for (i in 1:length(niter.test)){
    mcmc1 = as.mcmc(chain1[1:niter.test[i],])
    mcmc2 = as.mcmc(chain2[1:niter.test[i],])
    mcmc3 = as.mcmc(chain3[1:niter.test[i],])
    mcmc4 = as.mcmc(chain4[1:niter.test[i],])
    mcmc_chain_list = mcmc.list(list(mcmc1, mcmc2, mcmc3, mcmc4))
    gr.test[i] = gelman.diag(mcmc_chain_list)[2]
  }
}

# Heidelberger and Welch diagnostics

# TODO

# heidel.diag(as.mcmc(chain1), eps=0.1, pvalue=0.05)



#
#===============================================================================
# Once satisfied they are converged, write samples to netCDF file
#===============================================================================
#

parameters.posterior <- vector('list', length(types.of.model)); names(parameters.posterior) <- types.of.model
covjump.posterior <- vector('list', length(types.of.model)); names(covjump.posterior) <- types.of.model
for (model in types.of.model) {

 # TODO
 # TODO

  # for now, just using from amcmc_prelim[[model]]$samples
  parameters.posterior[[model]] <- amcmc_prelim[[model]]$samples[round(0.5*nrow(amcmc_prelim[[model]]$samples)):nrow(amcmc_prelim[[model]]$samples),]
  covjump.posterior[[model]] <- amcmc_prelim[[model]]$cov.jump

 # TODO
 # TODO

}

## Get maximum length of parameter name, for width of array to write to netcdf
## this code will write an n.parameters (rows) x n.ensemble (columns) netcdf file
## to get back into the shape we expect, just transpose it
lmax=0
for (model in types.of.model) {for (i in 1:length(parnames_all[[model]])){lmax=max(lmax,nchar(parnames_all[[model]][i]))}}

## Name the file
today <- Sys.Date(); today=format(today,format="%d%b%Y")
appen <- ''
filename.parameters <- paste('evt_models_calibratedParameters_',appen,'_',today,'.nc',sep='')

dim.parameters <- vector('list', length(types.of.model)); names(dim.parameters) <- types.of.model
dim.parnames   <- vector('list', length(types.of.model)); names(dim.parnames)   <- types.of.model
var.parameters <- vector('list', length(types.of.model)); names(var.parameters) <- types.of.model
var.parnames   <- vector('list', length(types.of.model)); names(var.parnames)   <- types.of.model
var.covjump    <- vector('list', length(types.of.model)); names(var.covjump)    <- types.of.model
dim.ensemble   <- vector('list', length(types.of.model)); names(dim.ensemble)   <- types.of.model
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

parameters <- vector('list', length(types.of.model)); names(parameters) <- types.of.model
parnames_all <- vector('list', length(types.of.model)); names(parnames_all) <- types.of.model
covjump <- vector('list', length(types.of.model)); names(covjump) <- types.of.model
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

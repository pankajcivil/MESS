#===============================================================================
# want to do this for each site, for each data length
# estimate the BMA weights by the fraction of time the converged chain spends in
# each of the models
#===============================================================================

setwd('~/codes/EVT/R')

library(extRemes)
library(MASS)

# start off with a bunch of stuff defined, but will want to edit to start this
# one from scratch (well, take in the previuos individual model MCMC results)

#load('../output/analysis_inprogress.RData')

source('read_data_temperature.R')

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
all.data <- rep(NA, nsites); names(all.data) <- site.names
data.lengths <- vector('list', nsites); names(data.lengths) <- site.names
data.sites <- vector('list', nsites); names(data.sites) <- site.names


# initialize with the maximum likelihood guys and transition covariance matrices
# from each of the other chains

mcmc.balboa <- '../output/everything_mcmc_ppgpd-experiments_balboa_normalgamma_28Jul2017.RData'
mcmc.norfolk <- '../output/everything_mcmc_ppgpd-experiments_norfolk_normalgamma_27Jul2017.RData'
mcmc.delfzijl <- '../output/everything_mcmc_ppgpd-experiments_delfzijl_normalgamma_29Jul2017.RData'

results <- calibration.data <- vector('list', nsites)
names(results) <- names(calibration.data) <- site.names
trans.cov <- mp.parameters <- list.init

for (site in site.names) {
  if(site=='Balboa') {load(mcmc.balboa)
  } else if(site=='Norfolk') {load(mcmc.norfolk)
  } else if(site=='Delfzijl') {load(mcmc.delfzijl)}
  calibration.data[[site]] <- data_calib
  trans.cov[[site]][[dd]][[model]] <- amcmc_out[[dd]][[model]][[1]]$cov.jump
  mp.parameters[[site]][[dd]][[model]] <- amcmc_out[[dd]][[model]][[1]]$samples[which.max(amcmc_out[[dd]][[model]][[1]]$log.p),]
}

#===============================================================================

# initialize things for the BMA MCMC

niter_mcmc <- 1e4
mcmc.model <- rep(NA, niter_mcmc)
mcmc.logpost <- rep(NA, niter_mcmc)
mcmc.parameters <- vector('list', nmodel); names(mcmc.parameters) <- types.of.gpd
mcmc.parameters.last <- vector('list', nmodel); names(mcmc.parameters.last) <- types.of.gpd
for (model in types.of.gpd) {mcmc.parameters[[model]] <- mat.or.vec(niter_mcmc, length(parnames[[model]])); names(mcmc.parameters[[model]]) <- parnames[[model]]}

auxiliary <- vector('list', nmodel); names(auxiliary) <- types.of.gpd
for (model in types.of.gpd) {
  if(model=='gpd3') {auxiliary[[model]] <- NULL
  } else {auxiliary[[model]] <- trimmed_forcing(calibration.data[[site]][[data.lengths[[site]][dd]]]$year, time_forc, temperature_forc)$temperature}
}

# start a chain from each of the candidate models
# TODO

# initialize the parameter iteration and model iteration

site <- 'Delfzijl'
#dd <- all.data[[site]]
dd <- 1
model.init <- model.old <- 2
model.numbers <- c(1,2,3,4)
mcmc.model[1] <- model.init

# define the model's neighbors, where can we transition to from each model?
model.nbrs <- vector('list',nmodel); names(model.nbrs) <- types.of.gpd
model.nbrs$gpd3 <- c(1,2)
model.nbrs$gpd4 <- c(1,2,3)
model.nbrs$gpd5 <- c(2,3,4)
model.nbrs$gpd6 <- c(3,4)

for (model in types.of.gpd) {
  mcmc.parameters.last[[model]] <- mvrnorm(n=1, mu=mp.parameters[[site]][[dd]][[model]], Sigma=trans.cov[[site]][[dd]][[model]])
  logpost.init <- log_post_ppgpd(parameters=mcmc.parameters.last[[model]], parnames=parnames[[model]],
                               data_calib=calibration.data[[site]][[data.lengths[[site]][dd]]],
                               priors=priors, model=model, auxiliary=auxiliary[[model]])
  while(is.na(logpost.init)) {
    mcmc.parameters.last[[model]] <- mvrnorm(n=1, mu=mp.parameters[[site]][[dd]][[model]], Sigma=trans.cov[[site]][[dd]][[model]])
    logpost.init <- log_post_ppgpd(parameters=mcmc.parameters.last[[model]], parnames=parnames[[model]],
                                 data_calib=calibration.data[[site]][[data.lengths[[site]][dd]]],
                                 priors=priors, model=model, auxiliary=auxiliary[[model]])
  }
}

mcmc.parameters[[model.init]][1,] <- parameters.old <- mcmc.parameters.last[[model.init]]

# initialize log-posterior vector (conditioned on the model structure)

#logpost.old <- log_post_ppgpd(parameters=parameters.old, parnames=parnames[[mcmc.model[1]]],
#                             data_calib=calibration.data[[site]][[data.lengths[[site]][dd]]],
#                             priors=priors, model=types.of.gpd[[mcmc.model[1]]], auxiliary=auxiliary[[types.of.gpd[model]]])
logpost.old <- sample(x=lpost[[site]][[dd]][[model.init]], size=1)
mcmc.logpost[1] <- logpost.old

# MCMC iteration

pb <- txtProgressBar(min=0,max=niter_mcmc,initial=0,style=3)
for (iter in 2:niter_mcmc) {

  # propose a model from {1,2,3,4} (corresponding to {ST, NS1, NS2, NS3})

  model.new <- sample(x=model.nbrs[[model.old]], size=1)

  # propose a set of parameters corresponding to that model structure

  #parameters.new <- mvrnorm(n=1, mu=mcmc.parameters.last[[model.new]], Sigma=trans.cov[[site]][[dd]][[model.new]])

  # calculate log-posterior with that model

  #logpost.new <- log_post_ppgpd(parameters=parameters.new, parnames=parnames[[model.new]],
  #                             data_calib=calibration.data[[site]][[data.lengths[[site]][dd]]],
  #                             priors=priors, model=types.of.gpd[model.new], auxiliary=auxiliary[[types.of.gpd[model.new]]])

  # instead, propose a log-posterior from the model (already 'ran' them, so do
  # not need to do that again)

  logpost.new <- sample(x=lpost[[site]][[dd]][[model.new]], size=1)

  # calcualte acceptance probability

  p_accept <- min(0, logpost.new - logpost.old)

  # make the transition... maybe!

  log.unif.rnd <- log(runif(1))
  if(log.unif.rnd < p_accept) {
    # accept the proposed parameters, and accept the proposed model
    mcmc.logpost[iter] <- logpost.new
    mcmc.model[iter] <- model.new
    #mcmc.parameters[[model.new]][iter,] <- mcmc.parameters.last[[model.new]] <- parameters.new
  } else {
    # reject the proposal and stick with the previous parameters and model
    mcmc.logpost[iter] <- logpost.old
    mcmc.model[iter] <- model.old
    #mcmc.parameters[[model.old]][iter,] <- parameters.old
  }

  # if we have taken enough iterates and accepted enough parameters along each
  # model structure, then start adapting the proposal covariance matrix?

  #TODO

  setTxtProgressBar(pb, iter)
}
close(pb)

#===============================================================================
# end
#===============================================================================

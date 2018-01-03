#===============================================================================
# Calibration of PP-GPD model(s) for storm surge at Delfzijl, The Netherlands;
# Balboa, Panama; and Norfolk, Virgina, USA.
# This version does the varying record length experiments and includes the
# experiment using all of the data for each site, eliminating the need for a
# separate script.
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================
#===============================================================================
# Copyright 2017 Tony Wong
#
# MESS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# MESS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# MESS.  If not, see <http://www.gnu.org/licenses/>.
#===============================================================================

rm(list=ls())

# vvv IMPORTANT SETTINGS YOU SHOULD MODIFY, DEPENDING ON THE EXPERIMENT vvv

station <- 'delfzijl'             # can be 'delfzijl', 'balboa', or 'norfolk'
type.of.priors <- 'normalgamma'      # can be either 'uniform' or 'normalgamma'
pot.threshold <- 0.95            # GPD threshold (percentile, 0-1)

niter_mcmc_prelim000 <- 5e3      # number of MCMC iterations (PRELIMINARY chains)
nnode_mcmc_prelim000 <- 1        # number of CPUs to use (PRELIMINARY chains)
niter_mcmc_prod000 <- 5e5        # number of MCMC iterations (PRODUCTION chains)
#nnode_mcmc_prod000 <- 10          # number of CPUs to use (PRODUCTION chains)
gamma_mcmc000 <- 0.5             # speed of adaptation (0.5=faster, 1=slowest)

filename.mles <- 'surge_MLEs_ppgpd_20Dec2017.rds'  # file holding the 'deoptim.all' object with the MLEs (for initial parameters)

output.dir <- '../output/'
dat.dir <- '../data/'


if(Sys.info()['nodename']=='Tonys-MBP') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('/Users/tony/codes/EVT/R')
  nnode_mcmc_prod000 <- 1          # number of CPUs to use (PRODUCTION chains)
} else {
  # assume on Napa cluster
  machine <- 'remote'
  setwd('/home/scrim/axw322/codes/EVT/R')
  nnode_mcmc_prod000 <- 10          # number of CPUs to use (PRODUCTION chains)
}



# ^^^ IMPORTANT SETTINGS YOU SHOULD MODIFY, DEPENDING ON THE EXPERIMENT ^^^

#
#===============================================================================
# Here, set the file names for the prior distribution RDS file from the 'fit
# priors' script and the calibration data files from the processing scripts.
#===============================================================================
#

filename.priors <- paste('surge_priors_',type.of.priors,'_ppgpd_20Dec2017.rds',sep='')

if (station=='delfzijl') {
  appen <- paste('ppgpd-experiments_delfzijl_',type.of.priors,'_pot',pot.threshold*100,sep='')
  filename.datacalib <- 'tidegauge_processed_deflzijl_decl3_pot95-annual_31Dec2017.rds' # file holding the calibration data object for Delfzijl
  ind.in.mles <- 29
} else if (station=='norfolk') {
  appen <- paste('ppgpd-experiments_norfolk_',type.of.priors,'_pot',pot.threshold*100,sep='')
  filename.datacalib <- 'tidegauge_processed_norfolk_decl3-pot95-annual_31Dec2017.rds' # file holding the calibration data object for Norfolk
  ind.in.mles <- 30
} else if (station=='balboa') {
  appen <- paste('ppgpd-experiments_balboa_',type.of.priors,'_pot',pot.threshold*100,sep='')
  filename.datacalib <- 'tidegauge_processed_balboa_decl3-pot95-annual_31Dec2017.rds' # file holding the calibration data object for Balboa
  ind.in.mles <- 10
}

# Name the calibrated parameters output file
today <- Sys.Date(); today=format(today,format="%d%b%Y")
filename.parameters <- paste(output.dir,'calibratedParameters_',appen,'_',today,'.nc',sep='')

# Name the saved progress RData workspace image file
filename.mcmc <- paste(output.dir,'mcmc_',appen,'_',today,'.RData', sep='')

#
#===============================================================================
# Nothing below here should need to be modified. And really, all you ought to
# need to do is adjust the name of the station to calibrate, or the set of prior
# distributions for a supplementary experiment.
#===============================================================================
#

# On with the show!

#
#===============================================================================
# relevant libraries - do 'install.packages([library name])' if you do not have
# one yet
#===============================================================================
#

library(extRemes)
library(zoo)
library(adaptMCMC)
library(lhs)
library(DEoptim)
library(ncdf4)

#
#===============================================================================
# read and process data for temperature (auxiliary covariate for nonstationary
# parameters)
# yields: temperature_forc, time_forc
#===============================================================================
#

print('reading temperature data...')

source('read_data_temperature.R')

# maximum temperature serves as an additinoal prior constraint on kappa0, kappa1
# that is, kappa1 > -kappa0/Tmax (otherwise, kappa = kappa0 + kappa1*T might be
# negative)
Tmax <- max(temperature_forc)

print('...done.')

#
#===============================================================================
# read processed data object
#===============================================================================
#

print('reading processed tide gauge data...')

data_calib <- readRDS(paste(dat.dir,filename.datacalib,sep=''))
gpd.experiments <- names(data_calib)[intersect(which(nchar(names(data_calib))>3) , grep('gpd', names(data_calib)))]

print('...done.')

#
#===============================================================================
# set up PP-GPD model parameters
#===============================================================================
#

print('setting up PP-GPD model parameters from DE optimization...')

source('parameter_setup_dayPOT.R')
source('likelihood_ppgpd.R')

priors <- readRDS(paste(output.dir,filename.priors,sep=''))
mle.fits <- readRDS(paste(output.dir,filename.mles,sep=''))

initial.values <- vector('list', nmodel); names(initial.values) <- types.of.gpd
for (model in types.of.gpd) {initial.values[[model]] <- mle.fits[[model]][ind.in.mles,]}

print('...done.')

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
amcmc_prelim <- vector('list', length(gpd.experiments)); names(amcmc_prelim) <- gpd.experiments
for (gpd.exp in gpd.experiments) {amcmc_prelim[[gpd.exp]] <- vector('list', nmodel); names(amcmc_prelim[[gpd.exp]]) <- types.of.model}

for (gpd.exp in gpd.experiments) {
  print(paste('Starting preliminary calibration for experiment ',gpd.exp,'...', sep=''))
  for (model in types.of.gpd) {
    print(paste('Starting preliminary calibration for model ',model,' (',nnode_mcmc,' cores x ',niter_mcmc,' iterations)...', sep=''))
    if(model=='gpd3') {auxiliary <- NULL
    } else {auxiliary <- trimmed_forcing(data_calib[[gpd.exp]]$year, time_forc, temperature_forc)$temperature}
    accept_mcmc <- accept_mcmc_many + (accept_mcmc_few - accept_mcmc_many)/length(parnames_all[[model]])
    step_mcmc <- as.numeric(0.05*apply(X=mle.fits[[model]], MARGIN=2, FUN=sd))
    tbeg=proc.time()
    amcmc_prelim[[gpd.exp]][[model]] = MCMC(log_post_ppgpd, niter_mcmc, initial.values[[model]],
                              adapt=TRUE, acc.rate=accept_mcmc, scale=step_mcmc,
                              gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                              parnames=parnames_all[[model]], data_calib=data_calib[[gpd.exp]],
                              priors=priors, auxiliary=auxiliary, model=model)
    tend=proc.time()

    print(paste('... done. Took ',round(as.numeric(tend-tbeg)[3]/60,2),' minutes', sep=''))
  }
}

# save progress
##print(paste('saving preliminary results as .RData file (',filename.mcmc,') to read and use later...',sep=''))
##save.image(file=filename.mcmc)
##print('...done.')

if(FALSE){
# preliminary plot: gpd3
bounds <- rbind( c(0.007, 0.014), c(240, 550), c(-0.25, 0.45))
par(mfrow=c(5,3))
for (gpd.exp in gpd.experiments) {
  for (p in 1:3) {hist(amcmc_prelim[[gpd.exp]]$gpd3$samples[(0.5*niter_mcmc):niter_mcmc,p], xlab=parnames_all$gpd3[p], main='', xlim=bounds[p,])}
}
# preliminary plot: gpd6
bounds <- rbind( c(0, 0.015), c(-0.01, 0.02), c(4.4, 7), c(-2, 2), c(-0.7, 1), c(-1.5, 2))
par(mfrow=c(5,6))
for (gpd.exp in gpd.experiments) {
  for (p in 1:6) {hist(amcmc_prelim[[gpd.exp]]$gpd3$samples[(0.5*niter_mcmc):niter_mcmc,p], xlab=parnames_all$gpd6[p], main='', xlim=bounds[p,])}
}
# preliminary plot: all as kernel density estimates
bounds.plot <- vector('list', nmodel); names(bounds.plot) <- types.of.gpd
bounds.plot$gpd3 <- rbind( c(0.007, 0.014), c(240, 550), c(-0.25, 0.45))
bounds.plot$gpd4 <- rbind( c(0, 0.015), c(-0.01, 0.02), c(240, 550), c(-0.25, 0.45))
bounds.plot$gpd5 <- rbind( c(0, 0.015), c(-0.01, 0.02), c(4.4, 7), c(-2, 2), c(-0.25, 0.45))
ylims <- vector('list', nmodel); names(ylims) <- types.of.gpd

# set ylims for each of these, probably based on values in one of the KDEs

par(mfrow=c(4,6))
model <- 'gpd3'; nnode <- 512
colors <- c('black','brown','goldenrod','orange','red')
for (p in 1:length(parnames_all[[model]])) {
    for (gpd.exp in gpd.experiments) {
        tmp <- density(x=amcmc_prelim[[gpd.exp]][[model]]$samples[(0.5*niter_mcmc):niter_mcmc,p], n=nnode, from=bounds[p,1], to=bounds[p,2])
        if(gpd.exp == 'gpd30') {
            plot(tmp$x, tmp$y, xlab=parnames_all[[model]][p], ylab='Density', main='', xlim=bounds[p,], col=colors[1], type='l', lwd=2, ylim=c(0, 2*max(tmp$y)))
        } else {
            lines(tmp$x, tmp$y, xlab=parnames_all[[model]][p], ylab='Density', main='', xlim=bounds[p,], col=colors[match(gpd.exp,gpd.experiments)], lwd=2)
        }
    }
    plot.new()
}
model <- 'gpd4'
for (p in 1:length(parnames_all[[model]])) {
    for (gpd.exp in gpd.experiments) {
        tmp <- density(x=amcmc_prelim[[gpd.exp]][[model]]$samples[3e4:5e4,p], n=nnode, from=bounds[p,1], to=bounds[p,2])
        if(gpd.exp == 'gpd30') {
            plot(tmp$x, tmp$y, xlab=parnames_all[[model]][p], ylab='Density', main='', xlim=bounds[p,], col=colors[1], type='l', lwd=2, ylim=c(0, 4.5*max(tmp$y)))
        } else {
            lines(tmp$x, tmp$y, xlab=parnames_all[[model]][p], ylab='Density', main='', xlim=bounds[p,], col=colors[match(gpd.exp,gpd.experiments)], lwd=2)
        }
    }
    if(p > 2) {plot.new()}
}
model <- 'gpd5'
for (p in 1:length(parnames_all[[model]])) {
    for (gpd.exp in gpd.experiments) {
        tmp <- density(x=amcmc_prelim[[gpd.exp]][[model]]$samples[3e4:5e4,p], n=nnode, from=bounds[p,1], to=bounds[p,2])
        if(gpd.exp == 'gpd30') {
            plot(tmp$x, tmp$y, xlab=parnames_all[[model]][p], ylab='Density', main='', xlim=bounds[p,], col=colors[1], type='l', lwd=2, ylim=c(0, 4.5*max(tmp$y)))
        } else {
            lines(tmp$x, tmp$y, xlab=parnames_all[[model]][p], ylab='Density', main='', xlim=bounds[p,], col=colors[match(gpd.exp,gpd.experiments)], lwd=2)
        }
    }
    if(p==5) {plot.new()}
}
model <- 'gpd6'; bounds <- rbind( c(0, 0.015), c(-0.01, 0.02), c(4.4, 7), c(-2, 2), c(-0.7, 1), c(-1.5, 2))
for (p in 1:length(parnames_all[[model]])) {
    for (gpd.exp in gpd.experiments) {
        tmp <- density(x=amcmc_prelim[[gpd.exp]][[model]]$samples[3e4:5e4,p], n=nnode, from=bounds[p,1], to=bounds[p,2])
        if(gpd.exp == 'gpd30') {
            plot(tmp$x, tmp$y, xlab=parnames_all[[model]][p], ylab='Density', main='', xlim=bounds[p,], col=colors[1], type='l', lwd=2, ylim=c(0, 4.5*max(tmp$y)))
        } else {
            lines(tmp$x, tmp$y, xlab=parnames_all[[model]][p], ylab='Density', main='', xlim=bounds[p,], col=colors[match(gpd.exp,gpd.experiments)], lwd=2)
        }
    }
}

} # end preliminary plots

#
#===============================================================================
# set up and run PRODUCTION MCMC calibration
#===============================================================================
#

# then use these initial estimates of step_mcmc to launch the parallel chains
# (from amcmc_prelim[[gpd.exp]][[model]]$cov.jump)

nnode_mcmc <- nnode_mcmc_prod000
niter_mcmc <- niter_mcmc_prod000
gamma_mcmc <- gamma_mcmc000
startadapt_mcmc <- max(500,round(0.05*niter_mcmc))
stopadapt_mcmc <- round(niter_mcmc*1.0)
accept_mcmc_few <- 0.44         # optimal for only one parameter
accept_mcmc_many <- 0.234       # optimal for many parameters
amcmc_out <- vector('list', length(gpd.experiments)); names(amcmc_out) <- gpd.experiments
for (gpd.exp in gpd.experiments) {amcmc_out[[gpd.exp]] <- vector('list', nmodel); names(amcmc_out[[gpd.exp]]) <- types.of.model}

for (gpd.exp in gpd.experiments) {
  print(paste('Starting production calibration for experiment ',gpd.exp,'...', sep=''))
  for (model in types.of.model) {
    print(paste('Starting production calibration for model ',model,' (',nnode_mcmc,' cores x ',niter_mcmc,' iterations)...', sep=''))

    if (model %in% types.of.gpd) {
      initial_parameters <- amcmc_prelim[[gpd.exp]][[model]]$samples[amcmc_prelim[[gpd.exp]][[model]]$n.sample,]
      if(model=='gpd3') {auxiliary <- NULL
      } else {auxiliary <- trimmed_forcing(data_calib[[gpd.exp]]$year, time_forc, temperature_forc)$temperature}
      accept_mcmc <- accept_mcmc_many + (accept_mcmc_few - accept_mcmc_many)/length(parnames_all[[model]])
      step_mcmc <- amcmc_prelim[[gpd.exp]][[model]]$cov.jump
      if(nnode_mcmc==1) {
        # do single chain
        tbeg=proc.time()
        amcmc_out[[gpd.exp]][[model]] = MCMC(log_post_ppgpd, niter_mcmc, initial_parameters,
                               adapt=TRUE, acc.rate=accept_mcmc, scale=step_mcmc,
                               gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                               parnames=parnames_all[[model]], data_calib=data_calib[[gpd.exp]],
                               priors=priors, auxiliary=auxiliary, model=model)
        tend=proc.time()
      } else if(nnode_mcmc > 1) {
        # do parallel chains
        tbeg <- proc.time()
        amcmc_out[[gpd.exp]][[model]] <- MCMC.parallel(log_post_ppgpd, niter_mcmc, initial_parameters,
                               n.chain=nnode_mcmc, n.cpu=nnode_mcmc, packages='extRemes',
                               scale=step_mcmc, adapt=TRUE, acc.rate=accept_mcmc,
                               gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                               parnames=parnames_all[[model]], data_calib=data_calib[[gpd.exp]],
                               priors=priors, auxiliary=auxiliary, model=model)
        tend <- proc.time()
      }
    } else {print('error - unknown model type')}
    print(paste('... done. Took ',round(as.numeric(tend-tbeg)[3]/60,2),' minutes', sep=''))
  }
##  print(paste('saving production results (so far) as .RData file (',filename.mcmc,') to read and use later...',sep=''))
##  save.image(file=filename.mcmc)
##  print('...done.')
}

# test plot
if (FALSE) {
model <- 'gpd3'; gpd.exp <- 'gpd90'
par(mfrow=c(3,2))
for (p in 1:length(parnames_all[[model]])) {
    plot(amcmc_out[[gpd.exp]][[model]][[1]]$samples[,p], type='l', ylab=parnames_all[[model]][p])
    hist(amcmc_out[[gpd.exp]][[model]][[1]]$samples[round(0.5*niter_mcmc):niter_mcmc,p], xlab=parnames_all[[model]][p], main='')
}
}

#
#===============================================================================
# convergence diagnostics
#===============================================================================
#

# Gelman and Rubin diagnostics - determine and chop off for burn-in
niter.test <- seq(from=round(0.1*niter_mcmc), to=niter_mcmc, by=round(0.05*niter_mcmc))
gr.test <- vector('list', length(gpd.experiments)); names(gr.test) <- gpd.experiments
for (gpd.exp in gpd.experiments) {
  gr.test[[gpd.exp]] <- mat.or.vec(length(niter.test), nmodel)
  colnames(gr.test[[gpd.exp]]) <- types.of.model
}
gr.tmp <- rep(NA, length(niter.test))

for (gpd.exp in gpd.experiments) {
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
          eval(parse(text=paste('mcmc',m,' <- as.mcmc(amcmc_out[[gpd.exp]][[model]][[m]]$samples[1:niter.test[i],])', sep='')))
        }
        eval(parse(text=paste('mcmc_chain_list = mcmc.list(list(', string.mcmc.list , '))', sep='')))

        gr.test[[gpd.exp]][i,model] <- as.numeric(gelman.diag(mcmc_chain_list)[2])
      }
    } else {print('error - nnode_mcmc < 1 makes no sense')}
  }
}

# Monitor posterior 5, 50 and 95% quantiles for drift
# Only checking for one of the chains
quant <- vector('list', length(gpd.experiments)); names(quant) <- gpd.experiments
for (gpd.exp in gpd.experiments) {quant[[gpd.exp]] <- vector('list', nmodel); names(quant[[gpd.exp]]) <- types.of.model}
names.monitor <- c('q05', 'q50', 'q95')
for (gpd.exp in gpd.experiments) {
  for (model in types.of.model) {
    quant[[gpd.exp]][[model]] <- vector('list', 3); names(quant[[gpd.exp]][[model]]) <- names.monitor
    for (q in names.monitor) {
      quant[[gpd.exp]][[model]][[q]] <- mat.or.vec(length(niter.test)-1, length(parnames_all[[model]]))
      for (i in 1:(length(niter.test)-1)) {
        if(nnode_mcmc==1) {
          quant[[gpd.exp]][[model]][[q]][i,] <- apply(X=amcmc_out[[gpd.exp]][[model]]$samples[niter.test[i]:niter_mcmc,], MARGIN=2, FUN=quantile, probs=as.numeric(substr(q, 2,3))*0.01)
        } else {
          quant[[gpd.exp]][[model]][[q]][i,] <- apply(X=amcmc_out[[gpd.exp]][[model]][[1]]$samples[niter.test[i]:niter_mcmc,], MARGIN=2, FUN=quantile, probs=as.numeric(substr(q, 2,3))*0.01)
        }
      }
    }
  }
}

if(FALSE) {
# examples monitoring of stability of quantiles:
model <- 'gpd3'; gpd.exp <- 'gpd30'
par(mfrow=c(3,2))
for (p in 1:length(parnames_all[[model]])) {
  ran <- max(quant[[gpd.exp]][[model]]$q95[,p])-min(quant[[gpd.exp]][[model]]$q05[,p])
  lb <- min(quant[[gpd.exp]][[model]]$q05[,p]) - 0.05*ran; ub <- max(quant[[gpd.exp]][[model]]$q95[,p]) + 0.05*ran
  plot(niter.test[1:(length(niter.test)-1)], quant[[gpd.exp]][[model]]$q50[,p], type='l',
    ylim=c(lb,ub), ylab=parnames_all[[model]][p], xlab='From HERE to end of chain')
  lines(niter.test[1:(length(niter.test)-1)], quant[[gpd.exp]][[model]]$q05[,p], lty=2); lines(niter.test[1:(length(niter.test)-1)], quant[[gpd.exp]][[model]]$q95[,p], lty=2);
}
# the thing to note is that these stabilize as you include more members (i.e.,
# as you move from right to left)
}


# Heidelberger and Welch diagnostics?
hw.diag <- vector('list', length(gpd.experiments)); names(hw.diag) <- gpd.experiments
for (gpd.exp in gpd.experiments) {
  hw.diag[[gpd.exp]] <- vector('list', nmodel); names(hw.diag[[gpd.exp]]) <- types.of.model
  for (model in types.of.model) {
    hw.diag[[gpd.exp]][[model]] <- heidel.diag(as.mcmc(amcmc_out[[gpd.exp]][[model]][[1]]$samples), eps=0.1, pvalue=0.05)
  }
}
# 30, 50 and 70-year experiments might need more iterations in order to pass
# these diagnostics. can do that, but the ensemble statistics do not change
# much (note the quantile monitoring results). so go with the shorter ensembles
# for a fair comparison with the longer record experiments.

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
# save a separate ifirst for each experiment
ifirst <- rep(NA, length(gpd.experiments)); names(ifirst) <- gpd.experiments
for (gpd.exp in gpd.experiments) {
  if(nnode_mcmc==1) {
    ifirst[[gpd.exp]] <- round(0.5*niter_mcmc)
  } else {
    gr.max <- 1.1
    lgr <- rep(NA, length(niter.test))
    for (i in 1:length(niter.test)) {lgr[i] <- all(gr.test[[gpd.exp]][i,] < gr.max)}
    for (i in seq(from=length(niter.test), to=1, by=-1)) {
      if( all(lgr[i:length(lgr)]) ) {ifirst[[gpd.exp]] <- niter.test[i]}
    }
  }
}

chains_burned <- vector('list', length(gpd.experiments)); names(chains_burned) <- gpd.experiments
for (gpd.exp in gpd.experiments) {
  chains_burned[[gpd.exp]] <- vector('list', nmodel); names(chains_burned[[gpd.exp]]) <- types.of.model
  for (model in types.of.model) {
    if(nnode_mcmc > 1) {
      chains_burned[[gpd.exp]][[model]] <- vector('list', nnode_mcmc)
      for (m in 1:nnode_mcmc) {
        chains_burned[[gpd.exp]][[model]][[m]] <- amcmc_out[[gpd.exp]][[model]][[m]]$samples[(ifirst[[gpd.exp]]+1):niter_mcmc,]
      }
    } else {
      chains_burned[[gpd.exp]][[model]] <- amcmc_out[[gpd.exp]][[model]]$samples[(ifirst[[gpd.exp]]+1):niter_mcmc,]
    }
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

n.sample <- 10000

for (gpd.exp in gpd.experiments) {
  for (model in types.of.model) {
    if(nnode_mcmc == 1) {
      ind.sample <- sample(x=1:nrow(chains_burned[[gpd.exp]][[model]]), size=n.sample, replace=FALSE)
      chains_burned_thinned[[gpd.exp]][[model]] <- chains_burned[[gpd.exp]][[model]][ind.sample,]
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
        ind.sample <- sample(x=1:nrow(chains_burned[[gpd.exp]][[model]][[m]]), size=n.sample.sub[m], replace=FALSE)
        chains_burned_thinned[[gpd.exp]][[model]][[m]] <- chains_burned[[gpd.exp]][[model]][[m]][ind.sample,]
      }
    }
  }
}

}#====================================


# Combine all of the chains from 'ifirst' to 'niter_mcmc' into a potpourri of
# [alleged] samples from the posterior. Only saving the transition covariance
# matrix for one of the chains (if in parallel).
parameters.posterior <- vector('list', length(gpd.experiments)); names(parameters.posterior) <- gpd.experiments
covjump.posterior <- vector('list', length(gpd.experiments)); names(covjump.posterior) <- gpd.experiments
for (gpd.exp in gpd.experiments) {
  parameters.posterior[[gpd.exp]] <- vector('list', nmodel); names(parameters.posterior[[gpd.exp]]) <- types.of.model
  covjump.posterior[[gpd.exp]]    <- vector('list', nmodel); names(covjump.posterior[[gpd.exp]])    <- types.of.model
  for (model in types.of.model) {
    if(nnode_mcmc==1) {
      parameters.posterior[[gpd.exp]][[model]] <- chains_burned_thinned[[gpd.exp]][[model]]
      covjump.posterior[[gpd.exp]][[model]] <- amcmc_out[[gpd.exp]][[model]]$cov.jump
    } else {
      parameters.posterior[[gpd.exp]][[model]] <- chains_burned_thinned[[gpd.exp]][[model]][[1]]
      covjump.posterior[[gpd.exp]][[model]] <- amcmc_out[[gpd.exp]][[model]][[1]]$cov.jump
      for (m in 2:nnode_mcmc) {
        parameters.posterior[[gpd.exp]][[model]] <- rbind(parameters.posterior[[gpd.exp]][[model]], chains_burned_thinned[[gpd.exp]][[model]][[m]])
      }
    }
  }
}

# save results in case you need to revisit later
##print(paste('saving MCMC results as .RData file (',filename.mcmc,') to read and use later...',sep=''))
##save.image(file=filename.mcmc)
##save(list=c('amcmc_out','ifirst'), file=filename.mcmc) # just do this later
##print('...done.')

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
var.parnames   <- vector('list', nmodel); names(var.parnames)   <- types.of.model
var.parameters <- vector('list', length(gpd.experiments)); names(var.parameters) <- gpd.experiments
var.covjump    <- vector('list', length(gpd.experiments)); names(var.covjump)    <- gpd.experiments
for (gpd.exp in gpd.experiments) {
  var.parameters[[gpd.exp]] <- vector('list', nmodel); names(var.parameters[[gpd.exp]]) <- types.of.model
  var.covjump[[gpd.exp]]    <- vector('list', nmodel); names(var.covjump[[gpd.exp]])    <- types.of.model
}
dim.ensemble   <- vector('list', nmodel); names(dim.ensemble)   <- types.of.model
dim.name <- ncdim_def('name.len', '', 1:lmax, unlim=FALSE)
dim.time <- ncdim_def('ntime', '', (time_forc), unlim=FALSE)
var.time <- ncvar_def('time', '', dim.time, -999)
var.temperature <- ncvar_def('temperature', 'Global mean temperature relative to 1901-2000', dim.time, -999)
for (model in types.of.model) {
  dim.parameters[[model]] <- ncdim_def(paste('n.parameters.',model,sep=''), '', 1:length(parnames_all[[model]]), unlim=FALSE)
  dim.ensemble[[model]]   <- ncdim_def(paste('n.ensemble.',model,sep=''), 'ensemble member', 1:nrow(parameters.posterior[[gpd.exp]][[model]]), unlim=FALSE)
  var.parnames[[model]]   <- ncvar_def(paste('parnames.',model,sep=''), '', list(dim.name,dim.parameters[[model]]), prec='char')
  for (gpd.exp in gpd.experiments) {
    var.parameters[[gpd.exp]][[model]] <- ncvar_def(paste('parameters.',gpd.exp,'.',model,sep=''), '', list(dim.parameters[[model]],dim.ensemble[[model]]), -999)
    var.covjump[[gpd.exp]][[model]] <- ncvar_def(paste('covjump.',gpd.exp,'.',model,sep=''), '', list(dim.parameters[[model]],dim.parameters[[model]]), -999)
  }
}
# length(gpd.experiments) * nmodel * 2 accounts for sizes of var.parmaeters and var.covjump
output.to.file <- vector('list', (length(gpd.experiments)*nmodel*2) + length(var.parnames) + 2)
output.to.file[[1]] <- var.time
output.to.file[[2]] <- var.temperature
# this counter will keep track of how many items are on the output.to.file list so far
cnt <- 3
for (model in types.of.model) {
  output.to.file[[cnt]] <- var.parnames[[model]]
  cnt <- cnt + 1
  for (gpd.exp in gpd.experiments) {
    output.to.file[[cnt]] <- var.parameters[[gpd.exp]][[model]]
    output.to.file[[cnt+1]] <- var.covjump[[gpd.exp]][[model]]
    cnt <- cnt + 2
  }
}
#outnc <- nc_create(filename.parameters, list(var.parameters, var.parnames, var.covjump))
outnc <- nc_create(filename.parameters, output.to.file)
ncvar_put(outnc, var.time, time_forc)
ncvar_put(outnc, var.temperature, temperature_forc)
for (model in types.of.model) {
  ncvar_put(outnc, var.parnames[[model]], parnames_all[[model]])
  for (gpd.exp in gpd.experiments) {
    ncvar_put(outnc, var.parameters[[gpd.exp]][[model]], t(parameters.posterior[[gpd.exp]][[model]]))
    ncvar_put(outnc, var.covjump[[gpd.exp]][[model]], covjump.posterior[[gpd.exp]][[model]])
  }
}
nc_close(outnc)

# save results in case you need to revisit later
print(paste('saving MCMC results as .RData file (',filename.mcmc,') to read and use later...',sep=''))
##save.image(file=filename.mcmc)
save(list=c('amcmc_out','ifirst','data_calib','priors'), file=filename.mcmc)
print('...done.')

#
#===============================================================================
# End
#===============================================================================
#

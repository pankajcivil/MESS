#===============================================================================
# This file contains the working script for estimating the          #
# marginal likelihood of different GPD surge models and data        #
# experiments.                                                      #
#                                                                   #
# This revised version calculates all likelihood estimates as a parallel
# 'foreach' loop. (TODO)
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

# import libraries
library(mvtnorm)  # used for multivariate normal samples and densities
#library(coda)     # used to compute spectral density of MCMC draws at frequency 0
library(extRemes) # used to compute GEV densities within posterior likelihood function
library(ncdf4)

if(Sys.info()['user']=='tony') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('/Users/tony/codes/EVT/R')
  path.R <- '/Users/tony/codes/EVT/R'
  # set data and save directories
  path.data <- '/Users/tony/codes/EVT/output'
  path.save <- '/Users/tony/codes/EVT/output/bma'
  nnode <- 2          # number of CPUs to use
} else {
  # assume on Napa cluster
  machine <- 'remote'
  setwd('/home/scrim/axw322/codes/EVT/R')
  path.R <- '/home/scrim/axw322/codes/EVT/R'
  path.data <- '/home/scrim/axw322/codes/EVT/output'
  path.save <- '/home/scrim/axw322/codes/EVT/output/bma'
  nnode <- 10          # number of CPUs to use
}

# import file containing the log likelihood calculations
source(paste(path.R,'likelihood_ppgpd.R',sep='/'))

source(paste(path.R,'read_data_temperature.R',sep='/'))

# set up table of experiment parameters
# some of these combinations (in terms of the data length and sites) don't exist, but they
# just throw errors anyway. There is probably a better way to construct this data frame.
# could include the type of prior in this as well

sites <- c('delfzijl','balboa','norfolk')
n_sites <- length(sites)
gpd.models <- c('gpd3','gpd4','gpd5','gpd6')
n_model <- length(gpd.models)

data.lengths <- vector('list', n_sites); names(data.lengths) <- sites
data.lengths$norfolk <- c('30','50','70','89')
data.lengths$balboa <- c('30','50','70','90','107')
data.lengths$delfzijl <- c('30','50','70','90','110','137')

# list to store the actual output
output <- vector('list', n_sites); names(output) <- sites
for (site in sites) {
  output[[site]] <- vector('list', n_model)
  names(output[[site]]) <- gpd.models
  for (model in gpd.models) {
    output[[site]][[model]] <- vector('list', length(data.lengths[[site]]))
    names(output[[site]][[model]]) <- data.lengths[[site]]
  }
}

# data frame to store the experimental details
experiments <- expand.grid(station = c('delfzijl','balboa','norfolk'),
                           gpd.model=c('gpd3','gpd4','gpd5','gpd6'),
                           data.length=c('30','50','70','89','90','107','110','137'))
irem <- NULL
irem <- c(irem, which((experiments[,'data.length']==89 | experiments[,'data.length']==107) & experiments[,'station']=='delfzijl'))
irem <- c(irem, which((experiments[,'data.length']==89 | experiments[,'data.length']==137) & experiments[,'station']=='balboa'))
irem <- c(irem, which((experiments[,'data.length']==107 | experiments[,'data.length']==137) & experiments[,'station']=='norfolk'))
n_experiments <- nrow(experiments)
output <- vector('list', n_experiments)

source('bridge_sample_functions.R')

for (ee in 1:n_experiments) {

  # get parameters for this particular experiment
  print(experiments[ee,])
  station <- experiments[ee,'station']
  gpd.model <- experiments[ee,'gpd.model']
  data.length <- experiments[ee,'data.length']

  # set output (saved as .RData; to be collected into a single output file later) file name
  filename.out <- paste(path.save,'/ml_',station,'_',gpd.model,'_',data.length,'.RData',sep='')
  if (file.exists(filename.out)) {
    #stop('Output file already exists!')
    print('Output file already exists - skipping...')
  } else {

    # read in calibration output file
    print('loading calibration file...')

    type.of.priors <- 'normalgamma'     # can be 'uniform' or 'normalgamma'
    # use this if multiple files exist for the same location and prior
    #calib_date <- '28Jul2017'
    #setwd(path.data)
    if (exists('calib_date')) {
      filename.calib <- paste(path.data,'/everything_mcmc_ppgpd-experiments_',station,'_',type.of.priors,'_',calib_date,'.RData',sep='')
    } else {
      filename.calib <- Sys.glob(paste(path.data,'/everything_mcmc_ppgpd-experiments_',station,'_',type.of.priors,'_*','.RData',sep=''))
    }
    load(filename.calib)

    print('done!')

    # set name for data experiment from data.length variable
    gpd.exp <- paste('gpd',data.length,sep='')

    # compute the burn-in length for sampling based on the chains_burned object
    full.length <- nrow(amcmc_out[[gpd.exp]][[gpd.model]][[1]]$samples)
    burned.length <- nrow(chains_burned[[gpd.exp]][[gpd.model]][[1]])
    burn.in <- full.length - burned.length

    # set number of samples to use for estimate
    post.samp.num <- 50000
    imp.samp.num <- 50000

    # burn in samples and log.p values
    post.samples <- amcmc_out[[gpd.exp]][[gpd.model]][[1]]$samples
    post.samples <- post.samples[(burn.in+1):full.length,]
    post.ll <- amcmc_out[[gpd.exp]][[gpd.model]][[1]]$log.p
    post.ll <- post.ll[(burn.in+1):full.length]

    # fit normal approximation to the posteriors
    post.mean <- colMeans(post.samples)
    post.cov <- cov(post.samples)

    # get posterior samples
    print('sampling from posterior distribution...')

    samp.names <- c('samples','log.imp','log.p')
    post.samp <- setNames(vector("list",length(samp.names)),samp.names)
    samp.idx <- sample(x=nrow(post.samples), size=post.samp.num, replace=TRUE)
    post.samp$samples <- post.samples
    post.samp$samples <- post.samples[samp.idx,]
    # get posterior log-likelihood of sampled posterior values
    post.samp$log.p <- post.ll
    post.samp$log.p <- post.ll[samp.idx]
    # get importance log-likelhood of posterior samples
    post.samp$log.imp <- dmvnorm(x=post.samp$samples, mean=post.mean, sigma=post.cov, log=TRUE)

    print('done!')

    # get importance samples and likelihood
    print('sampling from importance distribution...')

    imp.samp <- setNames(vector("list",length(samp.names)),samp.names)
    imp.samp$samples <- rmvnorm(n=imp.samp.num, mean=post.mean, sigma=post.cov)
    imp.samp$log.imp <- dmvnorm(x=imp.samp$samples, mean=post.mean, sigma=post.cov, log=TRUE)
    colnames(imp.samp$samples) <- colnames(post.samp$samples)
    # compute posterior log-likelihood of importance samples

    # set auxiliary parameters as trimmed_forcing for the relevant model
    if (gpd.model == 'gpd3') {aux <- NULL
    } else {aux <-  trimmed_forcing(data_calib[[gpd.exp]]$year, time_forc, temperature_forc)$forcing}

    imp.samp$log.p <- apply(imp.samp$samples, 1, log_post_ppgpd,
                            parnames=colnames(post.samp$samples),
                            priors=priors,
                            data_calib=data_calib[[gpd.exp]],
                            model = gpd.model,
                            auxiliary=aux)

    print('done!')

    print('beginning bridge sampling recursion...')

    # set tolerance for halting of iteration
    TOL <- 1e-10

    # initialize storage for estimates
    ml <- mat.or.vec(nr=1,nc=1)

    # initialize with starting value
    # we can't quite start with the reciprocal importance sampling estimate from
    # Gelfand and Dey (1994) due to numerics (we get 0 values when we exponentiate
    # the difference of the importance log-densities and posterior log-likelihoods), so we just
    # average the ratios on a log scale.
    ml[1] <- -mean(post.samp$log.imp - post.samp$log.p)
    ml[2] <- bridge.samp.iter(ml[1], post.samp[c('log.p','log.imp')], imp.samp[c('log.p','log.imp')])

    # iterate until within tolerance.
    t <- 2
    while (abs(ml[t] - ml[t-1]) >= TOL) {
      ml[t+1] <- bridge.samp.iter(ml[t], post.samp[c('log.p', 'log.imp')], imp.samp[c('log.p', 'log.imp')])
      t <- t+1
    }

    print('done!')

    print('computing relative standard error of estimate')

    # compute the relative standard error of the bridge sampling estimator
    # we can treat the posterior samples as iid due to re-sampling from the posterior,
    # so we use the error formula from Fruhwirth-Schnatter (2004) with the spectral density
    # at frequency 0 set equal to 1.

    re.sq <- bridge.samp.rel.err(ml[length(ml)], post.samp[c('log.p','log.imp')], imp.samp[c('log.p','log.imp')])

    # save result of run
    # if save directory doesn't exist, create it
    #ifelse(!dir.exists(path.save), dir.create(path.save), FALSE)
    #setwd(path.save)

    save(list=c('post.samp','imp.samp', 'ml', 're.sq', 'station', 'data.length', 'gpd.model'), file=filename.out)
  }
}
#data_many <- finalOutput
#names(data_many) <- names(data_set)


#===============================================================================
# end
#===============================================================================

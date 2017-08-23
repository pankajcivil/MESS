#####################################################################
# This file contains the working script for estimating the          #
# marginal likelihood of different GPD surge models and data        #
# experiments.                                                      #
#                                                                   #
# It is assumed that each likelihood estimate will be done on       #
# separate cluster nodes using a PBS array. This can be changed to  #
# an MPI job using the cluster package relatively simply.           #
#####################################################################

# import libraries
library(mvtnorm)  # used for multivariate normal samples and densities
#library(coda)     # used to compute spectral density of MCMC draws at frequency 0
library(extRemes) # used to compute GEV densities within posterior likelihood function

# import file containing the log likelihood calculations
path.R <- '/storage/home/vxs914/work/MESS/R'
source(paste(path.R,'likelihood_ppgpd.R',sep='/'))

# set data and save directories
path.data <- '/storage/home/vxs914/scratch/BMA'
path.save <- '/storage/home/vxs914/work/MESS/ml'

# if running from the command line, get arguments specifying the station, 
# model, and data length; otherwise, pass those in order out of the options 
# in the experiment table below
args = commandArgs(trailingOnly=TRUE)
# if no arguments are passed, it is assumed this is being run in parallel using PBS arrays
if (length(args) == 0) {
  # set up table of experiment parameters
  # some of these combinations (in terms of the data length and sites) don't exist, but they
  # just throw errors anyway. There is probably a better way to construct this data frame.
  # could include the type of prior in this as well
  experiments <- expand.grid(station = c('delfzijl','balboa','norfolk'),
                             gpd.model=c('gpd3','gpd4','gpd5','gpd6'),
                             data.length=c('30','50','70','89','90','107','110','137'))
  # get the array ID
  AI <- Sys.getenv("PBS_ARRAYID")
  NAI <- as.numeric(AI)
  
  # get parameters for this particular experiment
  print(experiments[NAI,])
  station <- experiments[NAI,'station']
  gpd.model <- experiments[NAI,'gpd.model']
  data.length <- experiments[NAI,'data.length']
  
} else {
  station <- args[1]
  gpd.model <- args[2]
  data.length <- args[3]
}

# set output (saved as .RData; to be collected into a single output file later) file name
filename.out <- paste('ml_',station,'_',gpd.model,'_',data.length,'.RData',sep='')
if (file.exists(paste(path.save, filename.out, sep='/'))) {
   stop('Output file already exists!')
}


# read in calibration output file
print('loading calibration file...')

type.of.priors <- 'normalgamma'     # can be 'uniform' or 'normalgamma'
# use this if multiple files exist for the same location and prior
#calib_date <- '28Jul2017'
setwd(path.data)
if (exists('calib_date')) {
  filename.calib <- paste('everything_mcmc_ppgpd-experiments_',station,'_',type.of.priors,'_',calib_date,'.RData',sep='')
} else {
  filename.calib <- Sys.glob(paste('everything_mcmc_ppgpd-experiments_',station,'_',type.of.priors,'_*','.RData',sep=''))
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
post.samp.num <- 10000
imp.samp.num <- 10000

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
} else {aux <-  trimmed_forcing(data_calib[[gpd.exp]]$year, time_forc, temperature_forc)$temperature}

imp.samp$log.p <- apply(imp.samp$samples, 1, log_post_ppgpd,
                        parnames=colnames(post.samp$samples),
                        priors=priors,
                        data_calib=data_calib[[gpd.exp]],
                        model = gpd.model,
                        auxiliary=aux)

print('done!')
  
# compute the bridge sampling estimate using the iterative procedure from Meng and Wong (1996)
# starting value is not quite the reciprocal importance sampling estimate 
# from Gelfand and Dey (1994)
# due to infinite values obtained when exponentiating, we start with the reciprocal mean of
recip.imp.samp <- function(log.p,log.imp) {
  log.ratio <- log.imp - log.p
  -log(mean(exp(log.ratio)))
}


print('beginning bridge sampling recursion...')

# function to update the bridge sampling estimate at each iteration
# norm.const is the log normalizing constant estimate from the previous iteration
# post, imp are lists with log.p and log.imp passing the associated log-likelihoods
bridge.samp.iter <- function(log.norm.const,
                        post, 
                        imp) {
  
  # normalize posterior likelihoods based on previous normalizing constant estimate
  # some samples (mostly importance samples) might have infinite posterior log-likelihoods
  post.log.p.norm <- post$log.p[is.finite(post$log.p)] - log.norm.const
  imp.log.p.norm <- imp$log.p[is.finite(imp$log.p)] - log.norm.const
  
  post.log.imp <- post$log.imp[is.finite(post$log.p)]
  imp.log.imp <- imp$log.imp[is.finite(imp$log.p)]

  # get number of samples
  post.num <- length(post.log.p.norm)
  imp.num <- length(imp.log.p.norm)
  
  # compute updated estimate numerator and denominator
  imp.mean <- mean(exp(imp.log.p.norm)/(imp.num*exp(imp.log.imp)+post.num*exp(imp.log.p.norm)))
  post.mean <- mean(exp(post.log.imp)/(imp.num*exp(post.log.imp)+post.num*exp(post.log.p.norm)))
  
  # return updated estimate
  log.norm.const + log(imp.mean) - log(post.mean)
}

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
bridge.samp.rel.err <- function(log.norm.const,
                                post,
                                imp) {
  
  # normalize posterior likelihoods based on previous normalizing constant estimate
  # some samples (mostly importance samples) might have infinite posterior log-likelihoods
  post.log.p.norm <- post$log.p[is.finite(post$log.p)] - log.norm.const
  imp.log.p.norm <- imp$log.p[is.finite(imp$log.p)] - log.norm.const
  
  post.log.imp <- post$log.imp[is.finite(post$log.p)]
  imp.log.imp <- imp$log.imp[is.finite(imp$log.p)]

  # get number of samples
  post.num <- length(post.log.p.norm)
  imp.num <- length(imp.log.p.norm)
  imp.factor <- imp.num/(imp.num+post.num)
  post.factor <- post.num/(imp.num+post.num)
    
  # compute bridging function process estimates for both sequences
  post.bridge.int <- exp(post.log.imp)/(imp.factor*exp(post.log.imp)+post.factor*(exp(post.log.p.norm)))
  imp.bridge.int <- exp(imp.log.p.norm)/(imp.factor*exp(imp.log.imp)+post.factor*(exp(imp.log.p.norm)))
  
  # return squared relative error estimate
  (var(imp.bridge.int)/mean(imp.bridge.int)^2)/imp.num + (var(post.bridge.int)/mean(post.bridge.int)^2)/post.num
}

re.sq <- bridge.samp.rel.err(ml[length(ml)], post.samp[c('log.p','log.imp')], imp.samp[c('log.p','log.imp')])

# save result of run
# if save directory doesn't exist, create it
ifelse(!dir.exists(path.save), dir.create(path.save), FALSE)
setwd(path.save)

save(list=c('post.samp','imp.samp', 'ml', 're.sq', 'station', 'data.length', 'gpd.model'), file=filename.out)
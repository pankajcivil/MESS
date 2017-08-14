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
library(coda)     # used to compute spectral density of MCMC draws at frequency 0

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
  # could include the type of prior in this as well
  experiments <- expand.grid(station = c('delfzijl','balboa','norfolk'),
                             gpd.model=c('gpd3','gpd4','gpd5','gpd6'),
                             data.length=c('30','50','70','90','110','130'))
  # get the array ID
  AI <- Sys.getenv("PBS_ARRAYID")
  NAI <- as.numeric(AI)
  
  # get parameters for this particular experiment
  station <- experiments[NAI]$station
  gpd.model <- experiments[NAI]$gpd.model
  data.length <- experiments[NAI]$data.length
  
} else {
  station <- args[1]
  gpd.model <- args[2]
  data.length <- args[3]
}

# set number of posterior and importance samples
imp.samp.num <- 5000
post.samp.num <- 5000

# read in calibration output file
type.of.priors <- 'normalgamma'     # can be 'uniform' or 'normalgamma'
calib_date <- '28Jul2017'
setwd(path.data)
filename.calib <- paste('everything_mcmc_ppgpd-experiments_',station,'_',type.of.priors,'_',calib_date,'.RData',sep='')
load(filename.calib)

# set name for data experiment from data.length variable
gpd.exp <- paste('gpd',data.length,sep='')

# compute the burn-in length for sampling based on the chains_burned object
full.length <- dim(amcmc_out[[gpd.exp]][[gpd.model]][[1]]$samples)[1]
burned.length <- dim(chains_burned[[gpd.exp]][[gpd.model]][[1]])[1]
burn.in <- full.length - burned.length
# burn in samples and log.p values
post.samples <- amcmc_out[[gpd.exp]][[gpd.model]][[1]]$samples
post.samples <- post.samples[(burn.in+1):full.length,]
post.ll <- amcmc_out[[gpd.exp]][[gpd.model]][[1]]$log.p
post.ll <- post.ll[(burn.in+1):full.length]

# fit normal approximation to the posteriors
post.means <- colMeans(post.samples)
post.cov <- cov(post.samples)

# get normal approximation samples and likelihood
print('sampling from importance distribution...')

samp.names <- c('samples','log.imp','log.p')
imp.samp <- setNames(vector("list",length(samp.names)),samp.names)
imp.samp$samples <- rmvnorm(n=imp.samp.num, mean=post.mean, sigma=post.cov)
imp.samp$log.imp <- dmvnorm(x=imp.samp$samples, mean=post.mean, sigma=post.cov, log=TRUE)
# compute posterior log-likelihood of importance samples
# is the log-likelihood function vectorized?
imp.samp$log.p <- log_post_ppgpd(parameters=imp.samp$samples,
                              parnames=parnames_all[[gpd.model]],
                              priors=priors[[gpd.model]],
                              data_calib=data_calib[[gpd.exp]])

print('done!')
  
# get posterior samples
print('sampling from posterior distribution...')

post.samp <- setNames(vector("list",length(samp.names)),samp.names)
samp.idx <- sample(x=nrow(post.samples), size=post.samp.num, replace=TRUE)
post.samp$samples <- post.samples[post.samp.ind,]
# get posterior log-likelihood of sampled posterior values
post.samp$log.p <- post.ll[post.samp.ind]
# get importance log-likelhood of posterior samples
post.samp$log.imp <- dmvnorm(x=post.samp$samples, mean=post.mean, sigma=post.cov, log=TRUE)

print('done!')

# compute the bridge sampling estimate using the iterative procedure from Meng and Wong (1996)
# starting value for the recursion is the reciprocal importance sampling estimate 
# from Gelfand and Dey (1994)
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
  # get number of samples
  post.num <- length(post$log.p)
  imp.num <- length(imp$log.p)

  # normalize posterior likelihoods based on previous normalizing constant estimate
  post$log.p.norm <- post$log.p - log.norm.const
  imp$log.p.norm <- imp$log.p - log.norm.const
  
  # compute updated estimate numerator and denominator
  imp.mean <- with(imp, mean(exp(log.p.norm)/(imp.num*exp(log.imp)+post.num*exp(log.p.norm))))
  post.mean <- with(post, mean(exp(log.imp)/(imp.num*exp(log.imp)+post.num*exp(log.p.norm))))
  
  # return updated estimate
  log.norm.const + log(imp.mean) - log(post.mean)
}

# set tolerance for halting of iteration
TOL <- 1e-3
# initialize storage for estimates
ml <- mat.or.vec(nr=1,nc=1)
# initialize with starting value
ml[1] <- recip.imp.samp(post.samp$log.lik,post.samp$log.imp)
# compute second iteration
ml[2] <- bridge.samp.iter(ml[1], post.samp[c('log.p','log.imp')], imp.samp[c('log.p','log.imp')])
# iterate until within tolerance
t <- 2
while (abs(ml[t] - ml[t-1]) >= TOL) {
  ml[t+1] <- bridge.samp.iter(ml[t], post.samp[c('log.p', 'log.imp')], imp.samp[c('log.p', 'log.imp')])
  t <- t+1
}

print('done!')

print('computing relative standard error of estimate')

# compute the relative standard error of the bridge sampling estimator
# since the posterior samples are not necessarily iid, we use the error formula from
# Fruhwirth-Schnatter (2004) which accounts for any autocorrelation
bridge.samp.rel.err <- function(log.norm.const,
                                post,
                                imp) {
  # get number of samples
  post.num <- length(post$log.p)
  imp.num <- length(imp$log.p)
  imp.factor <- imp.num/(imp.num+p.num)
  post.factor <- post.num/(imp.num+post.num)
  
  # normalize posterior likelihoods using normalizing constant estimate
  post$log.p.norm <- post$log.p - log.norm.const
  imp$log.p.norm <- imp$log.p - log.norm.const
  
  # compute bridging function process estimates for both sequences
  post$bridge.int <- with(post,exp(log.imp)/(imp.factor*exp(log.imp)+post.factor*(exp(log.post.norm))))
  imp$bridge.int <- with(imp,exp(log.p.norm)/(imp.factor*exp(log.imp)+post.factor*(exp(log.post.norm))))
  
  # compute frequency-0 autocorrelation of posterior process
  rho.0 <- spectrum0.ar(post$bridge.int)$spec
  
  # return squared relative error estimate
  (var(imp$bridge.int)/mean(imp$bridge.int)^2)/imp.num + (rho.f*var(post$bridge.int)/mean(post$bridge.int)^2)/post.num
}

re.sq <- bridge.samp.rel.err(ml[length(ml)], post.samp[c('log.p','log.imp')], imp.samp[c('log.p','log.imp')])

# save result of run
# if save directory doesn't exist, create it
ifelse(!dir.exists(path.save), dir.create(path.save), FALSE)
setwd(path.save)
# set output (saved as .RData; to be collected into a single output file later) file name
filename.out <- paste('ml_',station,'_',gpd.model,'_',data.length,'.RData',sep='')
save(list=c('post.samp','imp.samp', 'ml', 're.sq'), file=filename.out)
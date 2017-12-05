##==============================================================================
##  timeblock_motivating_experiment_fig1.R
##
##  Estimate storm surge PP/GPD parameters for one of the three sites using
##  30-year blocks and stationary model. To assess how sensitive the estimated
##  return levels are to data availability/length.
##
##  This version for submitting on HPC using:
##    qsub data_sensitivity_run.pbs
##
##  You will need to fix the 'setwd(...)' call below to match your local directory
##  structure, and perhaps need to 'install.packages(...)' for some of the libraries
##  that you may not have (all of the ones you should need are in the initial
##  section of code just below this note).
##
## Questions? Tony Wong (twong@psu.edu)
##==============================================================================

#setwd('~/codes/EVT/R')
setwd('/home/scrim/axw322/codes/EVT/R')

rm(list=ls())

output.dir <- '../output/'
dat.dir <- '../data/'

site <- 'Delfzijl'

nnode_mcmc000 <- 10
niter_mcmc000 <- 1e5
gamma_mcmc000 <- 0.5

library(adaptMCMC)
library(extRemes)
library(date)
library(Hmisc)
library(zoo)
library(ncdf4)

filename.priors   <- 'surge_priors_normalgamma_ppgpd_26Jul2017.rds'  # file holding the 'priors' object
filename.mles <- 'surge_MLEs_ppgpd_26Jul2017.rds'  # file holding the 'deoptim.all' object with the MLEs (for initial parameters)

# Name the saved progress RData workspace image file
today <- Sys.Date(); today=format(today,format="%d%b%Y")
filename.saveprogress <- paste(output.dir,'sensitivity_returnlevels_mcmc_',site,'_',today,'.RData', sep='')

# Name the calibrated parameters output file
today <- Sys.Date(); today=format(today,format="%d%b%Y")
filename.returnlevels <- paste(output.dir,'sensitivity_returnlevels_',site,'_',today,'.rds',sep='')

##==============================================================================

if(site=='Delfzijl') {data.site <- readRDS(paste(dat.dir,'tidegauge_processed_delfzijl_26Jul2017.rds',sep=''))
} else if(site=='Balboa') {data.site <- readRDS(paste(dat.dir,'tidegauge_processed_balboa_26Jul2017.rds',sep=''))
} else if(site=='Norfolk') {data.site <- readRDS(paste(dat.dir,'tidegauge_processed_norfolk_26Jul2017.rds',sep=''))
} else {print('ERROR - unrecognized site name')}

##==============================================================================
## Decompose into sets of 30-year blocks.
## Moving window, by 'block.offset' years each time.
##==============================================================================

block.size   <- 30 # how many years in each block?
block.offset <- 10 # how many years are the blocks shifted relative to neighboring blocks?
ind.block.right.endpt <- seq(from=length(data.site$gpd$year), to=block.size, by=-block.offset)
nblocks <- length(ind.block.right.endpt)
ind.block.left.endpt <- ind.block.right.endpt-block.size+1

# create some other data objects for calibration for each experiment block
block.names <- NULL
for (bb in 1:nblocks) {block.names <- c(block.names, paste('block',bb,sep=''))}
data.blocks <- vector('list', nblocks); names(data.blocks) <- block.names
for (bb in 1:nblocks) {
  # initialize each block with all of the data, then trim down to just that block
  data.blocks[[bb]] <- data.site$gpd
  ind.block <- ind.block.left.endpt[bb]:ind.block.right.endpt[bb]
  data.blocks[[bb]]$counts <- data.site$gpd$counts[ind.block]
  data.blocks[[bb]]$time_length <- data.site$gpd$time_length[ind.block]
  data.blocks[[bb]]$excesses <- data.site$gpd$excesses[ind.block]
  data.blocks[[bb]]$year <- data.site$gpd$year[ind.block]
  data.blocks[[bb]]$counts_all <- NULL
  data.blocks[[bb]]$time_length_all <- NULL
  data.blocks[[bb]]$excesses_all <- NULL
}

##==============================================================================
## Preliminary MCMC steps
##==============================================================================

source('parameter_setup_dayPOT.R')
parnames <- parnames_all$gpd3

priors <- readRDS(paste(output.dir,filename.priors,sep=''))

mle.fits <- readRDS(paste(output.dir,filename.mles,sep=''))
if(site=='Delfzijl') {initial.parameters <- mle.fits$gpd3['delfzijl',]
} else if(site=='Balboa') {initial.parameters <- mle.fits$gpd3[10,]
} else if(site=='Norfolk') {initial.parameters <- mle.fits$gpd3['norfolk',]}
step.mcmc <- apply(X=mle.fits$gpd3, MARGIN=2, FUN=sd) # mcmc step estimtae (will be adapted within the algorithm though)

## Set up prior distributions
parnames <- c('lambda','sigma','xi')

source('likelihood_ppgpd.R')

nnode_mcmc <- nnode_mcmc000
niter_mcmc <- niter_mcmc000
gamma_mcmc <- gamma_mcmc000
startadapt_mcmc <- max(500,round(0.05*niter_mcmc))
stopadapt_mcmc <- round(niter_mcmc*1.0)
accept_mcmc_few <- 0.44         # optimal for only one parameter
accept_mcmc_many <- 0.234       # optimal for many parameters
accept_mcmc <- accept_mcmc_many + (accept_mcmc_few - accept_mcmc_many)/length(parnames)

amcmc_out <- vector('list', nblocks); names(amcmc_out) <- block.names

##==============================================================================
## Actually run the MCMC - fit PP/GPD parameters (stationary model)
##==============================================================================

print(paste('Starting calibrations for ',site,': ',nnode_mcmc,' cores x ',niter_mcmc,' iterations)...', sep=''))
for (bb in 1:nblocks) {
  print(paste('-- block ',bb,' / ',nblocks, sep=''))
  tbeg <- proc.time()
  if(nnode_mcmc==1) {
    amcmc_out[[bb]] <- MCMC(log_post_ppgpd, niter_mcmc, initial.parameters,
                            adapt=TRUE, acc.rate=accept_mcmc, scale=step.mcmc,
                            gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                            parnames=parnames, data_calib=data.blocks[[bb]],
                            priors=priors, model='gpd3')
  } else if (nnode_mcmc > 1) {
    amcmc_out[[bb]] <- MCMC.parallel(log_post_ppgpd, niter_mcmc, initial.parameters,
                            n.chain=nnode_mcmc, n.cpu=nnode_mcmc, packages='extRemes',
                            scale=step.mcmc, adapt=TRUE, acc.rate=accept_mcmc,
                            gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                            parnames=parnames, data_calib=data.blocks[[bb]],
                            priors=priors, model='gpd3')
  } else {print('nnode_mcmc must be positive integer')}
  tend <- proc.time()
  print(paste('... done. Took ',round(as.numeric(tend-tbeg)[3]/60,2),' minutes', sep=''))
  print(paste('... Saving MCMC workspace results as .RData file (',filename.saveprogress,') to read and use later...',sep=''))
  save.image(file=filename.saveprogress)
  print('...done.')
}

##==============================================================================
## Convergence diagnostics and burn-in removal
##==============================================================================

## Gelman and Rubin diagnostics - determine and chop off for burn-in
niter.test <- seq(from=round(0.1*niter_mcmc), to=niter_mcmc, by=round(0.05*niter_mcmc))
gr.test <- mat.or.vec(length(niter.test), nblocks)
colnames(gr.test) <- block.names

for (bb in 1:nblocks) {
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
        eval(parse(text=paste('mcmc',m,' <- as.mcmc(amcmc_out[[bb]][[m]]$samples[1:niter.test[i],])', sep='')))
      }
      eval(parse(text=paste('mcmc_chain_list = mcmc.list(list(', string.mcmc.list , '))', sep='')))

      gr.test[i,bb] <- as.numeric(gelman.diag(mcmc_chain_list)[2])
    }
  } else {print('error - nnode_mcmc < 1 makes no sense')}
}

# 'ifirst' is the first spot where the GR stat gets to and stays below gr.max
# for all of the models. save a separate ifirst for each site
gr.max <- 1.1
if(nnode_mcmc==1) {
  ifirst <- round(0.5*niter_mcmc)
} else {
  lgr <- rep(NA, length(niter.test))
  for (i in 1:length(niter.test)) {lgr[i] <- all(gr.test[i,] < gr.max)}
  for (i in seq(from=length(niter.test), to=1, by=-1)) {
    if( all(lgr[i:length(lgr)]) ) {ifirst <- niter.test[i]}
  }
}

chains_burned <- vector('list', nblocks); names(chains_burned) <- block.names
for (bb in 1:nblocks) {
  if(nnode_mcmc > 1) {
    chains_burned[[bb]] <- vector('list', nnode_mcmc)
    for (m in 1:nnode_mcmc) {
      chains_burned[[bb]][[m]] <- amcmc_out[[bb]][[m]]$samples[(ifirst+1):niter_mcmc,]
    }
  } else {
    chains_burned[[bb]] <- amcmc_out[[bb]]$samples[(ifirst+1):niter_mcmc,]
  }
}

##==============================================================================
## Combine the parallel chains and thin to a manageable size?
##==============================================================================

# If no thinning, then this initialization will remain
chains_burned_thinned <- chains_burned

n.sample <- 100000

for (bb in 1:nblocks) {
  if(nnode_mcmc == 1) {
    ind.sample <- sample(x=1:nrow(chains_burned[[bb]]), size=n.sample, replace=FALSE)
    chains_burned_thinned[[bb]] <- chains_burned[[bb]][ind.sample,]
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
      ind.sample <- sample(x=1:nrow(chains_burned[[bb]][[m]]), size=n.sample.sub[m], replace=FALSE)
      chains_burned_thinned[[bb]][[m]] <- chains_burned[[bb]][[m]][ind.sample,]
    }
  }
}

# Combine all of the chains from 'ifirst' to 'niter_mcmc' into a potpourri of
# [alleged] samples from the posterior. Only saving the transition covariance
# matrix for one of the chains (if in parallel).
parameters.posterior <- vector('list', nblocks); names(parameters.posterior) <- block.names
for (bb in 1:nblocks) {
  if(nnode_mcmc==1) {
    parameters.posterior[[bb]] <- chains_burned_thinned[[bb]]
  } else {
    parameters.posterior[[bb]] <- chains_burned_thinned[[bb]][[1]]
    for (m in 2:nnode_mcmc) {
      parameters.posterior[[bb]] <- rbind(parameters.posterior[[bb]], chains_burned_thinned[[bb]][[m]])
    }
  }
}

##==============================================================================
## Calculate 100-year return level for each of the blocks, for each site.
##==============================================================================

returnperiod.of.interest <- 100 # in years

returnlevel <- vector('list', nblocks); names(returnlevel) <- block.names
for (bb in 1:nblocks) {
  returnlevel[[bb]] <- sapply(1:nrow(parameters.posterior[[bb]]), function(i) {
                              rlevd(returnperiod.of.interest,
                                    rate=parameters.posterior[[bb]][i,1],
                                    scale=parameters.posterior[[bb]][i,2],
                                    shape=parameters.posterior[[bb]][i,3],
                                    threshold=data.blocks[[bb]]$threshold,
                                    type='GP', npy=365.25)})
}

## Save progress to revisit later
print(paste('... Saving MCMC workspace results as .RData file (',filename.saveprogress,') to read and use later...',sep=''))
save.image(file=filename.saveprogress)
print('...done.')

##==============================================================================
## Get kernel density estimates for each of the distributions and plot
##==============================================================================

returnlevel.kde <- vector('list', nblocks); names(returnlevel.kde) <- block.names
for (bb in 1:nblocks) {
  returnlevel.kde[[bb]] <- density(returnlevel[[bb]], from=0, to=15000, n=512)
  returnlevel.kde[[bb]]$x <- returnlevel.kde[[bb]]$x/1000 # convert to m from mm
}

# note - can quantify the uncertainty in the PP-year return level by checking out
# the variance in the maximum likelihood posterior estimates?

# block years
block.years <- cbind(data.site$gpd$year[ind.block.left.endpt], data.site$gpd$year[ind.block.right.endpt])
names.block.years <- rep(NA, max(nblocks))
for (bb in 1:length(names.block.years)) {names.block.years[bb] <- paste(block.years[bb,1],block.years[bb,2],sep='-')}

##==============================================================================
## Save resulting return levels
##==============================================================================

saveRDS(returnlevel, file=filename.returnlevels)
save.image(file=filename.saveprogress)

##
##==============================================================================
## End
##==============================================================================
##

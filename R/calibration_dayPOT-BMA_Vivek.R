#===============================================================================
#
#===============================================================================

setwd('~/codes/EVT/R')

# each of these RData files has the raw MCMC output object 'amcmc_out' on it, for
# the site Balboa, Panama; Norfolk, Virgina, USA; or Delfzijl, the Netherlands.

mcmc.balboa <- '../output/everything_mcmc_ppgpd-experiments_balboa_normalgamma_28Jul2017.RData'
mcmc.norfolk <- '../output/everything_mcmc_ppgpd-experiments_norfolk_normalgamma_27Jul2017.RData'
mcmc.delfzijl <- '../output/everything_mcmc_ppgpd-experiments_delfzijl_normalgamma_29Jul2017.RData'

# each amcmc_out object has the following elements:
#   gpd30, gpd50, gpd70, ...
# where the number refers to the number of years of data used for the
# calibration (most recent years).
# We would like to determine the BMA weights for each of these data length experiments.
# The 'holy grail' plot will be showing that as more data become available, the
# weights shift toward non-stationary models.
# balboa goes up to gpd107, norfolk to gpd89, and delfzijl to gpd137.
#
# each of the amcmc_out[[data.length]] objects (e.g., amcmc_out$gpd30) has the
# following elements:
#   gpd3, gpd4, gpd5, gpd6
# where the number refers to the total number of parameters. gpd3 is a stationary
# model (poisson rate, and generalized pareto scale and shape parameters); gpd4
# has a non-stationary poisson rate parameter, gpd5 has non-stationary rate and
# scale parameters, and gpd6 has all non-stationary parameters.
# These are the 4 models we would like to determine the BMA weights for. So for
# a given site ('site'), and for a given data length experiment ('data.exp'),
# use the output from each of amcmc_out[[data.exp]]$gpd3, $gpd4, $gpd5 and $gpd6
# to determine the weights.
#
# each of amcmc_out[[data.exp]][[model]] has the results of 10 parallel
# MCMC chains of 500,000 iterations each ('niter_mcmc'). So the raw parameter
# output for one of the chains for a site for the experiment using 50 years of
# data, calibrating the stationary GPD model would be:
#   amcmc_out$gpd50$gpd3[[1]]$samples
#   amcmc_out$gpd50$gpd3[[2]]$samples
#   ...
#   amcmc_out$gpd50$gpd3[[10]]$samples
#
# If you want the log-posterior values from the chains, you can get them at:
#   amcmc_out$gpd50$gpd3[[1]]$log.p
# for example.
#
# We also have the transition covariance matrices, because we use an adaptive
# sampler, because we are not neanderthals:
#   amcmc_out$gpd50$gpd3[[1]]$cov.jump
#
# Note that I have not removed for burn-in, nor have I thinned and/or shuffled
# the results from all these different chains together.
# There is a quantity called 'ifirst' on the RData files which is the result of
# calculating the Gelman and Rubin diagnostics for each of the data length
# experiments, and is the first iteration (after some initial period) where the
# scale reduction factor for all of the GPD models for that data experiment, for
# that site, are less than and stay less than 1.1. So you could use 'ifirst' to
# pick burn-in (just take ifirst:niter_mcmc from each of the chains).




#===============================================================================
# end
#===============================================================================

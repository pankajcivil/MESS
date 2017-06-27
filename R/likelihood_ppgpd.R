#===============================================================================
# likelihood function(s), priors, posterior
# for generalized pareto/poisson process distribution
#
# The idea is that a GPD governs how excesses above a given threshold are
# distributed, but this is conditioned on the probability that we observe an
# excess. These exceedance observations are governed by a Poisson process, where
# the rate parameter gives (1/) the expected number of exceedances per year (or
# other time unit desired; here I use a year).
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================


#
#===============================================================================
# pdf and cdf for pp-gpd (that can handle arrays)
# note that this assumes you aren't a jerk and send in different length
# arrays. so... don't do that.
# This assumes that you only send in
#===============================================================================


# TODO
# TODO - modifying below here in progress!
# TODO

# log-likelihood given by Coles (2001), eq. 4.10:
# llike <- -length(data)*log(sigma) - (1+1/xi)*sum(log(1+xi*data/sigma))

# use eq. 11 of Martins and Stedinger (2001), for each year independently:
# llike <- log([likelihood of seeing exactly m exceedances in time interval dt with rate lambda[t]])
#           + log([joint GPD density for the m exceedances in time interval dt])
# ^^^ do this for all of the intervals

# need to send in for each interval:
# 1. length of intervals
# 2. number of exceedances each interval (assumed to be [level-threshold])
# 3. all parameters (rate, scale, shape) as time series (annual values/each block)
# 4. values of the exceedance levels (level-threshold)

# TODO
# TODO

# NOTE: don't use this one for now; calculate in likelhiood function using
# devd(type='GP', log=TRUE) and
# dpois(x=[counts this interval], lambda=[actual lambda]*[interval length], log=TRUE)

#===============================================================================


#
#===============================================================================
# project PP-GPD parameters
#===============================================================================
#
project_ppgpd <- function(parameters,
                          parnames,
                          auxiliary
){
  parameters_project <- mat.or.vec(length(auxiliary), 3)
  colnames(parameters_project) <- c('lambda','sigma','xi')
  n.param <- length(parnames)
  if(n.param==3) {
    # fit a standard stationary PP-GPD
    lambda <- rep(parameters[match('lambda',parnames)], length(auxiliary))
    sigma <- rep(parameters[match('sigma',parnames)], length(auxiliary))
    xi <- rep(parameters[match('xi',parnames)], length(auxiliary))
  } else {print('ERROR - invalid number of parameters for PP-GPD')}


# TODO
# TODO - add support for nonstationary models
# TODO


  parameters_project[,'lambda'] <- lambda
  parameters_project[,'sigma'] <- sigma
  parameters_project[,'xi'] <- xi

  return(parameters_project)
}
#===============================================================================


#
#===============================================================================
# log(prior) for PP-GPD model
#===============================================================================
#
log_prior_ppgpd <- function(parameters,
                            parnames,
                            priors,
                            model,
                            auxiliary=NULL
){
  lpri <- 0

  for (par in parnames) {
    parameter.value <- as.numeric(parameters[match(par,parnames)])
    if(priors[[model]][[par]]$type=='normal') {
      lpri <- lpri + dnorm(x=parameter.value, mean=priors[[model]][[par]]$mean, sd=priors[[model]][[par]]$sd, log=TRUE)
    } else if(priors[[model]][[par]]$type=='gamma') {
      lpri <- lpri + dgamma(x=parameter.value, shape=priors[[model]][[par]]$shape, rate=priors[[model]][[par]]$rate, log=TRUE)
    } else if(priors[[model]][[par]]$type=='uniform') {
      lpri <- lpri + dunif(x=parameter.value, min=priors[[model]][[par]]$lower, max=priors[[model]][[par]]$upper, log=TRUE)
    }
  }

  return(lpri)
}
#===============================================================================


#
#===============================================================================
# -log(likelihood) for PP-GPD model
#===============================================================================
#
neg_log_like_ppgpd <- function(parameters,
                               parnames,
                               data_calib,
                               auxiliary=NULL
){
  nll <- -1 * log_like_ppgpd(parameters, parnames, data_calib, auxiliary)
  return(nll)
}
#===============================================================================


#
#===============================================================================
# log(likelihood) for PP-GPD model
#===============================================================================
#
log_like_ppgpd <- function(parameters,
                           parnames,
                           data_calib,
                           auxiliary=NULL
){
  llik <- 0
  #print(parameters)
  n.param <- length(parnames)
  if(n.param==3) {
    # fit a standard stationary PP-GPD
    lambda <- parameters[match('lambda',parnames)]
    sigma <- parameters[match('sigma',parnames)]
    xi <- parameters[match('xi',parnames)]
  } else {print('ERROR - invalid number of parameters for PP-GPD')}


# TODO
# TODO - add support for nonstationary models
# TODO


# devd(type='GP', log=TRUE) and
# dpois(x=[counts this interval], lambda=[actual lambda]*[interval length], log=TRUE)


  # this way is working, but slow. speed up with an "apply" ?
  nbins <- length(data_calib$gpd$counts)
  llik.bin <- rep(NA, nbins)
  for (b in 1:nbins) {
    if(data_calib$gpd$counts[b] > 0) {
      llik.bin[b] <- dpois(x=data_calib$gpd$counts[b], lambda=(lambda*data_calib$gpd$time_length[b]), log=TRUE) +
                     sum(devd(data_calib$gpd$excesses[[b]]-data_calib$gpd$threshold, threshold=0, scale=sigma, shape=xi, log=TRUE, type='GP'))
    } else {
      llik.bin[b] <- dpois(x=data_calib$gpd$counts[b], lambda=(lambda*data_calib$gpd$time_length[b]), log=TRUE)
    }
  }
  llik <- sum(llik.bin)



# this way works, but cannot take nonstationarity in the poisson process
#if(data_calib$gpd$counts_all > 0) {
#  llik <- dpois(x=data_calib$gpd$counts_all, lambda=(lambda*data_calib$gpd$time_length_all), log=TRUE) +
#          sum(devd(data_calib$gpd$excesses_all-data_calib$gpd$threshold, threshold=0, scale=sigma, shape=xi, log=TRUE, type='GP'))
#} else {
#  llik <- dpois(x=data_calib$gpd$counts_all, lambda=(lambda*data_calib$gpd$time_length_all), log=TRUE)
#}


  return(llik)
}
#===============================================================================


#
#===============================================================================
# log(post) for PP-GPD model
#===============================================================================
#
log_post_ppgpd <- function(parameters,
                           parnames,
                           data_calib,
                           priors,
                           model,
                           auxiliary
){
  lpost <- 0
  llik <- 0
  lpri <- 0

  # calculate prior
  lpri <- log_prior_ppgpd(parameters=parameters,
                          parnames=parnames,
                          priors=priors,
                          model=model,
                          auxiliary=auxiliary)

  if(is.finite(lpri)){
    # calculate likelihood (only if parameters pass the prior test)
    llik <- log_like_ppgpd(parameters=parameters,
                           parnames=parnames,
                           data_calib=data_calib,
                           auxiliary=auxiliary)
  }

  lpost <- lpri + llik
  return(lpost)
}
#===============================================================================


#===============================================================================
# End
#===============================================================================

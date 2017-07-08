#===============================================================================
# likelihood function(s), priors, posterior
# for GEV distribution
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================


#
#===============================================================================
# pdf and cdf for gev (that can handle arrays)
# note that this assumes you aren't a jerk and send in different length
# arrays. so... don't do that.
#===============================================================================

gev_pdf <- function(x, loc, scale, shape, log=FALSE){
  p <- rep(NA,length(x))
  if(all(shape==0)) {
    if(log) {
      p <- (loc-x)/scale - exp((loc-x)/scale) - log(scale)
    } else {
      p <- (1/scale) * exp((loc-x)/scale) * exp( -( exp((loc-x)/scale) ) )
    }
  } else {
    if(log) {
      p <- (-1 - 1/shape)*log(1+shape*((x-loc)/scale)) - ((1+shape*((x-loc)/scale))^(-1/shape)) - log(scale)
      if(any(shape < 0)) {
        i1 <- which(shape < 0); i2 <- which(x > (loc-scale/shape)); i3 <- intersect(i1,i2)
        p[i3] <- -Inf
      } else if(any(shape > 0)) {
        i1 <- which(shape > 0); i2 <- which(x < (loc-scale/shape)); i3 <- intersect(i1,i2)
        p[i3] <- -Inf
      }
    } else {
      p <- (1/scale) * ( (1+shape*((x-loc)/scale))^(-1 - 1/shape) ) * exp( -(1+shape*((x-loc)/scale))^(-1/shape) )
      if(any(shape < 0)) {
        i1 <- which(shape < 0); i2 <- which(x > (loc-scale/shape)); i3 <- intersect(i1,i2)
        p[i3] <- 0
      } else if(any(shape > 0)) {
        i1 <- which(shape > 0); i2 <- which(x < (loc-scale/shape)); i3 <- intersect(i1,i2)
        p[i3] <- 0
      }
    }
  }
  return(p)
}

gev_cdf <- function(q, loc, scale, shape){
  p <- rep(NA,length(q))
  if(all(shape==0)) {
    p <- exp( -( exp(-(q-loc)/scale) ) )
  } else {
    p <- exp( -(1+shape*((q-loc)/scale))^(-1/shape) )
    if(any(shape < 0)) {
      i1 <- which(shape < 0); i2 <- which(q > (loc-scale/shape)); i3 <- intersect(i1,i2)
      # if shape < 0, i3 is all the places where q > theoretical upper bound
      # -> there ought to be 100% probability mass below here
      p[i3] <- 1
    } else if(any(shape > 0)) {
      i1 <- which(shape > 0); i2 <- which(q < (loc-scale/shape)); i3 <- intersect(i1,i2)
      # if shape > 0, i3 is all the places where q < theoretical upper bound
      # -> there ought to be 100% probability mass above here    }
      p[i3] <- 0
    }
  }
  return(p)
}
#===============================================================================


#
#===============================================================================
# project GEV parameters
#===============================================================================
#
project_gev <- function(parameters,
                        parnames,
                        auxiliary
){
  parameters_project <- mat.or.vec(length(auxiliary), 3)
  colnames(parameters_project) <- c('mu','sigma','xi')
  n.param <- length(parnames)
  if(n.param==3) {
    # fit a standard stationary GEV
    mu <- rep(parameters[match('mu',parnames)], length(auxiliary))
    sigma <- rep(parameters[match('sigma',parnames)], length(auxiliary))
    xi <- rep(parameters[match('xi',parnames)], length(auxiliary))
  } else if(n.param==4) {
    # location parameter nonstationary
    mu0 <- parameters[match('mu0',parnames)]
    mu1 <- parameters[match('mu1',parnames)]
    sigma <- rep(parameters[match('sigma',parnames)], length(auxiliary))
    xi <- rep(parameters[match('xi',parnames)], length(auxiliary))
    mu <- mu0 + mu1*auxiliary
  } else if(n.param==5) {
    # location and scale parameters nonstationary
    mu0 <- parameters[match('mu0',parnames)]
    mu1 <- parameters[match('mu1',parnames)]
    sigma0 <- parameters[match('sigma0',parnames)]
    sigma1 <- parameters[match('sigma1',parnames)]
    xi <- rep(parameters[match('xi',parnames)], length(auxiliary))
    mu <- mu0 + mu1*auxiliary
    sigma <- exp(sigma0 + sigma1*auxiliary)
  } else if(n.param==6) {
    # location, scale and shape all nonstationary
    mu0 <- parameters[match('mu0',parnames)]
    mu1 <- parameters[match('mu1',parnames)]
    sigma0 <- parameters[match('sigma0',parnames)]
    sigma1 <- parameters[match('sigma1',parnames)]
    xi0 <- parameters[match('xi0',parnames)]
    xi1 <- parameters[match('xi1',parnames)]
    mu <- mu0 + mu1*auxiliary
    sigma <- exp(sigma0 + sigma1*auxiliary)
    xi <- xi0 + xi1*auxiliary
  } else {print('ERROR - invalid number of parameters for GEV')}

  parameters_project[,'mu'] <- mu
  parameters_project[,'sigma'] <- sigma
  parameters_project[,'xi'] <- xi

  return(parameters_project)
}
#===============================================================================


#
#===============================================================================
# log(prior) for gev model
#===============================================================================
#
log_prior_gev <- function(parameters,
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
# -log(likelihood) for gev model
#===============================================================================
#
neg_log_like_gev <- function(parameters,
                             parnames,
                             data_calib,
                             auxiliary=NULL
){
  nll <- -1 * log_like_gev(parameters, parnames, data_calib, auxiliary)
  return(nll)
}
#===============================================================================


#
#===============================================================================
# log(likelihood) for gev model
#===============================================================================
#
log_like_gev <- function(parameters,
                         parnames,
                         data_calib,
                         auxiliary=NULL
){
  llik <- 0
  #print(parameters)
  n.param <- length(parnames)
  if(n.param==3) {
    # fit a standard stationary GEV
    mu <- parameters[match('mu',parnames)]
    sigma <- parameters[match('sigma',parnames)]
    xi <- parameters[match('xi',parnames)]
  } else if(n.param==4) {
    # location parameter nonstationary
    mu0 <- parameters[match('mu0',parnames)]
    mu1 <- parameters[match('mu1',parnames)]
    sigma <- parameters[match('sigma',parnames)]
    xi <- parameters[match('xi',parnames)]
    mu <- mu0 + mu1*auxiliary
  } else if(n.param==5) {
    # location and scale parameters nonstationary
    mu0 <- parameters[match('mu0',parnames)]
    mu1 <- parameters[match('mu1',parnames)]
    sigma0 <- parameters[match('sigma0',parnames)]
    sigma1 <- parameters[match('sigma1',parnames)]
    xi <- parameters[match('xi',parnames)]
    mu <- mu0 + mu1*auxiliary
    sigma <- exp(sigma0 + sigma1*auxiliary)
  } else if(n.param==6) {
    # location, scale and shape all nonstationary
    mu0 <- parameters[match('mu0',parnames)]
    mu1 <- parameters[match('mu1',parnames)]
    sigma0 <- parameters[match('sigma0',parnames)]
    sigma1 <- parameters[match('sigma1',parnames)]
    xi0 <- parameters[match('xi0',parnames)]
    xi1 <- parameters[match('xi1',parnames)]
    mu <- mu0 + mu1*auxiliary
    sigma <- exp(sigma0 + sigma1*auxiliary)
    xi <- xi0 + xi1*auxiliary
  } else {print('ERROR - invalid number of parameters for GEV')}
#  llik <- sum(gev_pdf(data_calib, loc=mu, scale=sigma, shape=xi, log=TRUE))
  llik <- sum(devd(data_calib, loc=mu, scale=sigma, shape=xi, log=TRUE, type='GEV'))
  return(llik)
}
#===============================================================================


#
#===============================================================================
# log(post) for gev model
#===============================================================================
#
log_post_gev <- function(parameters,
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
  lpri <- log_prior_gev(parameters=parameters,
                        parnames=parnames,
                        priors=priors,
                        model=model,
                        auxiliary=auxiliary)

  if(is.finite(lpri)){
    # calculate likelihood (only if parameters pass the prior test)
    llik <- log_like_gev(parameters=parameters,
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

#===============================================================================
# likelihood function(s), priors, posterior
# for Naveau, model (i)
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================


#
#===============================================================================
# pdf and cdf functions for Naveau (i) model (or Papastathopoulos and Tawn 2013
# model (iii))
#===============================================================================
#
naveau_pdf <- function(x, kappa, sigma, xi){
  prob <- NULL
  prob <- kappa*((1-(1+xi*x/sigma)^(-1/xi))^(kappa-1))*((1+xi*x/sigma)^(-(1+1/xi)))/sigma
  return(prob)
}

naveau_cdf <- function(x, kappa, sigma, xi){
  q <- NULL
  q <- (1-(1+xi*x/sigma)^(-1/xi))^kappa
  return(q)
}
#===============================================================================


#
#===============================================================================
# log(prior) for naveau (i) model
#===============================================================================
#
log_prior_naveau <- function(parameters,
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
    } else if(priors[[model]][[p]]$type=='gamma') {
      lpri <- lpri + dgamma(x=parameter.value, shape=priors[[model]][[par]]$shape, rate=priors[[model]][[par]]$rate, log=TRUE)
    }
  }

  return(lpri)
}
#===============================================================================


#
#===============================================================================
# -log(likelihood) for naveau (i) model
#===============================================================================
#
neg_log_like_naveau <- function(parameters,
                                parnames,
                                data_calib,
                                auxiliary=NULL
){
  nll <- -1 * log_like_naveau(parameters, parnames, data_calib, auxiliary)
  return(nll)
}
#===============================================================================


#
#===============================================================================
# log(likelihood) for naveau (i) model
#===============================================================================
#
log_like_naveau <- function(parameters,
                            parnames,
                            data_calib,
                            auxiliary=NULL
){
  llik <- 0
  n.param <- length(parnames)
  if(n.param==3) {
    # fit a standard stationary Naveau-(i)
    kappa <- parameters[match('kappa',parnames)]
    sigma <- parameters[match('sigma',parnames)]
    xi <- parameters[match('xi',parnames)]
  } else if(n.param==4) {
    # lower-tail parameter nonstationary
    kappa0 <- parameters[match('kappa0',parnames)]
    kappa1 <- parameters[match('kappa1',parnames)]
#    sigma <- parameters[match('sigma',parnames)]
#    xi <- parameters[match('xi',parnames)]
    sigma <- rep(parameters[match('sigma',parnames)], length(auxiliary))
    xi <- rep(parameters[match('xi',parnames)], length(auxiliary))
    kappa <- kappa0 + kappa1*auxiliary
  } else if(n.param==5) {
    # lower-tail and scale parameters nonstationary
    kappa0 <- parameters[match('kappa0',parnames)]
    kappa1 <- parameters[match('kappa1',parnames)]
    sigma0 <- parameters[match('sigma0',parnames)]
    sigma1 <- parameters[match('sigma1',parnames)]
    xi <- parameters[match('xi',parnames)]
    kappa <- kappa0 + kappa1*auxiliary
    sigma <- exp(sigma0 + sigma1*auxiliary)
  } else if(n.param==6) {
    # lower-tail, scale and shape all nonstationary
    kappa0 <- parameters[match('kappa0',parnames)]
    kappa1 <- parameters[match('kappa1',parnames)]
    sigma0 <- parameters[match('sigma0',parnames)]
    sigma1 <- parameters[match('sigma1',parnames)]
    xi0 <- parameters[match('xi0',parnames)]
    xi1 <- parameters[match('xi1',parnames)]
    kappa <- kappa0 + kappa1*auxiliary
    sigma <- exp(sigma0 + sigma1*auxiliary)
    xi <- xi0 + xi1*auxiliary
  } else {print('ERROR - invalid number of parameters for Naveau-(i)')}

  llik <- sum( log(naveau_pdf(x=data_calib, kappa=kappa, sigma=sigma, xi=xi)) )
  if(is.na(llik)) {llik <- -9e9}

  return(llik)
}
#===============================================================================


#
#===============================================================================
# log(post) for naveau (i) model
#===============================================================================
#
log_post_naveau <- function(parameters,
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
  lpri <- log_prior_naveau(parameters=parameters,
                           parnames=parnames,
                           priors=priors,
                           model=model,
                           auxiliary=auxiliary)

  if(is.finite(lpri)){
    # calculate likelihood (only run model if parameters are reasonable)
    llik <- log_like_naveau(parameters=parameters,
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

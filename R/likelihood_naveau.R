#===============================================================================
# likelihood function(s), priors, posterior
# for Naveau, model (i)
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================


#===============================================================================

log_prior_naveau <- function(parameters,
                             parnames,
                             priors,
                             auxiliary=NULL
){
  lpri <- 0

  # TODO

  # take in some list object "priors" for calculating this

  return(lpri)
}

#===============================================================================

neg_log_like_naveau <- function(parameters,
                                parnames,
                                data_calib,
                                auxiliary=NULL
){
  nll <- -1 * log_like(parameters, parnames, data_calib, auxiliary)
  return(nll)
}

#===============================================================================

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
    sigma <- parameters[match('sigma',parnames)]
    xi <- parameters[match('xi',parnames)]
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

# TODO CHECK!
  p1 <- devd(data_calib/sigma, scale=sigma, shape=xi, threshold=0, type="GP")
  p2 <- kappa*pevd(data_calib/sigma, scale=sigma, shape=xi, threshold=0, type="GP")^(kappa-1)
  p3 <- sigma^(-length(data_calib))
  llik <- sum(log(p1))+sum(log(p2))+log(p3)

  return(llik)
}

#===============================================================================

log_post_naveau <- function(parameters,
                            parnames,
                            data_calib,
                            priors,
                            auxiliary
){
  lpost <- 0
  llik <- 0
  lpri <- 0

  # TODO calculate prior

  if(is.finite(lpri)){
    # TODO calculate likelihood

  }

  lpost <- lpri + llik
  return(lpost)
}

#===============================================================================
# End
#===============================================================================

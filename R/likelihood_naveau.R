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
naveau_pdf <- function(x, kappa, sigma, xi, log=FALSE){
  prob <- NULL
  if(log) {prob <- log(kappa) - log(sigma) + (kappa-1)*log(1-(1+xi*x/sigma)^(-1/xi)) - (1+1/xi)*log(1+xi*x/sigma)
  } else  {prob <- kappa*((1-(1+xi*x/sigma)^(-1/xi))^(kappa-1))*((1+xi*x/sigma)^(-(1+1/xi)))/sigma}
  ibad <- which( ((x > (-sigma/xi)) & (xi < 0)) |  (x < 0) | (kappa < 0) )
  if(length(ibad)>0) {prob[ibad] <- 0; if(log) {prob[ibad] <- -Inf}}
  return(prob)
}

naveau_cdf <- function(x, kappa, sigma, xi){
  q <- NULL
  q <- (1-(1+xi*x/sigma)^(-1/xi))^kappa
  ibad.high <- which( ((x > (-sigma/xi)) & (xi < 0)) )
  if(length(ibad.high)>0) {q[ibad.high] <- 1}
  ibad.low  <- which( x < 0 )
  if(length(ibad.low )>0) {q[ibad.low]  <- 0}
  return(q)
}

naveau_invcdf <- function(q, kappa, sigma, xi){
  x <- NULL
  x <- sigma*((1-q^(1/kappa))^(-xi)-1)/xi
  return(x)
}
#===============================================================================


#
#===============================================================================
# project naveau parameters
#===============================================================================
#
project_naveau <- function(parameters,
                           parnames,
                           auxiliary
){
  parameters_project <- mat.or.vec(length(auxiliary), 3)
  colnames(parameters_project) <- c('kappa','sigma','xi')
  n.param <- length(parnames)
  if(n.param==3) {
    # fit a standard stationary Naveau-(i)
    if('kappa' %in% parnames)       {kappa <- rep(parameters[match('kappa',parnames)], length(auxiliary))}
    else if('lkappa' %in% parnames) {kappa <- rep(exp(parameters[match('lkappa',parnames)]), length(auxiliary))}
    sigma <- rep(parameters[match('sigma',parnames)], length(auxiliary))
    xi <- rep(parameters[match('xi',parnames)], length(auxiliary))
  } else if(n.param==4) {
    # lower-tail parameter nonstationary
    if('kappa0' %in% parnames)       {kappa0 <- rep(parameters[match('kappa0',parnames)], length(auxiliary))}
    else if('lkappa0' %in% parnames) {kappa0 <- rep(exp(parameters[match('lkappa0',parnames)]), length(auxiliary))}
    kappa1 <- rep(parameters[match('kappa1',parnames)], length(auxiliary))
    sigma <- rep(parameters[match('sigma',parnames)], length(auxiliary))
    xi <- rep(parameters[match('xi',parnames)], length(auxiliary))
    kappa <- kappa0 + kappa1*auxiliary
  } else if(n.param==5) {
    # lower-tail and scale parameters nonstationary
    if('kappa0' %in% parnames)       {kappa0 <- rep(parameters[match('kappa0',parnames)], length(auxiliary))}
    else if('lkappa0' %in% parnames) {kappa0 <- rep(exp(parameters[match('lkappa0',parnames)]), length(auxiliary))}
    kappa1 <- rep(parameters[match('kappa1',parnames)], length(auxiliary))
    sigma0 <- rep(parameters[match('sigma0',parnames)], length(auxiliary))
    sigma1 <- rep(parameters[match('sigma1',parnames)], length(auxiliary))
    xi <- rep(parameters[match('xi',parnames)], length(auxiliary))
    kappa <- kappa0 + kappa1*auxiliary
    sigma <- exp(sigma0 + sigma1*auxiliary)
  } else if(n.param==6) {
    # lower-tail, scale and shape all nonstationary
    if('kappa0' %in% parnames)       {kappa0 <- rep(parameters[match('kappa0',parnames)], length(auxiliary))}
    else if('lkappa0' %in% parnames) {kappa0 <- rep(exp(parameters[match('lkappa0',parnames)]), length(auxiliary))}
    kappa1 <- rep(parameters[match('kappa1',parnames)], length(auxiliary))
    sigma0 <- rep(parameters[match('sigma0',parnames)], length(auxiliary))
    sigma1 <- rep(parameters[match('sigma1',parnames)], length(auxiliary))
    xi0 <- rep(parameters[match('xi0',parnames)], length(auxiliary))
    xi1 <- rep(parameters[match('xi1',parnames)], length(auxiliary))
    kappa <- kappa0 + kappa1*auxiliary
    sigma <- exp(sigma0 + sigma1*auxiliary)
    xi <- xi0 + xi1*auxiliary
  } else {print('ERROR - invalid number of parameters for Naveau-(i)')}

  parameters_project[,'kappa'] <- kappa
  parameters_project[,'sigma'] <- sigma
  parameters_project[,'xi'] <- xi

  return(parameters_project)
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
                            auxiliary=NULL,
                            Tmax
){
  llik <- 0
  n.param <- length(parnames)
  if(n.param==3) {
    # fit a standard stationary Naveau-(i)
    if('kappa' %in% parnames)       {kappa <- rep(parameters[match('kappa',parnames)], length(data_calib))}
    else if('lkappa' %in% parnames) {kappa <- rep(exp(parameters[match('lkappa',parnames)]), length(data_calib))}
    sigma <- rep(parameters[match('sigma',parnames)], length(data_calib))
    xi <- rep(parameters[match('xi',parnames)], length(data_calib))
  } else if(n.param==4) {
    # lower-tail parameter nonstationary
    if('kappa0' %in% parnames)       {kappa0 <- rep(parameters[match('kappa0',parnames)], length(data_calib))}
    else if('lkappa0' %in% parnames) {kappa0 <- rep(exp(parameters[match('lkappa0',parnames)]), length(data_calib))}
    kappa1 <- rep(parameters[match('kappa1',parnames)], length(auxiliary))
    sigma <- rep(parameters[match('sigma',parnames)], length(auxiliary))
    xi <- rep(parameters[match('xi',parnames)], length(auxiliary))
    kappa <- kappa0 + kappa1*auxiliary
  } else if(n.param==5) {
    # lower-tail and scale parameters nonstationary
    if('kappa0' %in% parnames)       {kappa0 <- rep(parameters[match('kappa0',parnames)], length(data_calib))}
    else if('lkappa0' %in% parnames) {kappa0 <- rep(exp(parameters[match('lkappa0',parnames)]), length(data_calib))}
    kappa1 <- rep(parameters[match('kappa1',parnames)], length(auxiliary))
    sigma0 <- rep(parameters[match('sigma0',parnames)], length(auxiliary))
    sigma1 <- rep(parameters[match('sigma1',parnames)], length(auxiliary))
    xi <- rep(parameters[match('xi',parnames)], length(auxiliary))
    kappa <- kappa0 + kappa1*auxiliary
    sigma <- exp(sigma0 + sigma1*auxiliary)
  } else if(n.param==6) {
    # lower-tail, scale and shape all nonstationary
    if('kappa0' %in% parnames)       {kappa0 <- rep(parameters[match('kappa0',parnames)], length(data_calib))}
    else if('lkappa0' %in% parnames) {kappa0 <- rep(exp(parameters[match('lkappa0',parnames)]), length(data_calib))}
    kappa1 <- rep(parameters[match('kappa1',parnames)], length(auxiliary))
    sigma0 <- rep(parameters[match('sigma0',parnames)], length(auxiliary))
    sigma1 <- rep(parameters[match('sigma1',parnames)], length(auxiliary))
    xi0 <- rep(parameters[match('xi0',parnames)], length(auxiliary))
    xi1 <- rep(parameters[match('xi1',parnames)], length(auxiliary))
    kappa <- kappa0 + kappa1*auxiliary
    sigma <- exp(sigma0 + sigma1*auxiliary)
    xi <- xi0 + xi1*auxiliary
  } else {print('ERROR - invalid number of parameters for Naveau-(i)')}

  # check extra conditions that would otherwise 'break' the likelihood function
  if( any(kappa < 0) | any(sigma < 0) | any((1+xi*data_calib/sigma) < 0) | any((1-(1+xi*data_calib/sigma)^(-1/xi)) < 0) ) {llik <- -Inf}
  else {llik <- sum( naveau_pdf(x=data_calib, kappa=kappa, sigma=sigma, xi=xi, log=TRUE) )}

  # constraint on kappa0, kappa1 and Tmax
  if( exists('kappa0') & exists('kappa1') ) {
    if( kappa1[1] < (-kappa0[1]/Tmax) ) {llik <- -Inf}
  }

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
                            auxiliary,
                            Tmax
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
                            auxiliary=auxiliary,
                            Tmax=Tmax)
  }

  lpost <- lpri + llik
  return(lpost)
}
#===============================================================================


#===============================================================================
# End
#===============================================================================

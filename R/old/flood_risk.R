#===============================================================================
#
# questions? Tony Wong (twong@psu.edu)
#===============================================================================

#
#===============================================================================
# settings
#===============================================================================
#

filename.yearBM.parameters <- '../output/evt_models_calibratedParameters_gev-nav_21Jun2017.nc'
filename.dayPOT.parameters <- '../output/evt_models_calibratedParameters_ppgpd_28Jun2017.nc'
filename.sealevelrise <- '../../BRICK/output_model/BRICK-fastdyn_physical_gamma_01Jun2017.nc'
filename.datacalib <- 'datacalib_05Jul2017.rds'
filename.priors.dayPOT <- '../output/surge_priors_ppgpd_28Jun2017.rds'
filename.priors.yearBM <- '../output/surge_priors_gev_nav_19Jun2017.rds'

time.beg <- 2015           # inital year, "present"
time.end <- 2065           # final year (time horizon)
time.step <- 1             # time step in years
nmax.flood <- 183          # maximum number of flood events in a year considered
height.low <- 0            # lowest heightening considered [m]
height.high <- 10          # tallest heightening considered [m]
height.increment <- 0.05   # increments of heightening considered [m]
heightening <- seq(from=height.low, to=height.high, by=height.increment)
lat <- 51.9775
lon <- 4.12

# TODO
# TODO -- replace initial dike height with what ever it actually is, or calculate
# TODO -- from initial flood probability from Eijgenraam
# TODO -- (that would account for initial land below local sea level)
# TODO

height.initial <- 5        # initial dike ring height [m]

# use maximum posterior probability stationary GEV parameters to estimate H0
# consistent with Eijgenraam Table 1
#p0 <- 0.00137              # Eijgenraam et al (2012)
#gev.prelim <- amcmc_prelim$gev3$samples[which.max(amcmc_prelim$gev3$log.p),]
#height.initial <- qevd(p=1-p0, loc=gev.prelim[1], scale=gev.prelim[2], shape=gev.prelim[3])/1000

heightening <- seq(from=height.low, to=height.high, by=height.increment)

# useful
euro <- "\u20AC"

#
#===============================================================================
# additional parameters pertaining to flood risk
# (for dike ring 15, Eijgenraam et al 2012)
#===============================================================================
#

# TODO - sample, or hold constant?

subsidence.rate <- 0.002     # m/yr (Rietveld H. Land subsidence in the Netherlands. Pp. 455–465 in Proceedings of the 3rd International Symposium on Land Subsidence. Vol 151. 1986.)
discount.rate <- 0.04        # % (Eijgenraam et al 2012)
value.initial <- 11810.4     # million euro (Eijgenraam et al 2012)
econgrowth.rate <- 0.02      # /year (CPB Central Economic Plan 2017)
hgtdamage.rate <- 0.003764   # millino euro/cm heightening (Eijgenraam et al 2012)

#
#===============================================================================
# set up cost of dike heightening
#===============================================================================
#

cost.func <- vector('list',2); names(cost.func) <- c('exponential','quadratic')
cost.func$exponential <- vector('list',3); names(cost.func$exponential) <- c('c','b','lambda')
cost.func$exponential$c <- 125.6422       # dike ring 15 of Eijgenraam (Table 1)
cost.func$exponential$b <- 1.1268
cost.func$exponential$lambda <- 0.0098
cost.func$quadratic <- vector('list', 3); names(cost.func$quadratic) <- c('a0','b0','c0')
cost.func$quadratic$a0 <- 0.027           # dike ring 15 of Eijgenraam (Table 3)
cost.func$quadratic$b0 <- 3.779
cost.func$quadratic$c0 <- 67.699

# 'type' needs to be either 'quadratic' or 'exponential'
cost <- function(h0, dh, cost.func, type) {
  if(dh==0) {
    cost <- 0
  } else {
    if(type=='exponential') {
      cost <- (cost.func$exponential$c + cost.func$exponential$b*dh) * exp(cost.func$exponential$lambda*(h0+dh))
    } else if(type=='quadratic') {
      cost <- cost.func$quadratic$a0*(h0+dh)^2 + cost.func$quadratic$b0*dh + cost.func$quadratic$c0
    } else {
      print('error - unknown investment cost type')
    }
  }
  return(cost)
}

# 10*heightening because heightening in meters, but these cost functions are cm
expected_cost_exp  <- sapply(1:length(heightening), function(i) {cost(h0=0, dh=10*heightening[i], cost.func=cost.func, type='exponential')})
expected_cost_quad <- sapply(1:length(heightening), function(i) {cost(h0=0, dh=10*heightening[i], cost.func=cost.func, type='quadratic'  )})

#
#===============================================================================
# read a previous calibration output file, to use for flood risk
#===============================================================================
#

types.of.gev <- c('gev3','gev4','gev5','gev6')
types.of.nav <- c('nav3','nav4','nav5','nav6')
types.of.gpd <- c('gpd3','gpd4','gpd5','gpd6')
types.of.model <- c(types.of.gev, types.of.nav, types.of.gpd)
nmodel <- length(types.of.model)

# annual block maxima models
yearBM.models <- c(types.of.gev, types.of.nav)

# monthly block maxima models <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< TODO!
#monthBM.models <- c(types.of.gev, types.of.nav)

# daily peaks-over-thresholds models
dayPOT.models <- c(types.of.gpd) # c(types.of.gpd, types.of.nav)  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< TODO!


parameters <- vector('list', nmodel); names(parameters) <- types.of.model
parnames_all <- vector('list', nmodel); names(parnames_all) <- types.of.model
covjump <- vector('list', nmodel); names(covjump) <- types.of.model

# get the annual block maxima model parameters
ncdata <- nc_open(filename.yearBM.parameters)
  time_forc <- ncvar_get(ncdata, 'time')
  temperature_forc <- ncvar_get(ncdata, 'temperature')
  for (model in yearBM.models) {
    parameters[[model]]   <- t(ncvar_get(ncdata, paste('parameters.',model,sep='')))
    parnames_all[[model]] <- ncvar_get(ncdata, paste('parnames.',model,sep=''))
    covjump[[model]]      <- ncvar_get(ncdata, paste('covjump.',model,sep=''))
  }
nc_close(ncdata)

# get the daily peaks-over-thresholds model parameters
ncdata <- nc_open(filename.dayPOT.parameters)
  time_forc <- ncvar_get(ncdata, 'time')
  temperature_forc <- ncvar_get(ncdata, 'temperature')
  for (model in dayPOT.models) {
    parameters[[model]]   <- t(ncvar_get(ncdata, paste('parameters.',model,sep='')))
    parnames_all[[model]] <- ncvar_get(ncdata, paste('parnames.',model,sep=''))
    covjump[[model]]      <- ncvar_get(ncdata, paste('covjump.',model,sep=''))
  }
nc_close(ncdata)

# calibration data, including GPD threshold
data_calib <- readRDS(paste(filename.datacalib,sep=''))


#
#===============================================================================
# read sea-level rise realizations, fingerprinted to Delfzijl TG station
#===============================================================================
#

source('fingerprint_slr_delfzijl.R')

local_sea_level <- get_lsl(filename.sealevelrise, rcp='RCP85', lat=lat,
                           lon=lon, time.beg=time.beg, time.end=time.end,
                           dt=time.step)

#
#===============================================================================
# clip forcing (temperature), and sample subsidence
#===============================================================================
#

time_proj <- local_sea_level$time
lsl_proj <- local_sea_level$lsl
temperature_proj <- temperature_forc[which(time_forc==time_proj[1]):which(time_forc==time_proj[length(time_proj)])]
time_proj_rel <- time_proj - time_proj[1]
lsl_subsidence <- subsidence.rate * time_proj_rel


#
#===============================================================================
# flood risk analysis (Van Dantzig)
#===============================================================================
#

# Set some preliminary stuff so you don't have to keep calculating them

n.ensemble <- rep(NA, nmodel); names(n.ensemble) <- types.of.model
for (model in types.of.model) {n.ensemble[[model]] <- nrow(parameters[[model]])}

n.heightening <- length(heightening)
n.time <- length(time_proj)
names.vandantzig <- c('p_fail_tot','p_fail_avg','p_fail_max','expected_damage','expected_cost_exp','expected_cost_quad','total_loss_exp','total_loss_quad')
length.vandantzig <- length(names.vandantzig)
vandantzig.out <- vector('list', nmodel); names(vandantzig.out) <- types.of.model

h.eff0 <- height.initial - lsl_subsidence - lsl_proj

#===============================================================================
# define the functions for Van Dantzig (Flood risk) outcomes under each model
# structure (this saves a few minutes per model ensemble)
#===============================================================================

# GEV (annual block maxima)
outcome_gev <- function(h, h.eff0, heightening, par,
                        value.initial, hgtdamage.rate, econgrowth.rate, time_proj_rel,
                        discount.rate, expected_cost_exp, expected_cost_quad, names.vandantzig) {
    res <- NULL
    hgt <- 1000*(h.eff0+heightening[h])
    p_fail <- 1-gev_cdf(q=hgt, loc=par.tmp[,1], scale=par.tmp[,2], shape=par.tmp[,3])
    p_fail_tot <- 1 - prod(1-p_fail)
    p_fail_avg <- mean(p_fail)
    p_fail_max <- max(p_fail)
    # heightening costs are not discounted because they occur in present
    expected_damage <- mean( p_fail * value.initial * exp(hgtdamage.rate*heightening[h]*10) *
                                 exp(econgrowth.rate*time_proj_rel) / ((1+discount.rate)^time_proj_rel) )
    total_loss_exp  <- expected_damage + expected_cost_exp[h]
    total_loss_quad <- expected_damage + expected_cost_quad[h]
    res <- c(p_fail_tot, p_fail_avg, p_fail_max, expected_damage, expected_cost_exp[h], expected_cost_quad[h], total_loss_exp, total_loss_quad)
    names(res) <- names.vandantzig
    return(res)
}

# Naveau-i (annual block maxima)
outcome_nav <- function(h, h.eff0, heightening, par,
                        value.initial, hgtdamage.rate, econgrowth.rate, time_proj_rel,
                        discount.rate, expected_cost_exp, expected_cost_quad, names.vandantzig) {
    res <- NULL
    hgt <- 1000*(h.eff0+heightening[h])
    p_fail <- 1-naveau_cdf(x=hgt, kappa=par.tmp[,1], sigma=par.tmp[,2], xi=par.tmp[,3])
    p_fail_tot <- 1 - prod(1-p_fail)
    p_fail_avg <- mean(p_fail)
    p_fail_max <- max(p_fail)
    # heightening costs are not discounted because they occur in present
    expected_damage <- mean( p_fail * value.initial * exp(hgtdamage.rate*heightening[h]*10) *
                                 exp(econgrowth.rate*time_proj_rel) / ((1+discount.rate)^time_proj_rel) )
    total_loss_exp  <- expected_damage + expected_cost_exp[h]
    total_loss_quad <- expected_damage + expected_cost_quad[h]
    res <- c(p_fail_tot, p_fail_avg, p_fail_max, expected_damage, expected_cost_exp[h], expected_cost_quad[h], total_loss_exp, total_loss_quad)
    names(res) <- names.vandantzig
    return(res)
}

# PP-GPD
outcome_gpd <- function(h, h.eff0, heightening, par, threshold, nmax.flood, time.length,
                        value.initial, hgtdamage.rate, econgrowth.rate, time_proj_rel,
                        discount.rate, expected_cost_exp, expected_cost_quad, names.vandantzig) {
    res <- NULL
    hgt <- 1000*(h.eff0+heightening[h])
    p_fail <- ppgpd_overtop(h=hgt, lambda=par[,1], sigma=par[,2], xi=par[,3],
                            threshold=data_calib$gpd$threshold, nmax=nmax.flood, time.length=365.25)
    p_fail_tot <- 1 - prod(1-p_fail)
    p_fail_avg <- mean(p_fail)
    p_fail_max <- max(p_fail)
    # heightening costs are not discounted because they occur in present
    expected_damage <- mean( p_fail * value.initial * exp(hgtdamage.rate*heightening[h]*10) *
                                 exp(econgrowth.rate*time_proj_rel) / ((1+discount.rate)^time_proj_rel) )
    total_loss_exp  <- expected_damage + expected_cost_exp[h]
    total_loss_quad <- expected_damage + expected_cost_quad[h]
    res <- c(p_fail_tot, p_fail_avg, p_fail_max, expected_damage, expected_cost_exp[h], expected_cost_quad[h], total_loss_exp, total_loss_quad)
    names(res) <- names.vandantzig
    return(res)
}

# and get the functions to project the model parameters
source('likelihood_naveau.R')
source('likelihood_gev.R')
source('likelihood_ppgpd.R')

#===============================================================================
# actually do the flood risk analysis

for (model in types.of.model) {

  vandantzig.out[[model]] <- vector('list', n.ensemble[[model]])

  print(paste('starting flood risk assessment for model ',model,'...',sep=''))
  tbeg <- proc.time()
  pb <- txtProgressBar(min=0,max=nrow(parameters[[model]]),initial=0,style=3)

  for (i in 1:n.ensemble[[model]]) {

    if(substr(model,1,3) == 'nav') {
      # project the particular model's parameters across the time horizon
      par.tmp <- project_naveau(parameters=parameters[[model]][i,], parnames=parnames_all[[model]], auxiliary=temperature_proj)
      vandantzig.out[[model]][[i]] <- t(sapply(1:n.heightening, function(h) {outcome_nav(h, h.eff0=h.eff0, heightening=heightening,
                                      par=par.tmp,
                                      value.initial=value.initial, hgtdamage.rate=hgtdamage.rate, econgrowth.rate=econgrowth.rate,
                                      time_proj_rel=time_proj_rel, discount.rate=discount.rate, expected_cost_exp=expected_cost_exp,
                                      expected_cost_quad=expected_cost_quad, names.vandantzig=names.vandantzig)}))
    } else if(substr(model,1,3) == 'gev') {
      # project the particular model's parameters across the time horizon
      par.tmp <- project_gev(parameters=parameters[[model]][i,], parnames=parnames_all[[model]], auxiliary=temperature_proj)
      vandantzig.out[[model]][[i]] <- t(sapply(1:n.heightening, function(h) {outcome_gev(h, h.eff0=h.eff0, heightening=heightening,
                                      par=par.tmp,
                                      value.initial=value.initial, hgtdamage.rate=hgtdamage.rate, econgrowth.rate=econgrowth.rate,
                                      time_proj_rel=time_proj_rel, discount.rate=discount.rate, expected_cost_exp=expected_cost_exp,
                                      expected_cost_quad=expected_cost_quad, names.vandantzig=names.vandantzig)}))
    } else if(substr(model,1,3) == 'gpd') {
      # project the particular model's parameters across the time horizon
      par.tmp <- project_ppgpd(parameters=parameters[[model]][i,], parnames=parnames_all[[model]], auxiliary=temperature_proj)
      vandantzig.out[[model]][[i]] <- t(sapply(1:n.heightening, function(h) {outcome_gpd(h, h.eff0=h.eff0, heightening=heightening,
                                      par=par.tmp, threshold=data_calib$gpd$threshold, nmax.flood=nmax.flood, time.length=time.length,
                                      value.initial=value.initial, hgtdamage.rate=hgtdamage.rate, econgrowth.rate=econgrowth.rate,
                                      time_proj_rel=time_proj_rel, discount.rate=discount.rate, expected_cost_exp=expected_cost_exp,
                                      expected_cost_quad=expected_cost_quad, names.vandantzig=names.vandantzig)}))
    } else {print('error - unrecognized model type')}
    setTxtProgressBar(pb, i)
  }
  close(pb)
  tend <- proc.time()
  print(paste(' ... done. That took ',(tend[3]-tbeg[3])/60,' minutes', sep=''))
}


#
#===============================================================================
# write output file
#===============================================================================
#

# for now just save workspace image
save.image(file='../output/floodrisk_inprog.RData')

#
#===============================================================================
# analysis
#===============================================================================
#

# which cost function to use?
#cost.function <- 'exp'
cost.function <- 'quad'


# For each ensemble, for each ensemble member, calculate the optimal heightening
# Call P(s,x) the performance of strategy s in SOW x, and Popt(x) the best case
iopt <- vector('list', nmodel); names(iopt) <- types.of.model
Hopt <- vector('list', nmodel); names(Hopt) <- types.of.model
Popt <- vector('list', nmodel); names(Popt) <- types.of.model

# Also calculate the ensemble median total cost at each hegihtening, and minimize
# this, in order to develop an ensemble-optimal strategy to follow.
loss.ens <- vector('list', nmodel); names(loss.ens) <- types.of.model
loss.ens.med <- vector('list', nmodel); names(loss.ens.med) <- types.of.model
iopt.ens <- rep(NA, nmodel); names(iopt.ens) <- types.of.model
Hopt.ens <- rep(NA, nmodel); names(Hopt.ens) <- types.of.model
Popt.ens <- rep(NA, nmodel); names(Popt.ens) <- types.of.model

for (model in types.of.model) {
  iopt[[model]] <- rep(NA, n.ensemble[[model]])
  Hopt[[model]] <- rep(NA, n.ensemble[[model]])
  Popt[[model]] <- rep(NA, n.ensemble[[model]])
  loss.ens[[model]] <- mat.or.vec(n.heightening, n.ensemble[[model]])
  loss.ens.med[[model]] <- rep(NA, n.heightening)
  for (i in 1:n.ensemble[[model]]) {
    loss.ens[[model]][,i] <- vandantzig.out[[model]][[i]][,paste('total_loss_',cost.function,sep='')]
    iopt[[model]][i] <- which.min(vandantzig.out[[model]][[i]][,paste('total_loss_',cost.function,sep='')])
    Hopt[[model]][i] <- heightening[iopt[[model]][i]]
    Popt[[model]][i] <- vandantzig.out[[model]][[i]][iopt[[model]][i],paste('total_loss_',cost.function,sep='')]
  }
  loss.ens.med[[model]] <- apply(X=loss.ens[[model]], MARGIN=1, FUN=median)
  iopt.ens[[model]] <- which.min(loss.ens.med[[model]])
  Hopt.ens[[model]] <- heightening[iopt.ens[[model]]]
  Popt.ens[[model]] <- loss.ens.med[[model]][iopt.ens[[model]]]
}

# Regret for strategy s in SOW x is: R(s,x) = Popt(x) - P(s,x)
# So calculate the expected regret of each strategy (heightening) as
# R(s) = int_x{ R(s,x) p(x) dx}
# So following Hopt.ens, how much worse does this perform for each SOW than the
# optimal strategy for that SOW (Hopt[[model]][i]) ?

Reg <- vector('list', nmodel); names(Reg) <- types.of.model
Reg.med <- rep(NA, nmodel); names(Reg.med) <- types.of.model
for (model in types.of.model) {
  Reg[[model]] <- rep(NA, n.ensemble[[model]])
  for (i in 1:n.ensemble[[model]]) {
    Reg[[model]][i] <- vandantzig.out[[model]][[i]][iopt.ens[[model]],paste('total_loss_',cost.function,sep='')] -
                       Popt[[model]][i]
  }
  Reg.med[[model]] <- median(Reg[[model]])
}


# What is distribution of expected regret for each ensemble, if we heighten by
# the ensemble mean/median optimal heightening?

# TODO


#
#===============================================================================
# What is the regret of using model X if model Y is the 'truth'?
#===============================================================================
#

# heighten by Hopt[[model_X]] but calculate regret using model_Y

regret_matrix_avg <- matrix(nrow=nmodel, ncol=nmodel)
rownames(regret_matrix_avg) <- types.of.model
colnames(regret_matrix_avg) <- types.of.model
regret_matrix_med <- matrix(nrow=nmodel, ncol=nmodel)
rownames(regret_matrix_med) <- types.of.model
colnames(regret_matrix_med) <- types.of.model
RegX <- vector('list', nmodel); names(RegX) <- types.of.model
for (model_assumed in types.of.model) {
  RegX[[model_assumed]] <- vector('list', nmodel); names(RegX[[model_assumed]]) <- types.of.model
  for (model_truth in types.of.model) {
    RegX[[model_assumed]][[model_truth]] <- rep(NA, n.ensemble[[model]])
  }
}

for (model_assumed in types.of.model) {
  for (model_truth in types.of.model) {
    RegX[[model_assumed]][[model_truth]] <- loss.ens[[model_truth]][iopt.ens[[model_assumed]],] -
                                            Popt[[model_truth]]
    regret_matrix_avg[model_assumed, model_truth] <- mean(RegX[[model_assumed]][[model_truth]])
    regret_matrix_med[model_assumed, model_truth] <- median(RegX[[model_assumed]][[model_truth]])
  }
}

# need the priors?
priors.dayPOT <- readRDS(filename.priors.dayPOT)
priors.yearBM <- readRDS(filename.priors.yearBM)



# Calculate BIC for each model
# Note: this is only based on the thinned ensembles currently. Can use the full
# calibrated simulations from the MCMC if we want.

bic <- rep(NA, nmodel); names(bic) <- types.of.model
llik.mod <- rep(NA, nmodel); names(llik.mod) <- types.of.model
llik.mod.all <- vector('list', nmodel); names(llik.mod.all) <- types.of.model
lpost.mod.all <- vector('list', nmodel); names(lpost.mod.all) <- types.of.model
for (model in types.of.model) {
  if(substr(model,4,4)!='3') {auxiliary <- trimmed_forcing(data_calib$gev_year$year, time_forc, temperature_forc)$temperature}
  lpri.tmp <- rep(NA, n.ensemble[[model]])
  llik.tmp <- rep(NA, n.ensemble[[model]])
  llik.mod.all[[model]] <- rep(NA, n.ensemble[[model]])
  lpost.mod.all[[model]] <- rep(NA, n.ensemble[[model]])
  for (i in 1:n.ensemble[[model]]) {
    if(model %in% types.of.nav) {
      lpri.tmp[i] <- log_prior_naveau(parameters=parameters[[model]][i,], parnames=parnames_all[[model]], priors=priors.yearBM, model=model)
      llik.tmp[i] <- log_like_naveau(parameters=parameters[[model]][i,], parnames=parnames_all[[model]], data_calib=data_calib$gev_year$lsl_max, auxiliary=auxiliary, Tmax=Tmax)
      ndata <- length(data_calib$gev_year$year)
    } else if(model %in% types.of.gev) {
      lpri.tmp[i] <- log_prior_gev(parameters=parameters[[model]][i,], parnames=parnames_all[[model]], priors=priors.yearBM, model=model)
      llik.tmp[i] <- log_like_gev(parameters=parameters[[model]][i,], parnames=parnames_all[[model]], data_calib=data_calib$gev_year$lsl_max, auxiliary=auxiliary)
      ndata <- length(data_calib$gev_year$year)
    } else if(model %in% types.of.gpd) {
      lpri.tmp[i] <- log_prior_ppgpd(parameters=parameters[[model]][i,], parnames=parnames_all[[model]], priors=priors.dayPOT, model=model)
      llik.tmp[i] <- log_like_ppgpd(parameters=parameters[[model]][i,], parnames=parnames_all[[model]], data_calib=data_calib, auxiliary=auxiliary)
      ndata <- data_calib$gpd$counts_all
    }
    llik.mod.all[[model]][i] <- llik.tmp[i]
    lpost.mod.all[[model]][i] <- llik.tmp[i] + lpri.tmp[i]
  }
  imax <- which.max(llik.tmp)
  bic[[model]] <- -2*max(llik.tmp) + length(parnames_all[[model]])*log(ndata)
  llik.mod[[model]] <- mean(llik.tmp[is.finite(llik.tmp)])
}

# a few different model averaging schemes
# note that wgt_lik is actual BMA
wgt_sf <- rep(NA, nmodel); names(wgt_sf) <- types.of.model
wgt_bic <- rep(NA, nmodel); names(wgt_bic) <- types.of.model
wgt_lik <- rep(NA, nmodel); names(wgt_lik) <- types.of.model

reg_wgt_sf <- rep(NA, nmodel); names(reg_wgt_sf) <- types.of.model
reg_wgt_bic <- rep(NA, nmodel); names(reg_wgt_bic) <- types.of.model
reg_wgt_lik <- rep(NA, nmodel); names(reg_wgt_lik) <- types.of.model

llik.ref <- median(llik.mod)

wgt_bic <- min(bic)/bic; wgt_bic <- wgt_bic/sum(wgt_bic)
wgt_sf  <- exp(-2*(bic-bic[1])) / sum(exp(-2*(bic-bic[1])))
wgt_lik <- exp(llik.mod - llik.ref)/sum(exp(llik.mod - llik.ref))


for (model in types.of.model) {
  reg_wgt_bic[[model]] <- sum(regret_matrix_avg[match(model,types.of.model),]*wgt_bic)
  reg_wgt_lik[[model]] <- sum(regret_matrix_avg[match(model,types.of.model),]*wgt_lik)
  reg_wgt_sf[[model]] <- sum(regret_matrix_avg[match(model,types.of.model),]*wgt_sf)
}


regret_comparison <- cbind(regret_matrix_avg, reg_wgt_lik, reg_wgt_bic, reg_wgt_sf)

#
#===============================================================================
# save analysis output?
#===============================================================================
#

# for now just save workspace image
save.image(file='../output/floodrisk_inprog.RData')

#
#===============================================================================
# some plots
#                 - eventually, stick in the 'plotting.R' routine
#===============================================================================
#

plot.dir <- '../figures/'

#=====================
# BIC comparison, annual block maxima models
#=====================

pdf(paste(plot.dir,'bic_comparison.pdf',sep=''),width=9,height=3.5,colormodel='cmyk')
layout(matrix(c(1,2), nrow=1), c(5.4,2.5), c(3.5, 3.5))
par(mai=c(.8,.8,.5,.2))
barplot(height=bic[1:8], names.arg=types.of.model[1:8], ylim=c(2090,2110),
        main='BIC (lower is better)', ylab='BIC', xlab='Model', xpd=FALSE)
par(mai=c(.8,.2,.5,.2))
barplot(height=bic[9:12], names.arg=types.of.model[9:12], ylim=c(6340,6380),
        main='', ylab='BIC', xlab='', xpd=FALSE)
dev.off()


#=====================
# PP-GPD rate projection
#=====================

# plot of gpd4 rate (based on lambda) of exceedances projected
lambda.proj          <- vector('list', length(types.of.gpd)); names(lambda.proj)          <- types.of.gpd
lambda.proj.quantile <- vector('list', length(types.of.gpd)); names(lambda.proj.quantile) <- types.of.gpd
rate.proj.quantile   <- vector('list', length(types.of.gpd)); names(rate.proj.quantile)   <- types.of.gpd
for (model in types.of.gpd) {
  lambda.proj[[model]] <- matrix(rep(parameters[[model]][,1], n.time), ncol=n.time) +
                          matrix(rep(parameters[[model]][,2], n.time), ncol=n.time) *
                          t(matrix(rep(temperature_proj, n.ensemble[[model]]), ncol=n.ensemble[[model]]))
  lambda.proj.quantile[[model]] <- t(sapply(1:n.time, function(t) {quantile(lambda.proj[[model]][,t], c(0.05, 0.5, .95))}))
  rate.proj.quantile[[model]] <- lambda.proj.quantile[[model]] * 365.25
}


plot(time_proj, rate.quant[,2], type='l', lwd=2, ylim=c(0,12))
lines(time_proj, rate.quant[,1], type='l', lwd=2, lty=2)
lines(time_proj, rate.quant[,3], type='l', lwd=2, lty=2)


pdf(paste(plot.dir,'bic_comparison.pdf',sep=''),width=9,height=3.5,colormodel='cmyk')
layout(matrix(c(1,2), nrow=1), c(5.4,2.5), c(3.5, 3.5))
par(mai=c(.8,.8,.5,.2))
barplot(height=bic[1:8], names.arg=types.of.model[1:8], ylim=c(2090,2110),
        main='BIC (lower is better)', ylab='BIC', xlab='Model', xpd=FALSE)
par(mai=c(.8,.2,.5,.2))
barplot(height=bic[9:12], names.arg=types.of.model[9:12], ylim=c(6340,6380),
        main='', ylab='BIC', xlab='', xpd=FALSE)
dev.off()








# plot of the BMA weights and
par(mfrow=c(2,1))
barplot(height=wgt_lik, names.arg=types.of.model, ylim=c(0, .3), xpd=FALSE,
        ylab='Model weights', xlab='Model', main='BMA weights, from MCMC')
barplot(height=reg_wgt_lik, names.arg=types.of.model, ylim=c(0,70), xpd=FALSE,
        ylab=paste('BMA-weighted expected regret (M',euro,')', sep=''), xlab='Model')



#
#===============================================================================
# scratch below here
#===============================================================================
#

# some of the nav6 SOW have absurdly high regret
plot.dir <- '../figures/'
itmp <- which(Reg$nav6 > 1e5) # note that this is in the trillions of euros
model <- 'nav6'
pdf(paste(plot.dir,'parameter_correlations_',model,'.pdf',sep=''),width=10,height=10,colormodel='cmyk')
par(mfrow=c(5,5))
for (p1 in 1:(length(parnames_all[[model]])-1)) {
  if(p1 > 1) {for (i in 1:(p1-1)) {plot.new();}}
  for (p2 in (p1+1):length(parnames_all[[model]])) {
    plot(parameters[[model]][,p2], parameters[[model]][,p1], xlab=parnames_all[[model]][p2], ylab=parnames_all[[model]][p1])
    points(parameters[[model]][itmp,p2], parameters[[model]][itmp,p1], col='red')
  }
}
dev.off()

# filtering of the naveau models that project any kappa < 0
n.ensemble <- n.ensemble
parameters.filt <- parameters
ibad <- vector('list', length(types.of.nav)); names(ibad) <- types.of.nav
for (model in types.of.nav) {
  ibad.tmp <- NULL
  for (i in 1:n.ensemble[[model]]) {
    par.tmp <- project_naveau(parameters=parameters[[model]][i,], parnames=parnames_all[[model]], auxiliary=temperature_forc)
    if(any(par.tmp[,'kappa'] < 0)) {ibad.tmp <- c(ibad.tmp, i)}
  }
  ibad[[model]] <- ibad.tmp
  if(length(ibad[[model]]) > 0) {
    n.ensemble[[model]] <- n.ensemble[[model]] - length(ibad[[model]])
    parameters.filt[[model]] <- parameters[[model]][-ibad[[model]],]
  }
}


# testing GEV cdf timing
np <- 50
par.tmp <- parameters$gev3[1:np,]
niter <- 10000
tbeg <- proc.time()
pb <- txtProgressBar(min=0,max=niter,initial=0,style=3)
for ( i in 1:niter) {
  p_fail_1 <- 1-gev_cdf(q=1000*height.eff[1:np], loc=par.tmp[1:np,1], scale=par.tmp[1:np,2], shape=par.tmp[1:np,3])
  #p_fail_2 <- 1-sapply(1:n.time, function(t) {pevd(q=1000*height.eff[t], loc=par.tmp[t,1], scale=par.tmp[t,2], shape=par.tmp[t,3])})
  setTxtProgressBar(pb, i)
}
close(pb)
tend <- proc.time()
t_gevcdf <- tend-tbeg

tbeg <- proc.time()
pb <- txtProgressBar(min=0,max=niter,initial=0,style=3)
for ( i in 1:niter) {
  #p_fail_1 <- 1-gev_cdf(q=1000*height.eff, loc=par.tmp[,1], scale=par.tmp[,2], shape=par.tmp[,3])
  p_fail_2 <- 1-sapply(1:np, function(t) {pevd(q=1000*height.eff[t], loc=par.tmp[t,1], scale=par.tmp[t,2], shape=par.tmp[t,3])})
  setTxtProgressBar(pb, i)
}
close(pb)
tend <- proc.time()
t_sapply <- tend-tbeg


# test plot
if(FALSE) {
model <- 'nav5'
plot(heightening, vandantzig.out[[model]][[1]][,'total_loss_exp'], pch=16, ylim=c(0,300), xlim=c(0,4), xlab='Heightening [m]', ylab=paste('Expected costs [M',euro,']',sep=''))
points(heightening, vandantzig.out[[model]][[1]][,'expected_damage'], col='red', pch=16)
points(heightening[1], vandantzig.out[[model]][[1]][1,'total_loss_exp'], pch=16)
points(heightening, expected_cost_exp, col='blue', pch=16)
imin <- which.min(vandantzig.out[[model]][[1]][,'total_loss_exp'])
points(heightening[imin], vandantzig.out[[model]][[1]][imin,'total_loss_exp'], pch=8, col='purple', lwd=3, cex=1.5)
legend(1,300, c('total costs','build costs','damages'), pch=c(16,16,16), col=c('black','blue','red'), bty='n')
}



}



#===============================================================================
# End
#===============================================================================

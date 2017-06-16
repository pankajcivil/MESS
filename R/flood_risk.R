#===============================================================================
#
# questions? Tony Wong (twong@psu.edu)
#===============================================================================

#
#===============================================================================
# settings
#===============================================================================
#

filename.evt.parameters <- 'evt_models_calibratedParameters__13Jun2017.nc'
filename.sealevelrise <- '../../BRICK/output_model/BRICK-fastdyn_physical_gamma_01Jun2017.nc'

time.beg <- 2015           # inital year, "present"
time.end <- 2065           # final year (time horizon)
time.step <- 1             # time step in years
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
p0 <- 0.00137              # Eijgenraam et al (2012)
gev.prelim <- amcmc_prelim$gev3$samples[which.max(amcmc_prelim$gev3$log.p),]
height.initial <- qevd(p=1-p0, loc=gev.prelim[1], scale=gev.prelim[2], shape=gev.prelim[3])/1000

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
cost.exp  <- sapply(1:length(heightening), function(i) {cost(h0=0, dh=10*heightening[i], cost.func=cost.func, type='exponential')})
cost.quad <- sapply(1:length(heightening), function(i) {cost(h0=0, dh=10*heightening[i], cost.func=cost.func, type='quadratic'  )})

#
#===============================================================================
# read a previous calibration output file, to use for flood risk
#===============================================================================
#

types.of.gev <- c('gev3','gev4','gev5','gev6')
types.of.nav <- c('nav3','nav4','nav5','nav6')
types.of.model <- c(types.of.gev, types.of.nav)

parameters <- vector('list', length(types.of.model)); names(parameters) <- types.of.model
parnames_all <- vector('list', length(types.of.model)); names(parnames_all) <- types.of.model
covjump <- vector('list', length(types.of.model)); names(covjump) <- types.of.model
ncdata <- nc_open(filename.evt.parameters)
  time_forc <- ncvar_get(ncdata, 'time')
  temperature_forc <- ncvar_get(ncdata, 'temperature')
  for (model in types.of.model) {
    parameters[[model]]   <- t(ncvar_get(ncdata, paste('parameters.',model,sep='')))
    parnames_all[[model]] <- ncvar_get(ncdata, paste('parnames.',model,sep=''))
    covjump[[model]]      <- ncvar_get(ncdata, paste('covjump.',model,sep=''))
  }
nc_close(ncdata)


# TODO
# TODO
         # thinning for now, to test everything
         thin.len <- 1000
         for (model in types.of.model) {
           parameters[[model]] <- parameters[[model]][seq(from=1, to=nrow(parameters[[model]]), by=thin.len),]
         }
# TODO
# TODO


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

n.ensemble <- rep(NA, length(types.of.model)); names(n.ensemble) <- types.of.model
for (model in types.of.model) {n.ensemble[[model]] <- nrow(parameters[[model]])}
n.heightening <- length(heightening)
n.time <- length(time_proj)
names.vandantzig <- c('p_fail_tot','p_fail_avg','p_fail_max','expected_damage','expected_cost_exp','expected_cost_quad','total_loss_exp','total_loss_quad')
length.vandantzig <- length(names.vandantzig)
vandantzig.out <- vector('list', length(types.of.model)); names(vandantzig.out) <- types.of.model

for (model in types.of.model) {

  vandantzig.out[[model]] <- vector('list', n.ensemble[[model]])

  print(paste('starting flood risk assessment for model ',model,'...',sep=''))
  pb <- txtProgressBar(min=0,max=nrow(parameters[[model]]),initial=0,style=3)

  for (i in 1:n.ensemble[[model]]) {

    vandantzig.out[[model]][[i]] <- mat.or.vec(n.heightening, length.vandantzig)
    colnames(vandantzig.out[[model]][[i]]) <- names.vandantzig
    p_fail_tot <- rep(0, n.heightening)
    p_fail_tot <- rep(0, n.heightening)
    p_fail_tot <- rep(0, n.heightening)

# TODO -- make this loop over heightenings either an 'apply' or more efficient?

    if(substr(model,1,3) == 'nav') {
      # project the particular model's parameters across the time horizon
      par.tmp <- project_naveau(parameters=parameters[[model]][i,], parnames=parnames_all[[model]], auxiliary=temperature_proj)
      for (h in 1:n.heightening) {
        height.eff <- height.initial + heightening[h] - lsl_subsidence - lsl_proj
        p_fail <- 1-naveau_cdf(x=1000*height.eff, kappa=par.tmp[,1], sigma=par.tmp[,2], xi=par.tmp[,3])
        p_fail_tot[h] <- 1 - prod(1-p_fail)
        p_fail_avg[h] <- mean(p_fail)
        p_fail_max[h] <- max(p_fail)
        # heightening costs are not discounted because they occur in present
        expected_damage[h] <- mean( p_fail * value.initial * exp(hgtdamage.rate*heightening[h]*10) *
                                    exp(econgrowth.rate*time_proj_rel) / ((1+discount.rate)^time_proj_rel) )
        total_loss_exp[h]  <- expected_damage[h] + cost.exp[h]
        total_loss_quad[h] <- expected_damage[h] + cost.quad[h]
      }
    } else if(substr(model,1,3) == 'gev') {
      # project the particular model's parameters across the time horizon
      par.tmp <- project_gev(parameters=parameters[[model]][i,], parnames=parnames_all[[model]], auxiliary=temperature_proj)
      for (h in 1:n.heightening) {
        height.eff <- height.initial + heightening[h] - lsl_subsidence - lsl_proj
        p_fail <- 1-gev_cdf(q=1000*height.eff, loc=par.tmp[,1], scale=par.tmp[,2], shape=par.tmp[,3])
        p_fail_tot[h] <- 1 - prod(1-p_fail)
        p_fail_avg[h] <- mean(p_fail)
        p_fail_max[h] <- max(p_fail)
        # heightening costs are not discounted because they occur in present
        expected_damage[h] <- mean( p_fail * value.initial * exp(hgtdamage.rate*heightening[h]*10) *
                                    exp(econgrowth.rate*time_proj_rel) / ((1+discount.rate)^time_proj_rel) )
        total_loss_exp[h]  <- expected_damage[h] + cost.exp[h]
        total_loss_quad[h] <- expected_damage[h] + cost.quad[h]
      }
    } else {print('error - unrecognized model type')}

# TODO -- there is DEFINITELY a less embarassing way to do this...

    vandantzig.out[[model]][[i]][,'p_fail_tot'] <- p_fail_tot
    vandantzig.out[[model]][[i]][,'p_fail_avg'] <- p_fail_avg
    vandantzig.out[[model]][[i]][,'p_fail_max'] <- p_fail_max
    vandantzig.out[[model]][[i]][,'expected_damage'] <- expected_damage
    vandantzig.out[[model]][[i]][,'expected_cost_exp'] <- cost.exp
    vandantzig.out[[model]][[i]][,'expected_cost_quad'] <- cost.quad
    vandantzig.out[[model]][[i]][,'total_loss_exp'] <- total_loss_exp
    vandantzig.out[[model]][[i]][,'total_loss_quad'] <- total_loss_quad
    setTxtProgressBar(pb, i)
  }
  close(pb)
}

#
#===============================================================================
# write output file
#===============================================================================
#

# TODO

#
#===============================================================================
# analysis
#===============================================================================
#

# For each ensemble, for each ensemble member, calculate the optimal heightening
# Call P(s,x) the performance of strategy s in SOW x, and Popt(x) the best case

# TODO

# Regret for strategy s in SOW x is: R(s,x) = Popt(x) - P(s,x)
# So calculate the expected regret of each strategy (heightening) as
# R(s) = int_x{ R(s,x) p(x) dx}

# TODO

# What is distribution of expected regret for each ensemble, if we heighten by
# the ensemble mean/median optimal heightening?

# TODO

#
#===============================================================================
# save analysis output?
#===============================================================================
#

# TODO


#
#===============================================================================
# scratch below here
#===============================================================================
#

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
points(heightening, cost.exp, col='blue', pch=16)
imin <- which.min(vandantzig.out[[model]][[1]][,'total_loss_exp'])
points(heightening[imin], vandantzig.out[[model]][[1]][imin,'total_loss_exp'], pch=8, col='purple', lwd=3, cex=1.5)
legend(1,300, c('total costs','build costs','damages'), pch=c(16,16,16), col=c('black','blue','red'), bty='n')
}



}



#===============================================================================
# End
#===============================================================================

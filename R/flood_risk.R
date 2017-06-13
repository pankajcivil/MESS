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
#filename.sealevelrise <- TODO

time.beg <- 2015           # inital year, "present"
time.end <- 2065           # final year (time horizon)
height.low <- 0            # lowest heightening considered [m]
height.high <- 10          # tallest heightening considered [m]
height.increment <- 0.05   # increments of heightening considered [m]
heightening <- seq(from=height.low, to=height.high, by=height.increment)

# TODO
height.initial <- 5        # initial dike ring height [m]


#
#===============================================================================
# additional parameters pertaining to flood risk
# (for dike ring 15, Eijgenraam et al 2012)
#===============================================================================
#

# TODO - sample, or hold constant?

subsidence.rate <- 0.002     # m/yr (Rietveld H. Land subsidence in the Netherlands. Pp. 455–465 in Proceedings of the 3rd International Symposium on Land Subsidence. Vol 151. 1986.)
discount.rate <- 0.04        # % (Eijgenraam et al 2012)
value.initial <- 11810.4     #

# TODO - economic growth? increased potential damage with heightening?

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

#
#===============================================================================
# read sea-level rise realizations, fingerprinted to Delfzijl TG station
#===============================================================================
#

# TODO



#
#===============================================================================
# clip forcing (temperature) and project EV distribution parameters over time
# horizon
#===============================================================================
#

# TODO






#
#===============================================================================
# write output file
#===============================================================================
#


#===============================================================================
# End
#===============================================================================

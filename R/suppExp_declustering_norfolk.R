#===============================================================================
# suppExp_declustering_norfolk.R
#
# * Supporting experiment assessing the implications of the assumption of 3-day
#   declustering time scale
# * Process tide gauge data for Norfolk, for all 4 candidate models, using a range
#   of declustering time scales.
# * Calibrate (maximum likelihood)
#
# Questions? Tony Wong (anthony.e.wong@colorado.edu)
#===============================================================================
# Copyright 2017 Tony Wong
#
# MESS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# MESS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# MESS.  If not, see <http://www.gnu.org/licenses/>.
#===============================================================================

# some preliminaries
station <- 'norfolk'
output.dir <- '../output/'
dat.dir <- '../data/'

setwd('/Users/tony/codes/EVT/R')

today <- Sys.Date(); today=format(today,format="%d%b%Y")

appen <- paste('suppExp_declustering',station,today,sep='_')




# install some preliminary packages, or load the libraries if you already have them
l.installpackages <- FALSE
if(l.installpackages) {
  install.packages('date')
  install.packages('zoo')
  install.packages('Hmisc')
  install.packages('ncdf4')
  install.packages('extRemes')
  install.packages('adaptMCMC')
  install.packages('lhs')
  install.packages('DEoptim')
}
library(date)
library(zoo)
library(Hmisc)
library(ncdf4)
library(extRemes)
library(adaptMCMC)
library(DEoptim)


#===============================================================================
# read global mean temperature data, as covariate for PP/GPD parameters
source('read_data_temperature.R')
#===============================================================================


#===============================================================================
# read time series declustering routine (to filter events that are too close to
# one another temporally to be considered independent events)
source('decluster_timeseries.R')
#===============================================================================


#===============================================================================
# Use Norfolk data set (so global mean temperature is analogous assumption to
# Grinsted et al (2013)), and we create data sets using 1-, 3-, 5-, 7-, and 9-day
# declustering time scales
source('processing_norfolk.R')

dt.decluster.test <- c(1,3,5,7,9,11)

# check which files we already have
files.decluster <- list.files('../data/','tidegauge_processed_norfolk_decl')
have.decluster <- rep(FALSE, length(dt.decluster.test))
for (ii in 1:length(dt.decluster.test)){
  this.dt.decluster <- strtoi(substr(files.decluster[ii],33,33))
  ifound <- which(dt.decluster.test==this.dt.decluster)
  if(length(ifound)>0) {have.decluster[ifound]=TRUE}
}
need.decluster <- !have.decluster

for (dt.decluster in dt.decluster.test[need.decluster]) {
  processing_norfolk(dt.decluster)
}

files.decluster <- NULL
names.decluster <- NULL
for (ii in 1:length(dt.decluster.test)) {
  names.decluster <- c(names.decluster, paste('dt',dt.decluster.test[ii],sep=''))
  files.decluster <- c(files.decluster, list.files('../data/',paste('tidegauge_processed_norfolk_decl',dt.decluster.test[ii],sep='')))
}
#===============================================================================


#===============================================================================
# Run differential evolution optimization to obtain maximum likelihood parameter
# sets for each candidate model, using these declustering time scales

# the likelihood function for PP/GPD is in here
source('likelihood_ppgpd.R')
# DE optim can only minimize, so need to minimize neg_log_like_ppgpd (which is
# conveniently a function in likelihood_ppgpd.R)

# and set up the parameters for PP/GPD model
source('parameter_setup_dayPOT.R')

# parameters for DE optim
NP.deoptim <- 200
niter.deoptim <- 200
F.deoptim <- 0.8
CR.deoptim <- 0.9

# create data object to hold the bic, dic etc. error metrics (BIC for now)
evaluation_metrics <- c('bic','maxlik')
metrics <- vector('list', 2); names(metrics) <- evaluation_metrics
for (mm in evaluation_metrics) {
  # want a [num decluster] x [num model] table
  metrics[[mm]] <- mat.or.vec(length(dt.decluster.test), nmodel)
  rownames(metrics[[mm]]) <- names.decluster
  colnames(metrics[[mm]]) <- types.of.gpd
}

# create object to hold the calibrated parameters and (negative) log-likelihood
deoptim.decl <- vector('list', nmodel); names(deoptim.decl) <- types.of.model
for (gpd.type in types.of.model) {
  deoptim.decl[[gpd.type]] <- mat.or.vec(length(dt.decluster.test), length(parnames_all[[gpd.type]])+1)
  rownames(deoptim.decl[[gpd.type]]) <- names.decluster
  colnames(deoptim.decl[[gpd.type]]) <- c(parnames_all[[gpd.type]], 'nnlike')
}


# which declustering time scale experiment to use?
for (ii in 1:length(dt.decluster.test)) {
  filename.datacalib <- paste(dat.dir,files.decluster[ii],sep='')
  print(paste('Reading calibration data file: ',filename.datacalib,' ...',sep=''))
  data_calib <- readRDS(filename.datacalib)

  for (gpd.type in types.of.gpd) {
    print(paste(' - starting DE optimization for model ',gpd.type,'...',sep=''))
    if(gpd.type=='gpd3') {auxiliary <- NULL
    } else {auxiliary <- trimmed_forcing(data_calib$gpd$year, time_forc, temperature_forc)$forcing}
    # if tide gauge record starts before temperatures, clip it
    if(data_calib$gpd$year[1] < time_forc[1]) {
      irem <- which(data_calib$gpd$year < time_forc[1])
      auxiliary <- auxiliary[-irem]
      data_calib$gev_year$year <- data_calib$gev_year$year[-irem]
      data_calib$gev_year$lsl_max <- data_calib$gev_year$lsl_max[-irem]
      data_calib$gpd$year <- data_calib$gpd$year[-irem]
      data_calib$gpd$counts <- data_calib$gpd$counts[-irem]
      data_calib$gpd$excesses <- data_calib$gpd$excesses[-irem]
      data_calib$gpd$time_length <- data_calib$gpd$time_length[-irem]
      data_calib$gpd$time_length_all <- sum(data_calib$gpd$time_length)
      data_calib$gpd$counts_all <- sum(unlist(data_calib$gpd$counts), na.rm=TRUE)
      data_calib$gpd$excesses_all <- unlist(data_calib$gpd$excesses)[!is.na(unlist(data_calib$gpd$excesses))]
    }
    out.deoptim <- DEoptim(neg_log_like_ppgpd, lower=bound_lower_set[[gpd.type]], upper=bound_upper_set[[gpd.type]],
                         DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                         parnames=parnames_all[[gpd.type]], data_calib=data_calib$gpd, auxiliary=auxiliary)

    # save the output parameters, and likelihood values
    deoptim.decl[[gpd.type]][ii,] <- c(out.deoptim$optim$bestmem, out.deoptim$optim$bestval)

    metrics$maxlik[ii,gpd.type] <- -deoptim.decl[[gpd.type]][ii,length(parnames_all[[gpd.type]])+1]
    metrics$bic[ii,gpd.type] <- -2*metrics$maxlik[ii,gpd.type] + length(parnames_all[[gpd.type]])*log(length(data_calib$gpd$counts_all))
  }
}

# save progress
filename.saveprogress <- '../output/declustering_experiment_inprogress.RData'
save.image(file=filename.saveprogress)

#===============================================================================


#===============================================================================

# TODO
# TODO

# problem with BIC (and AIC, DIC, etc) is that the data set used is a bit
# different from one experiment to the next.

#===============================================================================
# calculate the empirical return levels, and compare against those from the
# MLE parameters for each model and

# TODO


# return level calculation (from extRemes package):     (100 is RP (years), xi, sigma, lambda and threshold are scalar)
#rl100[[site]][[data.len]][[model]][[year]][sow] <- rlevd(100, scale=sigma, shape=xi,
#                                                         threshold=data.sites[[site]]$gpd$threshold,
#                                                         type='GP',
#                                                         npy=365.25,
#                                                         rate=lambda)


# empirical return levels:
#esf.levels.eur <- data_set[[dd]]$lsl_max[order(data_set[[dd]]$lsl_max)]
#esf.values.eur <- seq(from=length(data_set[[dd]]$lsl_max), to=1, by=-1)/(length(data_set[[dd]]$lsl_max)+1)
#
# but that's for annual block maxima... what about the scattered, irregular data for PP/GPD model?

#===============================================================================



#===============================================================================
# Generate some figures for SOM



#
#===============================================================================
# End
#===============================================================================
#

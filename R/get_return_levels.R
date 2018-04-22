#===============================================================================
# analysis_and_plotting.R
#
# Analysis and plotting for wong, klufas and keller pp/gpd analysis.
# This assumes you've downloaded the Github repo and are running from the
# 'R' directory.
#
# If any this is vague/confusing and you have questions, email me. Seriously.
#
# Questions? Tony Wong (anthony.e.wong@colorado.edu)
#===============================================================================
#
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

rm(list=ls())

# Some preliminaries

setwd('~/codes/EVT/R')

library(ncdf4)
library(extRemes)
library(Bolstad)

# get some colorblind-friendly colors to plot with
source('colorblindPalette.R')

# set useful directories -- assumes you are in the 'R' directory within the repo
# directory structure
plot.dir <- '../figures/'
output.dir <- '../output/'

# calibrated parameter sets (samples; all should be same size)
# these are the results from 'calibration_dayPOT-experiments_driver.R'
filename.norfolk.normalgamma <-  paste(output.dir,'calibratedParameters_ppgpd-experiments_norfolk_normalgamma_decl3-pot99_31Mar2018.nc', sep='')
filename.delfzijl.normalgamma <- paste(output.dir,'calibratedParameters_ppgpd-experiments_delfzijl_normalgamma_decl3-pot99_31Mar2018.nc',sep='')
filename.norfolk.uniform <-      paste(output.dir,'calibratedParameters_ppgpd-experiments_norfolk_uniform_decl3-pot99_31Mar2018.nc',     sep='')
filename.delfzijl.uniform <-     paste(output.dir,'calibratedParameters_ppgpd-experiments_delfzijl_uniform_decl3-pot99_31Mar2018.nc',    sep='')

filename.norfolk.data <- '../data/tidegauge_processed_norfolk_decl3-pot99-annual_06Dec2017.rds'
filename.delfzijl.data <- '../data/tidegauge_processed_delfzijl_decl3-pot99-annual_20Dec2017.rds'

filename.priors <- '../output/surge_priors_normalgamma_ppgpd_decl3-pot99_28Mar2018.rds'

# file to save progress as you run
filename.saveprogress <- '../output/returnperiod_analysis.RData'

#===============================================================================



#===============================================================================
#===============================================================================
# ANALYSIS
#===============================================================================
#===============================================================================




#===============================================================================
# read NAO data for projecting surge levels covarying with NAO index (DJF)
#===============================================================================

# yields 'nao_forc' and 'time_forc', between 1850 and 2100
# both historical and future NAO index standardized relative to 2001-2016 mean/stdev

source('read_data_naoindex.R')

#===============================================================================
# read parameter sets
#===============================================================================

# create list objects to store the parameters
site.names <- c('Delfzijl','Norfolk'); nsites <- length(site.names)
types.of.gpd <- c('gpd3','gpd4','gpd5','gpd6'); nmodel <- length(types.of.gpd)
# return periods to calculate
n_rp <- 100
return_periods <- 10^seq(from=log10(2), to=log10(500), length.out=n_rp)

# years to grab the return levels
rl.years <- c(2016, 2065); nyears <- length(rl.years)
year.names <- rep(NA, length(rl.years)); for (year in 1:length(rl.years)) {year.names[year] <- paste('y',rl.years[year],sep='')}

# take the 11-year mean, centered at the year in question
nao.years <- NULL
for (year in rl.years) {
  nao.years <- c(nao.years, mean(nao_forc[which(abs(time_forc-year)<=5)]))
}
names(nao.years) <- year.names

gpd.parameters <- vector('list', nsites); names(gpd.parameters) <- site.names
for (site in site.names) {
  gpd.parameters[[site]] <- vector('list', nmodel); names(gpd.parameters[[site]]) <- types.of.gpd
}

# parameter nmaes
parnames <- vector('list', nmodel); names(parnames) <- types.of.gpd

data.sites <- vector('list', nsites); names(data.sites) <- site.names

for (site in site.names) {
  if(site=='Norfolk') {
    ncdata <- nc_open(filename.norfolk.normalgamma)
    for (model in types.of.gpd) {
      gpd.parameters[[site]][[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd89.',model,sep='')))
      parnames[[model]] <- ncvar_get(ncdata, paste('parnames.',model,sep=''))
      data.sites[[site]] <- readRDS(filename.norfolk.data)
    }
    n.ensemble <- nrow(gpd.parameters[[site]]$gpd3)
  } else if(site=='Delfzijl') {ncdata <- nc_open(filename.delfzijl.normalgamma)
    for (model in types.of.gpd) {
      gpd.parameters[[site]][[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd137.',model,sep='')))
      data.sites[[site]] <- readRDS(filename.delfzijl.data)
    }
    n.ensemble <- nrow(gpd.parameters[[site]]$gpd3)
  }
}

# initialize a giant data structure to hold all the return levels, for all the
# return periods, for both sites, for both time horizons, for all 4 candidate
# models
returnlevels <- vector('list', nsites); names(returnlevels) <- site.names
for (site in site.names) {
  returnlevels[[site]] <- vector('list', nmodel); names(returnlevels[[site]]) <- types.of.gpd
  for (model in types.of.gpd) {
    returnlevels[[site]][[model]] <- vector('list', nyears); names(returnlevels[[site]][[model]]) <- year.names
    for (year in year.names) {
      returnlevels[[site]][[model]][[year]] <- mat.or.vec(n_rp, n.ensemble)
    }
  }
}


#===============================================================================
# calculate return levels
#===============================================================================

# each element of returnlevels[[site]][[model]][[year]][[p]] is a vector of length n_ensemble
for (site in site.names) {
  for (model in types.of.gpd) {
    for (year in year.names) {
      print(paste('calculating return levels...',site,model,year,sep=' - '))
      tbeg <- proc.time()
      pb <- txtProgressBar(min=0,max=n.ensemble, initial=0,style=3)
      for (sow in 1:n.ensemble) {
        if (length(parnames[[model]])==3) {
          lambda <- gpd.parameters[[site]][[model]][sow,match('lambda', parnames[[model]])]
          sigma <- gpd.parameters[[site]][[model]][sow,match('sigma', parnames[[model]])]
          xi <- gpd.parameters[[site]][[model]][sow,match('xi', parnames[[model]])]
        } else if(length(parnames[[model]])==4) {
          lambda <- gpd.parameters[[site]][[model]][sow,match('lambda0', parnames[[model]])] +
                    gpd.parameters[[site]][[model]][sow,match('lambda1', parnames[[model]])]*nao.years[[year]]
          sigma <- gpd.parameters[[site]][[model]][sow,match('sigma', parnames[[model]])]
          xi <- gpd.parameters[[site]][[model]][sow,match('xi', parnames[[model]])]
        } else if(length(parnames[[model]])==5) {
          lambda <- gpd.parameters[[site]][[model]][sow,match('lambda0', parnames[[model]])] +
                    gpd.parameters[[site]][[model]][sow,match('lambda1', parnames[[model]])]*nao.years[[year]]
          sigma <- exp(gpd.parameters[[site]][[model]][sow,match('sigma0', parnames[[model]])] +
                       gpd.parameters[[site]][[model]][sow,match('sigma1', parnames[[model]])]*nao.years[[year]])
          xi <- gpd.parameters[[site]][[model]][sow,match('xi', parnames[[model]])]
        } else if(length(parnames[[model]])==6) {
          lambda <- gpd.parameters[[site]][[model]][sow,match('lambda0', parnames[[model]])] +
                    gpd.parameters[[site]][[model]][sow,match('lambda1', parnames[[model]])]*nao.years[[year]]
          sigma <- exp(gpd.parameters[[site]][[model]][sow,match('sigma0', parnames[[model]])] +
                       gpd.parameters[[site]][[model]][sow,match('sigma1', parnames[[model]])]*nao.years[[year]])
          xi <- gpd.parameters[[site]][[model]][sow,match('xi0', parnames[[model]])] +
                gpd.parameters[[site]][[model]][sow,match('xi1', parnames[[model]])]*nao.years[[year]]
        }
        returnlevels[[site]][[model]][[year]][,sow] <- rlevd(return_periods, scale=sigma, shape=xi,
                                                                   threshold=data.sites[[site]]$gpd$threshold,
                                                                   type='GP',
                                                                   npy=365.25,
                                                                   rate=lambda)
        setTxtProgressBar(pb, sow)
      }
      close(pb)
    }
  }
}

#===============================================================================
# Bayesian model averaging ensemble
#===============================================================================

# bma weight results object
bma_weights <- readRDS('../output/bma_weights_threshold99.rds')

# initialize BMA ensemble
returnlevels.bma <- vector('list', nsites); names(returnlevels.bma) <- site.names
for (site in site.names) {
  returnlevels.bma[[site]] <- vector('list', nyears); names(returnlevels.bma[[site]]) <- year.names
  for (year in year.names) {
      returnlevels.bma[[site]][[year]] <- mat.or.vec(n_rp, n.ensemble)
  }
}


# calculate BMA-weighted ensemble using the ensemble members already drawn
# have also checked this making larger (up to 1e6) and separate samples for each
# of the 4 models, and the quantiles do not change by more than a millimeter
for (site in site.names) {
  if (site=='Delfzijl') {dd <- 6}
  if (site=='Norfolk') {dd <- 4}
  for (year in year.names) {
    returnlevels.bma[[site]][[year]] <- bma_weights[[site]][[dd]]['gpd3']*returnlevels[[site]]$gpd3[[year]] +
                                        bma_weights[[site]][[dd]]['gpd4']*returnlevels[[site]]$gpd4[[year]] +
                                        bma_weights[[site]][[dd]]['gpd5']*returnlevels[[site]]$gpd5[[year]] +
                                        bma_weights[[site]][[dd]]['gpd6']*returnlevels[[site]]$gpd6[[year]]
  }
}


#===============================================================================
# Calculate some quantiles of interest for plotting
#===============================================================================

# calculate particular quantiles
quantiles_to_get <- c(0.05, 0.5, 0.95)
quantile.names <- rep(NA, length(quantiles_to_get))
for (q in 1:length(quantiles_to_get)) {quantile.names[q] <- paste('q',quantiles_to_get[q]*100,sep='')}

# names of models and bma
types.of.model <- c(types.of.gpd, 'bma')


returnlevels.q <- vector('list', nsites); names(returnlevels.q) <- site.names
for (site in site.names) {
  returnlevels.q[[site]] <- vector('list', nmodel+1); names(returnlevels.q[[site]]) <- types.of.model
  for (model in types.of.model) {
    returnlevels.q[[site]][[model]] <- vector('list', nyears); names(returnlevels.q[[site]][[model]]) <- year.names
    for (year in year.names) {
      returnlevels.q[[site]][[model]][[year]] <- mat.or.vec(n_rp, length(quantile.names))
      colnames(returnlevels.q[[site]][[model]][[year]]) <- quantile.names
      rownames(returnlevels.q[[site]][[model]][[year]]) <- return_periods
    }
  }
}

#returnlevels.q$Delfzijl$gpd3$y2016

for (site in site.names) {
  for (model in types.of.gpd) {
    for (year in year.names) {
      for (rp in 1:n_rp) {
        returnlevels.q[[site]][[model]][[year]][rp,] <- quantile(returnlevels[[site]][[model]][[year]][rp,], quantiles_to_get)
      }
    }
  }
  model <- 'bma'
  for (year in year.names) {
    for (rp in 1:n_rp) {
      returnlevels.q[[site]][[model]][[year]][rp,] <- quantile(returnlevels.bma[[site]][[year]][rp,], quantiles_to_get)
    }
  }
}


# save progress
save.image(file=filename.saveprogress)

#===============================================================================





#===============================================================================
# End
#===============================================================================

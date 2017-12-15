#===============================================================================
# processing_script.R
#
# Process tide gauge data for three experiments:
# 1. GEV and Naveau-(i) (aka Papastathopoulos and Tawn iii) to annual block maxima
# 2. GEV and Naveau-(i) to monthly block maxima (preprocessing needed)
# 3. PP-GPD and Naveau-(i) to POT daily maxima (preprocessing needed for Naveau,
#    and included/not included in two separate experiments for POT/GPD)
#    (POT = Peaks Over Thresholds)
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================
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

#===============================================================================
# set things

.dt.decluster <- 3           # declustering time scale (days)
.detrend.method <- 'annual'  # annual means or linear detrending?
.pot.threshold <- .99        # percentile for POT threshold


if(Sys.info()['nodename']=='Tonys-MacBook-Pro.local') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('/Users/tony/codes/EVT/R')
  .Ncore <- 0
} else {
  # assume on Napa cluster
  machine <- 'remote'
  setwd('/home/scrim/axw322/codes/EVT/R')
  .Ncore <- 15  # use multiple cores to process the many tide gauge stations
}
#===============================================================================

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
  install.packages('foreach')
  install.packages('doParallel')
}
library(date)
library(zoo)
library(Hmisc)
library(ncdf4)
library(extRemes)
library(adaptMCMC)
library(DEoptim)
library(foreach)
library(doParallel)

# read global mean temperature data, as covariate for PP/GPD parameters
source('read_data_temperature.R')

# read declustering routine that you'll need to process the data files
source('decluster_timeseries.R')


#
#===============================================================================
# Read tide gauge data for in-depth analyses at three sites, and for 27 other
# stations that will be used to derive prior probability distributions for the
# model parameters.
#===============================================================================
#


# many long record stations
if(machine=='local') {
  source('processing_many_stations.R')
} else if(machine=='remote') {
  source('processing_many_stations_parallel.R')
} else {
  print('ERROR: unknown machine')
}
processing_many_stations(dt.decluster=.dt.decluster, detrend.method='annual', pot.threshold=0.99, Ncore=.Ncore)

# Delfzijl, the Netherlands
source('processing_delfzijl.R')
processing_delfzijl(dt.decluster=.dt.decluster, detrend.method='annual', pot.threshold=0.99)

# Norfolk, Virgina, United State
source('processing_norfolk.R')
processing_norfolk(dt.decluster=.dt.decluster, detrend.method='annual', pot.threshold=0.99)

# Balboa, Panama
source('processing_balboa.R')
processing_balboa(dt.decluster=.dt.decluster, detrend.method='annual', pot.threshold=0.99)

#
#===============================================================================
# End
#===============================================================================
#

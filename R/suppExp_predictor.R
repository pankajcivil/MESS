#===============================================================================
# suppExp_predictor.R
#
# * Supporting experiment assessing the implications of using the global annual
#   mean temperature as the pp/gpd model covariate.
# * Grinsted et al (2013) found that local temperature had higher odds ratio
#   (their Table 2), so we check that.
# * Process tide gauge data for Norfolk site, for all 4 candidate models,
#   and look at both monthly mean time series (local) and
# * Calibrate (maximum likelihood), and calculate the BIC for each of the
#   temperature time series
# * Experiments:
#     1) global annual mean temperature (control)
#     2) global monthly mean temperature
#     3) local annual mean temperature
#     4) local monthly mean temperature
#     5) time
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

rm(list=ls())

# some preliminaries
station <- 'norfolk'
output.dir <- '../output/'
dat.dir <- '../data/'

setwd('/Users/tony/codes/EVT/R')

today <- Sys.Date(); today=format(today,format="%d%b%Y")

appen <- paste('suppExp_temperatures_',station,today,sep='_')

filename.saveprogress <- paste('../output/suppExp_predictor_norfolk_',today,'.RData',sep='')


# data for calibration:
data_calib <- readRDS('../data/tidegauge_processed_norfolk_decl3-pot99-annual_06Dec2017.rds')

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





#===============================================================================
# read in temperature/time forcing data for the experiments
# * Experiments:
#     1) global annual mean temperature (control)
#     2) local annual mean temperature
#     3) time
#     4) global monthly mean temperature
#     5) local monthly mean temperature

# read in local monthly average temperatures; construct time series of annual
# means for Norfolk
months <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
files.norfolk <- list.files(path='../data/NOAA_temperature_data/local_temperature_norfolk/', pattern='csv')

wgts_norm <- (1/365)*c(31,28,31,30,31,30,31,31,30,31,30,31)
wgts_leap <- (1/366)*c(31,29,31,30,31,30,31,31,30,31,30,31)
wgts_2017 <- (1/(365-30-31))*c(31,28,31,30,31,30,31,31,30,31,0,0)

# read one just to get the appropriate size
filename <- paste('../data/NOAA_temperature_data/local_temperature_norfolk/',
                  files.norfolk[which( substr(files.norfolk,1,3)==months[1])],
                  sep='')
tmp <- read.table(filename, header = TRUE, sep=',', skip=4)
years <- tmp$Year
n_year <- length(years)
temp_monthly <- mat.or.vec(n_year, length(months))
temp_annual <- mat.or.vec(n_year, 1)
rownames(temp_monthly) <- years
colnames(temp_monthly) <- months

# note that for november and december, you'll have a 0 fill row at the end (2017)
# and averaging weights for that year will be different
for (mm in 1:length(months)) {
  filename <- paste('../data/NOAA_temperature_data/local_temperature_norfolk/',
                    files.norfolk[which( substr(files.norfolk,1,3)==months[mm])],
                    sep='')
  new_temps <- read.table(filename, header = TRUE, sep=',', skip=4)
  if(max(new_temps[,1])==2017) {
    temp_monthly[,mm] <- new_temps[,2]
  } else if(max(new_temps[,1])==2016) {
    temp_monthly[1:length(new_temps[,1]),mm] <- new_temps[,2]
  }
}

# Useful note:  if [year] %% 4 == 0, then you're in a leap year
for (yy in 1:n_year) {
  if (years[yy] %% 4 == 0) {
    # leap year
    temp_annual[yy] <- sum(temp_monthly[yy,] * wgts_leap)
  } else if(years[yy] < 2017) {
    # normal year
    temp_annual[yy] <- sum(temp_monthly[yy,] * wgts_norm)
  } else if(years[yy]==2017) {
    # 2017 year (missing november and december)
    temp_annual[yy] <- sum(temp_monthly[yy,] * wgts_2017)
  }
}

# put together the 3 alternative predictor data sets

names_predictor <- c('global_annual', 'local_annual', 'time')
year_auxiliary <- vector('list', 3); names(year_auxiliary) <- names_predictor
forc_auxiliary <- vector('list', 3); names(forc_auxiliary) <- names_predictor

# 1. global annual mean temperature

source('read_data_temperature.R')
year_auxiliary$global_annual <- years
forc_auxiliary$global_annual <- temperature_forc[which(time_forc==min(years)):which(time_forc==max(years))]

# 2. local annual mean temperature

year_auxiliary$local_annual <- years
forc_auxiliary$local_annual <- temp_annual

# 3. time  -- rescaled to same scale as the temperatures

year_auxiliary$time <- years
# scales 'time' between 0 and 1
forc_auxiliary$time <- (years - years[1])/(max(years)-years[1])
# scales 'time' between 0 and max(temp)-min(temp)
forc_auxiliary$time <- forc_auxiliary$time * (max(forc_auxiliary$global_annual)-min(forc_auxiliary$global_annual))

# All normalized to same time period? (1901-2000 mean)
year_norm <- c(1901,2000)
ind_norm <- which(years==year_norm[1]):which(years==year_norm[2])
for (pp in names_predictor) {
  forc_auxiliary[[pp]] <- forc_auxiliary[[pp]] - mean(forc_auxiliary[[pp]][ind_norm])
}

#===============================================================================





#===============================================================================
# calculate maximum likelihood parameter estimates for each experiment


# the likelihood function for PP/GPD is in here
source('likelihood_ppgpd.R')

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
  metrics[[mm]] <- mat.or.vec(length(names_predictor), nmodel)
  rownames(metrics[[mm]]) <- names_predictor
  colnames(metrics[[mm]]) <- types.of.gpd
}

# create object to hold the calibrated parameters and (negative) log-likelihood
deoptim.decl <- vector('list', nmodel); names(deoptim.decl) <- types.of.model
for (gpd.type in types.of.model) {
  deoptim.decl[[gpd.type]] <- mat.or.vec(length(names_predictor), length(parnames_all[[gpd.type]])+1)
  rownames(deoptim.decl[[gpd.type]]) <- names_predictor
  colnames(deoptim.decl[[gpd.type]]) <- c(parnames_all[[gpd.type]], 'nnlike')
}


for (predictor in names_predictor) {

  for (gpd.type in types.of.gpd) {
    print(paste(' - starting DE optimization for model ',gpd.type,'...',sep=''))
    if(gpd.type=='gpd3') {auxiliary <- NULL
    } else {auxiliary <- trimmed_forcing(data_calib$gpd$year, year_auxiliary[[predictor]], forc_auxiliary[[predictor]])$forcing}
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
    deoptim.decl[[gpd.type]][predictor,] <- c(out.deoptim$optim$bestmem, out.deoptim$optim$bestval)

    metrics$maxlik[predictor,gpd.type] <- -deoptim.decl[[gpd.type]][predictor,length(parnames_all[[gpd.type]])+1]
    metrics$bic[predictor,gpd.type] <- -2*metrics$maxlik[predictor,gpd.type] + length(parnames_all[[gpd.type]])*log(length(data_calib$gpd$counts_all))
  }
}
#===============================================================================


#===============================================================================
# Save results; can come back later
save.image(file=filename.saveprogress)
#===============================================================================


#===============================================================================
# write table with the metrics

likrat <- metrics$maxlik
for (i in 1:4) {
    likrat[,i] <- exp(likrat[,i]-likrat[1,i])
}

to_file <- rbind( cbind(rep('bic',nrow(metrics$bic)), round(metrics$bic,2)),
                  cbind(rep('maxlik',nrow(metrics$maxlik)), round(metrics$maxlik,2)),
                  cbind(rep('likrat',nrow(likrat)), round(likrat,2)))

write.csv(x=t(to_file), file='../output/predictors_norfolk_SOM.csv')
#===============================================================================


#want to get ratio lik1/lik2
#have llik1 and llik2
#llik1 - llik2 = log(lik1/lik2)
#<=>  exp(llik1 - llik2) = lik1/lik2


#
#===============================================================================
# End
#===============================================================================
#

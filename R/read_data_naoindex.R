#===============================================================================
# read_data_naoindex.R
#
# read NAO index forcing data
# historical from Jones et al 1997 (updated since then):
#
#
# future projections from... TODO
#
#
# Use Winter (DJFM) index.
#
# Questions? Tony Wong (twong@psu.edu)
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
# read historical NAO index data
nao_dat <- read.table('../data/nao_3dp.dat')
colnames(nao_dat) <- c('year','jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec','ann')
ibeg <- which(nao_dat['year']==1850)
iend <- max(which(nao_dat[,'ann']!=-99.99))
time_hist <- nao_dat[ibeg:iend, 'year']

# get DJFM means
nao_hist <- rep(-999, length(time_hist))
for (y in 1:length(time_hist)) {
  nao_hist[y] <- mean( c(nao_dat$dec[ibeg+y-1], nao_dat$jan[ibeg+y], nao_dat$feb[ibeg+y], nao_dat$mar[ibeg+y]) )
}
time_hist <- nao_dat[ibeg:iend, 'year']

# get DJFM means
nao_hist <- rep(-999, length(time_hist))
for (y in 1:length(time_hist)) {
  nao_hist[y] <- mean( c(nao_dat$dec[ibeg+y-1], nao_dat$jan[ibeg+y], nao_dat$feb[ibeg+y], nao_dat$mar[ibeg+y]) )
}
#===============================================================================


#===============================================================================
# get SLP for historical period (1958-2000) to normalize
slp_ice <- read.table('../data/nao_ice.dat')
slp_azo <- read.table('../data/nao_azo.dat', skip=1)
colnames(slp_ice) <- c('year','jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec')
colnames(slp_azo) <- c('year','jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec')

ibeg_ice <- which(slp_ice$year==1958)
iend_ice <- which(slp_ice$year==2000)
ibeg_azo <- which(slp_azo$year==1958)
iend_azo <- which(slp_azo$year==2000)

slp_ice_norm <- rep(-999, 12)
slp_azo_norm <- rep(-999, 12)
for (m in 1:12) {
  slp_ice_norm[m] <- mean(as.matrix(slp_ice)[ibeg_ice:iend_ice, 1+m])
  slp_azo_norm[m] <- mean(as.matrix(slp_azo)[ibeg_azo:iend_azo, 1+m])
}
#===============================================================================


#===============================================================================
# read projections NAO index
library(ncdf4)

file.in <- '../data/DMIEH5_SRA1B_4_MM_psl.1-1200.nc'
ncdata <- nc_open(file.in)
psl <- ncvar_get(ncdata, 'psl')
lon <- ncvar_get(ncdata, 'lon')
lat <- ncvar_get(ncdata, 'lat')
time <- ncvar_get(ncdata, 'time')  # hours after 2001-01-31 (2001-2100 data)
nc_close(ncdata)

n_month <- length(time)

# As in Li and Wang, 2003 (http://www.lasg.ac.cn/staff/ljp/paperE/ljp_2003NAO.pdf)
# use proxy for station-based, since there is disagreement over how exactly
# to do the PCA (which EOFs to rotate...) and loss of direct physical interpretation
# of the results
ilat_azores <- which.min(abs(lat - 35))
ilon_azores <- which(lon >=(360-80) | lon <= 30)

ilat_iceland <- which.min(abs(lat - 65))
ilon_iceland <- which(lon >=(360-80) | lon <= 30)

psl_azores <- apply(X=psl[ilon_azores, ilat_azores, ], MARGIN=2, FUN=mean)
psl_iceland <- apply(X=psl[ilon_iceland, ilat_iceland, ], MARGIN=2, FUN=mean)

# normalize relative to long-term mean
for (m in 1:dim(time)) {
#  mon <- m%%12
#  if (mon==0) {mon <- 12}
#  psl_azores[m] <- psl_azores[m] - slp_azo_norm[mon]*10
#  psl_iceland[m] <- psl_iceland[m] - slp_ice_norm[mon]*10
  psl_azores[m] <- psl_azores[m] - mean(psl[ilon_azores, ilat_azores, seq(from=m%%12, to=n_month, by=12)])
  psl_iceland[m] <- psl_iceland[m] - mean(psl[ilon_iceland, ilat_iceland, seq(from=m%%12, to=n_month, by=12)])
}

# standardize relative to standard deviations
psl_azores <- psl_azores / sd(psl_azores)
psl_iceland <- psl_iceland / sd(psl_iceland)

# take winter (DJFM mean)
nao_monthly <- psl_azores - psl_iceland
nao_proj <- rep(-999, 99)
for (y in 1:99) {
  nao_proj[y] <- mean(nao_monthly[(y-1)*12 + 12:15])
#  nao_proj[y] <- mean(nao_monthly[(y-1)*12 + 1:12])
}
time_proj <- 2001:2099
#===============================================================================


#===============================================================================
# now note that the historical and projections would have assigned to year 2001,
# for example, December 2001 and January/February/March of 2002. Both are
# consistent with this convention, so it works to link the two together as-is.

# when to begin using the projections of NAO index
tbeg <- max(time_hist) + 1
ibeg <- which(time_proj==tbeg)

time_forc <- c(time_hist, time_proj[ibeg:length(time_proj)])
nao_forc <- c(nao_hist, nao_proj[ibeg:length(time_proj)])

forc_max <- max(nao_forc)
#===============================================================================


#===============================================================================
# Useful:

# function to trim forcing field to fit TG record unique years
# there might be missing years in TG record, so need to match each year and not
# just plop down an evenly spaced sequence
trimmed_forcing <- function(year_tidegauge, year_forcing, forcing) {
  output <- vector('list', 2); names(output) <- c('time','forcing')
  # check the beginning
  if(year_forcing[1] > year_tidegauge[1]) {print('ERROR - tide gauge record starts before forcing; add support for this situation')}
  # check the end
  if(max(year_forcing) < max(year_tidegauge)) {print('ERROR - tide gauge record ends after forcing; add support for this situation')}
  # match the indices of year_tidegauge within year_forcing
  imatch <- match(year_tidegauge, year_forcing)
  output$time <- year_forcing[imatch]
  output$forcing <- forcing[imatch]
  return(output)
}

#===============================================================================
# End
#===============================================================================

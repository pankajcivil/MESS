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

# get DJF means
nao_hist <- rep(-999, length(time_hist))
for (y in 1:length(time_hist)) {
  nao_hist[y] <- mean( c(nao_dat$dec[ibeg+y-1], nao_dat$jan[ibeg+y], nao_dat$feb[ibeg+y]) )
}
time_hist <- nao_dat[ibeg:iend, 'year']

# re-normalize (mean and stdev) relative to 2001-2016 mean/stdev, so it is
# consistent with the projections
ind_norm <- which(time_hist==2001):which(time_hist==2016)
nao_hist <- (nao_hist - mean(nao_hist[ind_norm]))/sd(nao_hist[ind_norm])
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
n_year <- n_month/12 -1

# As in Stephenson et al, 2006 (doi:  10.1007/s00382-006-0140-x)
# use regional, since there is disagreement over how exactly to do the PCA
# (which EOFs to rotate, e.g.) and a physical interpretation is direct this way
ilat_azores <- which(lat >= 20 & lat <= 55)
ilon_azores <- which(lon >=(360-90) | lon <= 60)

ilat_iceland <- which(lat >= 55 & lat <= 90)
ilon_iceland <- which(lon >=(360-90) | lon <= 60)

psl_azores <- psl[ilon_azores, ilat_azores, ]
psl_iceland <- psl[ilon_iceland, ilat_iceland, ]

# Don't do area-weighting as per explanation from Jesse Nusbaumer 27 March 2018
# email. Jesse says:
#   Don't area-weight the results.  Instead, just take the average of the SLP
#   over the specified region, treating every grid box the same.  The reason is
#   because in GCM papers they will usually specifically state "area-weighted"
#   if the actual meters-squared area is being taken into account, otherwise it
#   is assumed to just be a non-weighted average.  Also when it comes to a lot
#   of these indices the actual physical values (e.g. the total amount of
#   atmospheric mass producing the surface pressure) isn't really important,
#   they are just using the values to produce a strong statistical relationship.

# take the bulk average over each area, for each month
psl_azores_spavg <- apply(psl_azores, 3, mean)
psl_iceland_spavg <- apply(psl_iceland, 3, mean)

# Standardize each site separately (as discussed in Jones et al 1997). Do relative
# to 2001-2016 mean/stdev so we can be consistent between projections and the
# historical record
psl_azores_std <- rep(-999, n_month)
psl_iceland_std <- rep(-999, n_month)
for (m in 1:12) {
  ind_this_month <- seq(from=m, to=n_month, by=12)
  # first 16 are 2001-2016
  psl_tmp <- psl_azores_spavg[ind_this_month]
  psl_azores_std[ind_this_month] <- (psl_tmp - mean(psl_tmp[1:16]))/sd(psl_tmp[1:16])
  psl_tmp <- psl_iceland_spavg[ind_this_month]
  psl_iceland_std[ind_this_month] <- (psl_tmp - mean(psl_tmp[1:16]))/sd(psl_tmp[1:16])
}

# get SLP difference
nao_monthly <- psl_azores_std - psl_iceland_std

# get winter mean
nao_proj <- rep(-999, n_year)
for (y in 1:n_year) {
  nao_proj[y] <- mean(nao_monthly[(y-1)*12 + 12:14])  # DJF
#  nao_proj[y] <- mean(nao_monthly[(y-1)*12 + 12:15])  # DJFM
#  nao_proj[y] <- mean(nao_monthly[(y-1)*12 + 1:12])   # annual
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

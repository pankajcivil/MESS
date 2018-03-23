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

# read historical NAO index data
nao_hist <- read.table('../data/nao_3dp.dat')
colnames(nao_hist) <- c('year','jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec','ann')
time_hist <- nao_hist[,'year']
ibeg <- which(time_hist==1850)
iend <- max(which(nao_hist[,'ann']!=-99.99))
time_hist <- time_hist[ibeg:iend]
nao_hist <- nao_hist[ibeg:iend, 'ann']


# read projections NAO index

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< todo here now todo
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< todo here now todo
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< todo here now todo

# modify once the projections are in there
#auxiliary_forc <- c(nao_hist, nao_proc)
#time_forc <- c(time_hist, time_proj)
nao_forc <- nao_hist
time_forc <- time_hist

forc_max <- max(nao_forc)


if(FALSE) {
# read the projections temperature forcing from CRNM (CMIP5)
# note: these are in K, but going to normalize, so will take a difference and
# it's same as celsius
ncdata <- nc_open('../data/global.tas.aann.CNRM-CM5.historical+rcp85.r1i1p1.18500101-21001231.nc')
   temperature_proj <- ncvar_get(ncdata, 'tas')
   time_proj <- ncvar_get(ncdata, 'time')
nc_close(ncdata)

# 'time' on the netcdf ile is YYYYMMDD, where MMDD is July 2 each year (becuase
# of the averaging). so pluck off just the years
time_proj <- floor(time_proj/10000)

# read the historical forcing from NOAA
data.tmp <- read.table('../data/noaa_temperature_1880-2017.csv', header = TRUE, sep=',')
time_hist <- data.tmp$Year
temperature_hist <- data.tmp$Value

# extend historical back to 1850 with HadCRUT4
data.tmp = read.table('../data/HadCRUT.4.4.0.0.annual_ns_avg.txt')
time_hadcrut = data.tmp[,1]
temperature_hadcrut = data.tmp[,2]

# normalize all to 1901-2000
ind_norm <- which(time_hist==1901):which(time_hist==2000)
temperature_hist <- temperature_hist - mean(temperature_hist[ind_norm])

ind_norm <- which(time_hadcrut==1901):which(time_hadcrut==2000)
temperature_hadcrut <- temperature_hadcrut - mean(temperature_hadcrut[ind_norm])

ind_norm <- which(time_proj==1901):which(time_proj==2000)
temperature_proj <- temperature_proj - mean(temperature_proj[ind_norm])

# set up the forcing, wherein the historical starts, then projections finish
# 1850-1880 -- hadcrut4
# 1880-2016 -- NOAA historical
# 2017-2100 -- CRNM projection
time_forc <- min(time_hadcrut):max(time_proj)
temperature_forc <- rep(NA, length(time_forc))

ind_hadcrut <- which(time_hadcrut==time_forc[1]):(which(time_hadcrut==time_hist[1])-1)
ind_hist    <- 1:length(time_hist)
ind_proj    <- which(time_proj==(max(time_hist)+1)):which(time_proj==max(time_forc))

temperature_forc <- c(temperature_hadcrut[ind_hadcrut],
                      temperature_hist[ind_hist]      ,
                      temperature_proj[ind_proj]      )

# maximum temperature serves as an additinoal prior constraint on kappa0, kappa1
# that is, kappa1 > -kappa0/Tmax (otherwise, kappa = kappa0 + kappa1*T might be
# negative)
Tmax <- max(temperature_forc)
}





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

if(FALSE){
# checking if the first 9 months are a decent estimate of the entire year's index
  tmp <- read.table('../data/nao_3dp.dat')
  est_9mo = NULL
  est_ann = NULL
  for (yi in 1:nrow(tmp)) {
    if(!any(tmp[yi,]==-99.99)) {
      s = 0; for (mi in 2:10) {s = s+tmp[yi,mi]}
      est_9mo = c(est_9mo, s/9)
      est_ann = c(est_ann, tmp[yi,14])
    }
  }
}

#===============================================================================
# read temperature data
# historical from NOAA:
#   NOAA National Centers for Environmental information, Climate at a Glance:
#   Global Time Series, published May 2017, retrieved on June 7, 2017 from
#   http://www.ncdc.noaa.gov/cag/
# future projections from CRNM under RCP8.5, part of CMIP5:
#
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================

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

# normalize both to 1901-2000
ind_norm <- which(time_hist==1901):which(time_hist==2000)
temperature_hist <- temperature_hist - mean(temperature_hist[ind_norm])

ind_norm <- which(time_proj==1901):which(time_proj==2000)
temperature_proj <- temperature_proj - mean(temperature_proj[ind_norm])

# set up the forcing, wherein the historical starts, then projections finish
# 1880-2016 -- historical
# 2017-2100 -- projection
time_forc <- min(time_hist):max(time_proj)
temperature_forc <- rep(NA, length(time_forc))
temperature_forc[1:length(time_hist)] <- temperature_hist
temperature_forc[(length(time_hist)+1):length(time_forc)] <- temperature_proj[(which(time_proj==max(time_hist))+1):length(temperature_proj)]

#===============================================================================
# End
#===============================================================================

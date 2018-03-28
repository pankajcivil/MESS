## nao_calculation_projections.R
##
## Data obtained from
##  https://cera-www.dkrz.de/WDCC/ui/cerasearch/downloadForm?acronym=DMIEH5_SRA1B_4_MM_psl
## on 26 March 2018
##
## Calculation of NAO index for each month (indexed by m) proceeds as in Jianping
## and Wang 2003 (doi: 10.1007/BF02915394)
##
##
## Questions? Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================

#install.packages('geosphere')

library(ncdf4)
library(geosphere)

file.in <- '../data/DMIEH5_SRA1B_4_MM_psl.1-1200.nc'

ncdata <- nc_open(file.in)
psl <- ncvar_get(ncdata, 'psl')
lon <- ncvar_get(ncdata, 'lon')
lat <- ncvar_get(ncdata, 'lat')
time <- ncvar_get(ncdata, 'time')
nc_close(ncdata)

n_month <- length(time)


#===============================================================================
# get SLP for historical period (1951-1980) to normalize, following Jones et al
# (1997) to be consistent with the historical data
slp_ice <- read.table('../data/nao_ice.dat')
slp_azo <- read.table('../data/nao_azo.dat', skip=1)
colnames(slp_ice) <- c('year','jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec')
colnames(slp_azo) <- c('year','jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec')

ibeg_ice <- which(slp_ice$year==1951)
iend_ice <- which(slp_ice$year==1980)
ibeg_azo <- which(slp_azo$year==1951)
iend_azo <- which(slp_azo$year==1980)

slp_ice_norm <- rep(-999, 12)
slp_azo_norm <- rep(-999, 12)
slp_ice_sd <- rep(-999, 12)
slp_azo_sd <- rep(-999, 12)
for (m in 1:12) {
  slp_ice_norm[m] <- mean(10*as.matrix(slp_ice)[ibeg_ice:iend_ice, 1+m])
  slp_azo_norm[m] <- mean(10*as.matrix(slp_azo)[ibeg_azo:iend_azo, 1+m])
  slp_ice_sd[m] <- sd(10*as.matrix(slp_ice)[ibeg_ice:iend_ice, 1+m])
  slp_azo_sd[m] <- sd(10*as.matrix(slp_azo)[ibeg_azo:iend_azo, 1+m])
}
#===============================================================================


##==============================================================================
if(TRUE){
# As in Stephenson et al, 2006 (doi:  10.1007/s00382-006-0140-x)
# use proxy for station-based, since there is disagreement over how exactly
# to do the PCA (which EOFs to rotate...) and interpretation is easier this way
ilat_azores <- which(lat >= 20 & lat <= 55)
ilon_azores <- which(lon >=(360-90) | lon <= 60)

ilat_iceland <- which(lat >= 55 & lat <= 90)
ilon_iceland <- which(lon >=(360-90) | lon <= 60)

psl_azores <- psl[ilon_azores, ilat_azores, ]
psl_iceland <- psl[ilon_iceland, ilat_iceland, ]

# take the bulk average over each area, for each month
psl_azores_mean <- apply(psl_azores, 3, mean)
psl_iceland_mean <- apply(psl_iceland, 3, mean)

# normalize each site separately (as discussed in Jones et al 1997). Do relative
# to 2001-2016 mean/stdev so we can be consistent between projections and the
# historical record
psl_azores_norm <- rep(-999, n_month)
psl_iceland_norm <- rep(-999, n_month)
for (m in 1:12) {
  ind_this_month <- seq(from=m, to=n_month, by=12)
  # first 16 are 2001-2016
  psl_tmp <- psl_azores_mean[ind_this_month]
  psl_azores_norm[ind_this_month] <- (psl_tmp - mean(psl_tmp[1:16]))/sd(psl_tmp[1:16])
  psl_tmp <- psl_iceland_mean[ind_this_month]
  psl_iceland_norm[ind_this_month] <- (psl_tmp - mean(psl_tmp[1:16]))/sd(psl_tmp[1:16])
}

# get SLP difference
nao_monthly <- psl_azores_norm - psl_iceland_norm

# get winter mean
nao_proj <- rep(-999, 99)
for (y in 1:99) {
  nao_proj[y] <- mean(nao_monthly[(y-1)*12 + 12:14])  # DJF
#  nao_proj[y] <- mean(nao_monthly[(y-1)*12 + 12:15])  # DJFM
#  nao_proj[y] <- mean(nao_monthly[(y-1)*12 + 1:12])   # annual
}

}


##==============================================================================
if(FALSE){
# As in Li and Wang, 2003 (http://www.lasg.ac.cn/staff/ljp/paperE/ljp_2003NAO.pdf)
# use proxy for station-based, since htere is disagreement over how exactly
# to do the PCA (which EOFs to rotate...)
ilat_azores <- which.min(abs(lat - 35))
ilon_azores <- which(lon >=(360-80) | lon <= 30)

ilat_iceland <- which.min(abs(lat - 65))
ilon_iceland <- which(lon >=(360-80) | lon <= 30)

psl_azores <- apply(X=psl[ilon_azores, ilat_azores, ], MARGIN=2, FUN=mean)
psl_iceland <- apply(X=psl[ilon_iceland, ilat_iceland, ], MARGIN=2, FUN=mean)

# normalize relative to long-term mean
for (m in 1:dim(time)) {
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
}


##==============================================================================
if(FALSE){
# use proxy for station-based, since htere is disagreement over how exactly
# to do the PCA (which EOFs to rotate...)
ilat_azores <- which.min(abs(lat - 37.74))
ilon_azores <- which.min(abs(lon - (360-25.68)))

ilat_iceland <- which.min(abs(lat - 65.07))
ilon_iceland <- which.min(abs(lon - (360-22.74)))

psl_azores <- psl[ilon_azores, ilat_azores, ]
psl_iceland <- psl[ilon_iceland, ilat_iceland, ]

# normalize relative to long-term mean
for (m in 1:dim(time)) {
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
}
}

##==============================================================================
if (FALSE) {
# calculate based on atlantic region first principle component...
# azores:
ilat_azores <- which(lat>=25 & lat<=55)
ilon_azores <- which(lon >= 310 | lon <= 10)

# iceland:
ilat_iceland <- which(lat>=55 & lat<=85)
ilon_iceland <- which(lon >= 320 | lon <= 20)

# atlantic (PCA/EOF-based)
ilat_atlantic <- which(lat>=20 & lat<=70)
ilon_atlantic <- which(lon>=270 | lon<=40)

# subtract off long-term mean for each month
psl_norm <- psl
for (ilon in 1:dim(lon)) {
  for (ilat in 1:dim(lat)) {
    for (m in 1:dim(time)) {
      psl_norm[ilon, ilat, m] <- psl[ilon, ilat, m] - mean(psl[ilon, ilat, seq(from=m%%12, to=n_month, by=12)])
    }
  }
}

# only looking in these two domains
psl_azores <- psl_norm[ilon_azores, ilat_azores, ]
psl_iceland <- psl_norm[ilon_iceland, ilat_iceland, ]
psl_atlantic <- psl_norm[ilon_atlantic, ilat_atlantic, ]

# scratch?

# calculate seasonal means
tmp <- array(rep(-999,length(ilat_atlantic)*length(ilon_atlantic)*100), c(length(ilon_atlantic), length(ilat_atlantic), 100))
for (ilon in 1:length(ilon_atlantic)) {
  for (ilat in 1:length(ilat_atlantic)) {
    for (y in 1:99) {
      tmp[ilon, ilat, y] <- mean(psl_atlantic[ilon, ilat, (y-1)*12+12:15])
    }
  }
}
}

##==============================================================================

if(FALSE) { # don't do area-weighting as per explanation from Jesse Nusbaumer 27 March 2018 email
# area-weighting
# everything at same latitude should have same weight

deltaLat <- abs(median(diff(lat)))
deltaLon <- abs(median(diff(lon)))

lat_azores <- lat[ilat_azores]
wgt_azores <- rep(-999, length(lat_azores))
lon_tmp <- lon[ilon_azores[1]]
R_earth <- 6378137 # meters

# assumes the lat/lon given are the center of the gridcell
for (l in 1:length(wgt_azores)) {
  corners <- rbind( c(lon_tmp-0.5*deltaLon, lat_azores[l]-0.5*deltaLat),
                    c(lon_tmp-0.5*deltaLon, lat_azores[l]+0.5*deltaLat),
                    c(lon_tmp+0.5*deltaLon, lat_azores[l]+0.5*deltaLat),
                    c(lon_tmp+0.5*deltaLon, lat_azores[l]-0.5*deltaLat))
  wgt_azores[l] <- areaPolygon(corners, a=R_earth)
}

lat_iceland <- lat[ilat_iceland]
wgt_iceland <- rep(-999, length(lat_iceland))
lon_tmp <- lon[ilon_iceland[1]]
R_earth <- 6378137 # meters

# assumes the lat/lon given are the center of the gridcell
for (l in 1:length(wgt_iceland)) {
  corners <- rbind( c(lon_tmp-0.5*deltaLon, lat_iceland[l]-0.5*deltaLat),
                    c(lon_tmp-0.5*deltaLon, lat_iceland[l]+0.5*deltaLat),
                    c(lon_tmp+0.5*deltaLon, lat_iceland[l]+0.5*deltaLat),
                    c(lon_tmp+0.5*deltaLon, lat_iceland[l]-0.5*deltaLat))
  wgt_iceland[l] <- areaPolygon(corners, a=R_earth)
}

# normalize so all lons within a lat band sum to 1
wgt_azores <- wgt_azores/sum(wgt_azores)
wgt_iceland <- wgt_iceland/sum(wgt_iceland)
# normalize again because there are length(ilon_azores/iceland) lon bands
wgt_azores <- wgt_azores/length(ilon_azores)
wgt_iceland <- wgt_iceland/length(ilon_iceland)

# take the area-weighted mean
psl_azores_mean <- rep(0, n_month)
psl_iceland_mean <- rep(0, n_month)
for (t in 1:n_month) {
  for (ilon in 1:length(ilon_azores)) {
    psl_azores_mean[t] <- psl_azores_mean[t] + sum(wgt_azores*psl_azores[ilon, , t])
  }
  for (ilon in 1:length(ilon_iceland)) {
    psl_iceland_mean[t] <- psl_iceland_mean[t] + sum(wgt_iceland*psl_iceland[ilon, , t])
  }
}
}


if (FALSE) {
## Check against the historical NAO index data
# read historical NAO index data
nao_dat <- read.table('../data/nao_3dp.dat')
colnames(nao_dat) <- c('year','jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec','ann')
ibeg <- which(nao_dat['year']==1850)
iend <- max(which(nao_dat[,'ann']!=-99.99))
time_hist <- nao_dat[ibeg:iend, 'year']

# get winter means
nao_hist <- rep(-999, length(time_hist))
for (y in 1:length(time_hist)) {
  nao_hist[y] <- mean( c(nao_dat$dec[ibeg+y-1], nao_dat$jan[ibeg+y], nao_dat$feb[ibeg+y]) )  # DJF
#  nao_hist[y] <- mean( c(nao_dat$dec[ibeg+y-1], nao_dat$jan[ibeg+y], nao_dat$feb[ibeg+y], nao_dat$mar[ibeg+y]) )  # DJFM
}
time_hist <- nao_dat[ibeg:iend, 'year']

# re-normalize (mean and stdev) relative to 2001-2016 mean/stdev, so it is
# consistent with the projections
ind_norm <- which(time_hist==2001):which(time_hist==2016)
nao_hist <- (nao_hist - mean(nao_hist[ind_norm]))/sd(nao_hist[ind_norm])

icomp <- which(time_hist==2001):which(time_hist==2016)
}

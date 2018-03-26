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

library(ncdf4)

file.in <- '../data/DMIEH5_SRA1B_4_MM_psl.1-1200.nc'

ncdata <- nc_open(file.in)
psl <- ncvar_get(ncdata, 'psl')
lon <- ncvar_get(ncdata, 'lon')
lat <- ncvar_get(ncdata, 'lat')
time <- ncvar_get(ncdata, 'time')
nc_close(ncdata)

##==============================================================================
if(TRUE){
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
n_months <- length(time)
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

if (FALSE) {
## Check against the historical NAO index data
nao_hist <- read.table('../data/nao_3dp.dat')
colnames(nao_hist) <- c('year','jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec','ann')
ibeg <- which(time_hist==1850)
iend <- max(which(nao_hist[,'ann']!=-99.99))
time_hist <- nao_hist[ibeg:iend, 'year']

# get DJFM means
nao <- rep(-999, length(time_hist))
for (y in 1:length(time_hist)) {
  nao[y] <- mean( c(nao_hist$dec[ibeg+y-1], nao_hist$jan[ibeg+y], nao_hist$feb[ibeg+y], nao_hist$mar[ibeg+y]) )
}
nao_hist <- nao
}

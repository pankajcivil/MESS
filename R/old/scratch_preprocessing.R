#===============================================================================
# preprocessing of tide gauge data to fit GEV, GPD, Naveau, etc. models at
# sub-annual block maxima. need to achieve independence among samples.
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================

#rm(list=ls())
NP.deoptim000 <- 200      # number of DE population members (at least 10*[# parameters])
niter.deoptim000 <- 200   # number of DE iterations
#n_node000 <- 1            # number of CPUs to use
output.dir <- '../output/'
#setwd('/storage/home/axw322/work/codes/EVT/R')
setwd('/Users/axw322/codes/EVT/R')
appen <- 'gev_nav'

#
#===============================================================================
# relevant libraries - do 'install.pacakges([library name])' if you do not have
# one yet
#===============================================================================
#

library(extRemes)
library(zoo)
library(adaptMCMC)
library(lhs)
library(DEoptim)
library(ncdf4)
library(date)

#
#===============================================================================
# read and process data for temperature (auxiliary covariate for nonstationary
# parameters)
# yields: temperature_forc, time_forc
#===============================================================================
#

print('reading temperature data...')

source('read_data_temperature.R')

print('...done.')

#
#===============================================================================
# read and process data for Dutch station (Delfzijl)
#===============================================================================
#

print('reading Delfzijl, Netherlands, tide gauge data...')

source('read_data_tidegauge_delfzijl.R')

print('...done.')

#
#===============================================================================
# read and process data for set of European stations
#===============================================================================
#

print('reading lots of other European tide gauge data...')

filetype='csv'
dat.dir <- '~/codes/EVT/data/tide_gauge_Europe/'
files.tg <- list.files(path=dat.dir, pattern=filetype)

data_set <- vector('list', length(files.tg))
for (dd in 1:length(files.tg)) {
  names(data_set)[dd] <- substr(files.tg[dd], start=1, stop=7)
  data.tmp <- read.table(paste(dat.dir,files.tg[dd],sep=''), header=FALSE, sep=',')
  data_set[[dd]] <- vector('list', 5)
  names(data_set[[dd]]) <- c('year','month','day','hour','sl')
  data_set[[dd]]$year <- data.tmp$V1
  data_set[[dd]]$month <- data.tmp$V2
  data_set[[dd]]$day <- data.tmp$V3
  data_set[[dd]]$hour <- data.tmp$V4
  data_set[[dd]]$sl <- data.tmp$V5
}

# process data

# annual block maxima
for (dd in 1:length(data_set)) {
  data_set[[dd]]$sl_norm     <- rep(NA, length(data_set[[dd]]$year))
  data_set[[dd]]$year_unique <- unique(data_set[[dd]]$year)
  data_set[[dd]]$year_unique <- data_set[[dd]]$year_unique[order(data_set[[dd]]$year_unique)]
  data_set[[dd]]$lsl_max     <- rep(NA, length(data_set[[dd]]$year_unique))
  for (t in 1:length(data_set[[dd]]$year_unique)) {
    ind_this_year <- which(data_set[[dd]]$year==data_set[[dd]]$year_unique[t])
    data_set[[dd]]$sl_norm[ind_this_year] <- data_set[[dd]]$sl[ind_this_year] - mean(data_set[[dd]]$sl[ind_this_year])
    data_set[[dd]]$lsl_max[t] <- max(data_set[[dd]]$sl_norm[ind_this_year])
  }
}

# trim before 1850 (or whenever is time_forc[1]), which is when the temperature
# forcing starts. also make a note of how long each record is
nyear <- rep(NA, length(data_set))
for (dd in 1:length(data_set)) {
  if(data_set[[dd]]$year_unique[1] < time_forc[1]) {
    icut <- which(data_set[[dd]]$year_unique < time_forc[1])
    data_set[[dd]]$year_unique <- data_set[[dd]]$year_unique[-icut]
    data_set[[dd]]$lsl_max <- data_set[[dd]]$lsl_max[-icut]
  }
  nyear[dd] <- length(data_set[[dd]]$year_unique)
}

print('...done.')

#
#===============================================================================
# CAUTION - UNDER CONSTRUCTION - CAUTION - UNDER CONSTRUCTION - CAUTION
#===============================================================================
#



#===
#=== convert Delfzijl tide gauge data to hourly series
#===

dat.dir <- '~/codes/EVT/data/tide_gauge_Europe/Delfzijl_Oddo_data/'
data <- read_data(dat.dir=dat.dir, filetype='txt', septype='\t')

# convert sea levels from cm to mm, consistent with the other European data
data$sl <- 10* data$sl

data$year   <- as.numeric(substr(as.character(data$date), start=1, stop=4))
data$month  <- as.numeric(substr(as.character(data$date), start=5, stop=6))
data$day    <- as.numeric(substr(as.character(data$date), start=7, stop=8))
data$hour   <- as.numeric(substr(data$time, 1,2))
data$minute <- as.numeric(substr(data$time, 4,5))

# time in days since 01 January 1960
data$time.days <- as.numeric(mdy.date(month=month, day=day, year=year)) + hour/24 + minute/(24*60)

# difference between time stamps (in units of days)
time.diff <- diff(time.days)

# check that they are in the proper order, ascending
print(paste('Are there any times out of order? ',any(time.diff < 0), sep=''))

# put everything back in order - make sure you do this for the sea levels and
# other fields too, so do as a matrix. also recalculate the time differences,
# which you will need for averaging
data <- data[order(data$time.days),]
time.diff <- diff(data$time.days)

# what is three hours? in units of days
three.hours <- 3/24

# where are there gaps longer than three hours? (+10sec for precision)
igap <- which(time.diff > (three.hours+10/(24*60*60)))

# where are there gaps shorter than three hours? (-10sec)
# should find lots of places. need to average these up to 3 hours.
iokay <- which( (time.diff < (three.hours+10/(24*60*60))) & (time.diff > (three.hours-10/(24*60*60))) )
ishort <- which(time.diff < (three.hours-10/(24*60*60)))

i10min <- which( (time.diff < (10/60/24+10/(24*60*60))) & (time.diff > (10/60/24-10/(24*60*60))) )
i1hour <- which( (time.diff < (1/24+10/(24*60*60))) & (time.diff > (1/24-10/(24*60*60))) )
i3hour <- which( (time.diff < (3/24+10/(24*60*60))) & (time.diff > (3/24-10/(24*60*60))) )

tnew.1hour <- seq(from=data$time.days[max(i3hour)]+three.hours, to=data$time.days[max(i1hour)], by=three.hours)
tnew.10min <- seq(from=data$time.days[max(i1hour)]+three.hours, to=data$time.days[max(i10min)], by=three.hours)
tnew.3hour <- data$time.days[i3hour]

# there is one missing 3-hour data point, between indices 239599 and 239600
# the 3-hourly data get averaged up to daily maxima time series; just ignore
# this and any other gaps, since the largest one is 8 hours in November 2014.
# still have enough data to get a meaningful daily maximum.

time.days.3hour <- c(tnew.3hour, tnew.1hour, tnew.10min)

# average up to three hourly time series

sl.3hour <- rep(NA, length(time.days.3hour))
year.3hour <- rep(NA, length(time.days.3hour))

# plug in the data that are already 3-hourly
sl.3hour[1:length(tnew.3hour)] <- data$sl[i3hour]
year.3hour[1:length(tnew.3hour)] <- data$year[i3hour]

# average up the data that are 1-hourly
print('Averaging hourly measurements up to 3-hourly...')
pb <- txtProgressBar(min=0,max=length(tnew.1hour),initial=0,style=3)
sl.tmp <- rep(NA, length(tnew.1hour))
year.tmp <- rep(NA, length(tnew.1hour))
for (t in 1:length(tnew.1hour)) {
  itmp <- which( (data$time.days > (tnew.1hour[t]-three.hours)) & (data$time.days <= tnew.1hour[t]) )
  sl.tmp[t] <- mean(data$sl[itmp], na.rm=TRUE)
  year.tmp[t] <- median(data$year[itmp], na.rm=TRUE)
  setTxtProgressBar(pb, t)
}
close(pb)
print('   ... done.')
sl.3hour[(length(tnew.3hour)+1):(length(tnew.3hour)+length(tnew.1hour))] <- sl.tmp
year.3hour[(length(tnew.3hour)+1):(length(tnew.3hour)+length(tnew.1hour))] <- year.tmp

# average up the data that are 10-minutely
print('Averaging 10-minutely measurements up to 3-hourly...')
pb <- txtProgressBar(min=0,max=length(tnew.10min),initial=0,style=3)
sl.tmp <- rep(NA, length(tnew.10min))
year.tmp <- rep(NA, length(tnew.10min))
for (t in 1:length(tnew.10min)) {
  itmp <- which( (data$time.days > (tnew.10min[t]-three.hours)) & (data$time.days <= tnew.10min[t]) )
  sl.tmp[t] <- mean(data$sl[itmp], na.rm=TRUE)
  year.tmp[t] <- median(data$year[itmp], na.rm=TRUE)
  setTxtProgressBar(pb, t)
}
close(pb)
print('   ... done.')
sl.3hour[(length(tnew.3hour)+length(tnew.1hour)+1):length(sl.3hour)] <- sl.tmp
year.3hour[(length(tnew.3hour)+length(tnew.1hour)+1):length(sl.3hour)] <- year.tmp

# that takes a long time, so save the workspace image
save.image(file='../output/preprocessing.RData')


#===
#=== subtract linear sea-level trend (from fit to monthly means)
#===

# calcualte monthly means

dates.new <- date.mdy(time.days.3hour)
date.beg <- date.mdy(min(time.days.3hour))
date.end <- date.mdy(max(time.days.3hour))

# how many months are there?
months.firstyear <- 12-date.beg$month+1
if(date.end$day > 15) {months.lastyear  <- date.end$month
} else {months.lastyear <- date.end$month-1}
# the -1 makes sure we don't count the first or last years
n.months <- 12*(date.end$year - date.beg$year - 1) + months.firstyear + months.lastyear

sl.month <- rep(NA, length(n.months))
time.month <- rep(NA, length(n.months))
cnt <- 1
print('getting time series of monthly means, to calculate sea-level rise trend...')
pb <- txtProgressBar(min=0,max=n.months,initial=0,style=3)
for (yy in date.beg$year:date.end$year) {
  ind.this.year <- which(dates.new$year == yy)
  for (mm in min(dates.new$month[ind.this.year]):max(dates.new$month[ind.this.year])) {
    ind.this.month <- which( (dates.new$year == yy) & (dates.new$month == mm) )
    sl.month[cnt] <- mean(sl.3hour[ind.this.month], na.rm=TRUE)
    time.month[cnt] <- mean(time.days.3hour[ind.this.month], na.rm=TRUE)
    cnt <- cnt+1
    setTxtProgressBar(pb, cnt)
  }
}
close(pb)

# fit trend
slr.trend <- lm(sl.month ~ time.month)
slr.trend.3hour <- slr.trend$coefficients[1] + slr.trend$coefficients[2]*time.days.3hour

# subtract off from the 3-hourly data
sl.3hour.detrended <- sl.3hour - slr.trend.3hour
print('  ... done.')


#===
#=== daily block maxima; calculate 99% quantile as GPD threshold
#===

n.days <- ceiling(max(time.days.3hour) - min(time.days.3hour))
sl.daily.max <- rep(NA, n.days)
year.daily.max <- rep(NA, n.days) # save which year each came from to organize into blocks
time.daily <- seq(from=0.5+floor(min(time.days.3hour)), to=0.5+floor(max(time.days.3hour))
cnt <- 1
print('getting time series of daily maxima...')
pb <- txtProgressBar(min=0,max=n.days,initial=0,style=3)
for (ttmp in time.daily) {
  ind.today <- which( abs(time.days.3hour-ttmp) < 0.5 )
  sl.daily.max[cnt] <- max(sl.3hour[ind.today], na.rm=TRUE)
  year.daily.max[cnt] <- median(year.3hour[ind.today], na.rm=TRUE)
  setTxtProgressBar(pb, cnt)
  cnt <- cnt+1
}
close(pb)
print('  ... done.')

# that takes a long time, so save the workspace image
save.image(file='../output/preprocessing.RData')


#===
#=== find all the excesses, "declustering" = if two are within a day of each
#=== other, take only the maximum of the two (so make sure you save the times
#=== of each excess)
#===

# threshold is 99% quantile of tide gauge's observed values.
# Buchanan et al (2016) use the 99% quantile of the daily maximum time series.
#gpd.threshold <- as.numeric(quantile(sl.3hour.detrended, .99, na.rm=TRUE))
gpd.threshold <- as.numeric(quantile(sl.daily.max, .99, na.rm=TRUE))

ind.exceed <- which(sl.daily.max >= gpd.threshold)
time.exceed <- time.daily[ind.exceed]
sl.exceed <- sl.daily.max[ind.exceed]
year.exceed <- year.daily.max[ind.exceed]

# whoops, extRemes package has a function for this, and it works better for
# the EVT analysis. ah well.
decluster_timeseries <- function(time, year, time.series, min.dt) {
  decluster <- vector('list',3)
  names(decluster) <- c('time','year','time.series')
  tdiff <- diff(time)
  ind.too.close <- which(tdiff <= min.dt)
  if(length(ind.too.close) > 0) {
    # go through the places where there are time.series values too close together,
    # and set to throw the smaller value away.
    # need to account for the possibility that there are more than two values in
    # a 'cluster'.
    # indices where new clusters begin: (tack on first one manually)
    ind.clusters <- c(ind.too.close[1] , ind.too.close[which(diff(ind.too.close) > 1) + 1])
    if(length(ind.clusters) > 1) {
      # case where there are multiple clusters to consider
      # initialize by fixing to remove all of the spots that are too close. then
      # we will remove the indices of each cluster's maximum
      irem <- unique(c(ind.too.close, ind.too.close + 1))
      for (i in 1:length(ind.clusters)) {
        # how long is this cluster?
        if(i < length(ind.clusters)) {
          first <- ind.clusters[i]
          last  <- ind.too.close[which(ind.too.close==ind.clusters[i+1]) -1]+1
        } else {
          first <- ind.clusters[i]
          if (length(which(ind.too.close > first))==0) {last <- first + 1
          } else {last <- first+length(which(ind.too.close > first)) + 1}
        }
        ind.this.cluster <- first:last
        isave <- ind.this.cluster[which.max(time.series[first:last])]
        # remove this cluster's maximum from irem; they might be out of order
        isave <- which(irem==isave)
        irem <- irem[-isave]
      }
    } else {
      # case with only one cluster
      # initialize by fixing to remove all of the cluster. then we will remove
      # the index of the cluster maximum from this
      irem <- unique(c(ind.too.close, ind.too.close + 1))
      irem <- irem[-which.max(time.series[irem])]
    }
    decluster$time <- time[-irem]
    decluster$year <- year[-irem]
    decluster$time.series <- time.series[-irem]
  } else {
    decluster$time <- time
    decluster$year <- year
    decluster$time.series <- time.series
  }
  return(decluster)
}

declustered.exceed <- decluster_timeseries(time=time.exceed, year=year.exceed, time.series=sl.exceed, min.dt=1)
time.exceed.decl <- declustered.exceed$time
year.exceed.decl <- declustered.exceed$year
sl.exceed.decl <- declustered.exceed$time.series


# TODO - decluster and then reinsert
sl.daily.max.declustered <- decluster(sl.daily.max, threshold=gpd.threshold, r=1 )
data.tmp <- data.frame( cbind(1:n.days, time.daily, sl.daily.max.declustered))
data.tmp$time.daily <- data.tmp$time.daily/365.25
names(data.tmp) <- c('obs','year','sl')

# try initial MLE fit for GP-PP model with extRemes package
fit0 <- fevd(sl, data.tmp, threshold=gpd.threshold, type="PP", units="mm", time.units='days')

# note that 'rate' is giving you estimate of
# [total number of exceedances] / [total time span]

# also note that devd does not have support for type='PP', so need to code the
# likelihood functions yourself. but let's be honest, you'd do that anyway.

# and finally, note that this MLE fit is pretty terrible! so code your own
# likelihood functions and MLE fit and see if you and old friend DE optimization
# can do better. probably some weird assumption in the model...

# non-stationary GP-PP model: lambda, scale, shape

data.gpd <- data.frame(cbind(1:length(sl.exceed.decl), time.exceed.decl/365.25, sl.exceed.decl))
names(data.gpd) <- c('obs','year_exceedance','sl_exceedance')

# note that data.gpd$year_exceedance has units of years after 1 Jan 1960

#===
#=== sub-list object on 'data_calib' to store what is needed to calbirate PP-GPD
#===

data_calib$gpd <- vector('list', 4)
names(data_calib$gpd) <- c('counts','time_length','excesses','threshold')
# initialize
data_calib$gpd$threshold <- gpd.threshold
data_calib$gpd$counts <- data_calib$gpd$time_length <- rep(NA, length(year.unique))
data_calib$gpd$excesses <- vector('list', length(year.unique))
for (y in 1:length(year.unique)) {
  hits.this.year <- which(year.exceed.decl==year.unique[y])
  data_calib$gpd$counts[y] <- length(hits.this.year)
    data_calib$gpd$time_length[y] <- max(time.daily[which(year.daily.max==year.unique[y])]) -
                                     min(time.daily[which(year.daily.max==year.unique[y])]) + 1
  if(length(hits.this.year) > 0) {data_calib$gpd$excesses[[y]] <- sl.exceed.decl[hits.this.year]
  } else                         {data_calib$gpd$excesses[[y]] <- NA}
}

# alternatively, could bin em all together. but this won't allow for potential
# non-stationary behavior in the poisson process
data_calib$gpd$excesses_all <- sl.exceed.decl
data_calib$gpd$counts_all <- length(sl.exceed.decl)
data_calib$gpd$time_length_all <- max(time.daily) - min(time.daily) + 1


## this works!
out.deoptim <- DEoptim(neg_log_like_ppgpd, lower=c(0,0,-2), upper=c(1,15000,2),
                    DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                    parnames=c('lambda','sigma','xi'), data_calib=data_calib, auxiliary=NULL)











##
## same as above but for arbitrary European station, and with preprocessing
##

# pick a single data set to work on the pre-processing. using something already
# an hourly time series, and long record
itmp <- which.min( abs(nyear-length(data_calib$lsl_max)))
data.tg <- data_set[[itmp]]

# convert to days since 01 January 1960
time.days <- as.numeric(mdy.date(month=data.tg$month, day=data.tg$day, year=data.tg$year)) + data.tg$hour/24

# difference between time stamps (in units of days)
time.diff <- diff(time.days)

# check that they are in the proper order, ascending
print(paste('Are there any times out of order? ',any(time.diff < 0), sep=''))

# what is an hour? in units of days
one.hour <- 1/24

# where are there gaps longer than 1 hour? (+/- 1 second, for numerical precision)
iokay <- which( (time.diff < (one.hour+1/(24*60*60))) & (time.diff > (one.hour-1/(24*60*60))) )
igap  <- 1:length(time.diff); igap <- igap[-iokay]


# TODO
# TODO -- GAP FILLING?
# TODO


# get series of daily averages
days <- floor(time.days)
days.unique <- unique(days)
lsl.daily <- rep(NA, length(days.unique))
for (d in 1:length(days.unique)) {
  ind.today <- which(days == days.unique[d])
  lsl.daily[d] <- mean(data.tg$sl[ind.today])
}

# for now, focus only on the longest stretch with no gaps
igap2 <- igap2 <- which(diff(days.unique) > 1)
igap.max <- which.max(diff(igap2))
lsl.daily.tmp <- lsl.daily[(igap2[igap.max]+1):igap2[igap.max+1]]
days.tmp <- days.unique[(igap2[igap.max]+1):igap2[igap.max+1]]

# pre-processing according to Ben and Beckett:

# 1. Run a Fast Fourier Transform (spectral decomposition) on your time series and then plot the periodogram
lsl.daily.fft <- abs(fft(lsl.daily)/sqrt(length(lsl.daily)))^2

f.data <- GeneCycle::periodogram(lsl.daily.tmp)
harmonics <- 1:4500
plot(f.data$freq[harmonics]*length(lsl.daily.tmp),
     f.data$spec[harmonics]/sum(f.data$spec),
     xlab="Harmonics", ylab="Amplitute Density", type="h")

amplitudes <- f.data$spec[harmonics]/sum(f.data$spec)
isort <- rev(order(amplitudes))


# alternatively...

# need to do a detrending
trend <- lm(lsl.daily.tmp ~ days.tmp)
lsl.daily.tmp.detrended <- trend$residuals

# this one won't look like anything:
fft.tmp = abs(fft(lsl.daily.tmp)/sqrt(length(lsl.daily.tmp)))^2
P = (4/length(lsl.daily.tmp))*fft.tmp[1:1100] # Only need the first (n/2)+1 values of the FFT result.
f = 1:1100 # this creates harmonic frequencies from 0 to .5 in steps of 1/128.
plot(f, P, type="l") # This plots the periodogram;

# this one (detrended) will
fft.tmp = abs(fft(lsl.daily.tmp.detrended)/sqrt(length(lsl.daily.tmp)))^2
P = (4/length(lsl.daily.tmp))*fft.tmp[1:1100] # Only need the first (n/2)+1 values of the FFT result.
f = 1:1100 # this creates harmonic frequencies from 0 to .5 in steps of 1/128.
plot(f, P, type="l") # This plots the periodogram;

isort <- rev(order(P))

# TODO
# TODO - in progress here now
# TODO


# 2. Select the top frequencies using the periodogram.
imodes <- which(P > P[isort[1]]*0.1)
imodes <- as.numeric(imodes[rev(order(P[imodes]))])

# 3. Run a Harmonic Regression where the predictors are sines and cosines of the chosen frequencies
hreg.tmp <- harmonic.regression(inputts=lsl.daily.tmp.detrended, inputtime=days.tmp, Tau=f[imodes[1]])

plot(days.tmp , lsl.daily.tmp.detrended, type='l')
lines(days.tmp, tmp$fit.vals, col='blue')


# playing around
play.data <- 60*cos(2*pi*days.tmp/5000) + 100*sin( 2*pi*days.tmp / 5000) + rnorm(n=length(days.tmp), mean=0, sd=50)
play.data2 <- 60*cos(2*pi*days.tmp/5000) + 100*sin( 2*pi*days.tmp / 5000) + rnorm(n=length(days.tmp), mean=0, sd=50) +
              80*cos(2*pi*days.tmp/500) + 50*sin( 2*pi*days.tmp / 500)
hreg.tmp <- harmonic.regression(inputts=play.data, inputtime=days.tmp, Tau=5000)

# this works! the harmonic analysis is just a linear model with the appropriate
# sin/cos predictors, at the modes of interest (imodes)
pred1c <- cos(2*pi*days.tmp/5000); pred1s <- sin(2*pi*days.tmp/5000); pred2c <- cos(2*pi*days.tmp/500); pred2s <- sin(2*pi*days.tmp/500);
fit <- lm(play.data2 ~ pred1c + pred1s + pred2c + pred2s)

acf1 <- acf(lsl.daily.tmp.detrended, plot=FALSE)

# using it for real:
predictors <- cbind(cos(2*pi*days.tmp/f[imodes[1]]), sin(2*pi*days.tmp/f[imodes[1]]))
if(length(imodes)>1) {
  for (i in 2:length(imodes)) {
    predictors <- cbind(predictors, cbind(cos(2*pi*days.tmp/f[imodes[i]]), sin(2*pi*days.tmp/f[imodes[i]])))
  }
}

harmonic.fit <- lm(lsl.daily.tmp.detrended ~ predictors)
#harmonic.fit$residuals
acf2 <- acf(harmonic.fit$residuals, plot=FALSE)

# 4. After the harmonic regression, we are left with the residuals. Fit an ARMA (2,2) model to the residuals.
arima.fit <- arima(x=harmonic.fit$residuals, order=c(2,0,2), include.mean=FALSE)
acf3 <- acf(arima.fit$residuals, plot=FALSE)

# 5. After fitting an ARMA (2,2) model to the residuals, we obtain the residuals for this second model.
#arima.fit$residuals

# simulate from the fitted ARMA process:
#tmp <- arima.sim(n=length(days.tmp), list(ar=c(arima.fit$coef['ar1'], arima.fit$coef['ar2']),
                  ma=c(arima.fit$coef['ma1'], arima.fit$coef['ma2'])), sd = sqrt(arima.fit$sigma2))

# 6. The 2nd set of residuals is the data that we'll be using for our GEV model fitting.


# 7. Check the ACF and PACF plots to ensure that there is little to no autocorrelation in the data.
#acf.tmp <- acf(arima.fit$residuals, plot=FALSE)
# I checked with imodes taking the top 5; the acf immediately falls to 1e-6. Good!

# 8. We can take weekly or monthly block maxima here.





#===============================================================================
#===============================================================================
#===============================================================================
# BELOW HERE IS FROM ORIGINAL fit_priors.R ROUTINE
#===============================================================================
#===============================================================================
#===============================================================================








#
#===============================================================================
# set up GEV and Naveau (i) model parameters
#===============================================================================
#

print('setting up GEV and Naveau-i model parameters for DE optimization...')

source('parameter_setup.R')

print('...done.')

#
#===============================================================================
# parameters for DE optim (for maximum likelihood/minimum negative likelihoood)
#===============================================================================
#

#NP.deoptim <- 200
#niter.deoptim <- 200
NP.deoptim <- NP.deoptim000
niter.deoptim <- niter.deoptim000
F.deoptim <- 0.8
CR.deoptim <- 0.9

#
#===============================================================================
# fit MLE GEV and Naveau models (nav3,gev3=all stationary (i.e., 3 free parameters),
# gev/nav4=location/lower tail parameter is nonstationary, gev/nav5=lower tail
# and scale nonstationary, gev/nav6= all three nonstationary)
# for Dutch station
#===============================================================================
#

print('starting DE optimization for MLE GEV and Naveau-i parameters at Delfzijl, Netherlands...')

source('likelihood_gev.R')
source('likelihood_naveau.R')

deoptim.delfzijl <- vector('list', 8); names(deoptim.delfzijl) <- types.of.model
for (i in 1:length(types.of.model)) {
  deoptim.delfzijl[[types.of.model[i]]] <- mat.or.vec(1, length(parnames_all[[types.of.model[i]]]))
  rownames(deoptim.delfzijl[[types.of.model[i]]]) <- 'delfzijl'
  colnames(deoptim.delfzijl[[types.of.model[i]]]) <- parnames_all[[types.of.model[i]]]
}
bic.delfzijl <- mat.or.vec(1, 8)
colnames(bic.delfzijl) <- types.of.model; rownames(bic.delfzijl) <- 'delfzijl'

# GEV model fitting
for (gev.type in types.of.gev) {
  if(gev.type=='gev3') {auxiliary <- NULL
  } else {auxiliary <- trimmed_forcing(data_calib$year_unique, time_forc, temperature_forc)$temperature}
  out.deoptim <- DEoptim(neg_log_like_gev, lower=bound_lower_set[[gev.type]], upper=bound_upper_set[[gev.type]],
                       DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                       parnames=parnames_all[[gev.type]], data_calib=data_calib$lsl_max, auxiliary=auxiliary)
  deoptim.delfzijl[[gev.type]][1,] <- out.deoptim$optim$bestmem
  colnames(deoptim.delfzijl[[gev.type]]) <- parnames_all[[gev.type]]
  bic.delfzijl[1, gev.type] <- 2*out.deoptim$optim$bestval + length(parnames_all[[gev.type]])*log(length(data_calib$lsl_max))
}

# Naveau (i) model fitting
for (nav.type in types.of.nav) {
  if(nav.type=='nav3') {auxiliary <- NULL
  } else {auxiliary <- trimmed_forcing(data_calib$year_unique, time_forc, temperature_forc)$temperature}
  out.deoptim <- DEoptim(neg_log_like_naveau, lower=bound_lower_set[[nav.type]], upper=bound_upper_set[[nav.type]],
                       DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                       parnames=parnames_all[[nav.type]], data_calib=data_calib$lsl_max, auxiliary=auxiliary)
  deoptim.delfzijl[[nav.type]][1,] <- out.deoptim$optim$bestmem
  colnames(deoptim.delfzijl[[nav.type]]) <- parnames_all[[nav.type]]
  bic.delfzijl[1, nav.type] <- 2*out.deoptim$optim$bestval + length(parnames_all[[nav.type]])*log(length(data_calib$lsl_max))
}

print('...done.')

# plot this survival function against the empirical one
if(FALSE) {
parameters <- deoptim.delfzijl$nav3
x <- seq(0, 10000, 10)
nav.cdf <- naveau_cdf(x=x, kappa=parameters[1], sigma=parameters[2], xi=parameters[3])
plot(esf.levels, log10(esf.values), ylim=c(-2.5,0), xlim=c(0,6000), xlab='Level [cm]', ylab='Survival function [1-cdf]', yaxt='n')
axis(2, at=seq(-3, 0, by=1), label=parse(text=paste("10^", seq(-3,0), sep="")))
lines(x, log10(1-nav.cdf), col='red')
}

#
#===============================================================================
# fit MLE GEV and Naveau models (nav3,gev3=all stationary (i.e., 3 free parameters),
# gev/nav4=location/lower tail parameter is nonstationary, gev/nav5=lower tail
# and scale nonstationary, gev/nav6= all three nonstationary)
# for all tide gauge stations within +/-15 degrees lat/lon of the Dutch station
#===============================================================================
#

print('starting DE optimization for MLE GEV and Naveau-i parameters for all European stations in set...')

deoptim.eur <- vector('list', 8); names(deoptim.eur) <- types.of.model
for (i in 1:length(types.of.model)) {
  deoptim.eur[[types.of.model[i]]] <- mat.or.vec(length(data_set), length(parnames_all[[types.of.model[i]]]))
  rownames(deoptim.eur[[types.of.model[i]]]) <- files.tg
  colnames(deoptim.eur[[types.of.model[i]]]) <- parnames_all[[types.of.model[i]]]
}
bic.eur <- mat.or.vec(length(data_set), 8)
colnames(bic.eur) <- types.of.model; rownames(bic.eur) <- files.tg

for (dd in 1:length(data_set)) {

  data_set[[dd]]$bic.deoptim <- rep(NA, length(types.of.model))
  data_set[[dd]]$deoptim <- vector('list', length(types.of.model))
  names(data_set[[dd]]$deoptim) <- types.of.model

  # GEV model fitting
  for (gev.type in types.of.gev) {
    if(gev.type=='gev3') {auxiliary <- NULL
    } else {auxiliary <- trimmed_forcing(data_set[[dd]]$year_unique, time_forc, temperature_forc)$temperature}
    out.deoptim <- DEoptim(neg_log_like_gev, lower=bound_lower_set[[gev.type]], upper=bound_upper_set[[gev.type]],
                         DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                         parnames=parnames_all[[gev.type]], data_calib=data_set[[dd]]$lsl_max, auxiliary=auxiliary)
    deoptim.eur[[gev.type]][dd,] <- out.deoptim$optim$bestmem
    colnames(deoptim.eur[[gev.type]]) <- parnames_all[[gev.type]]
    bic.eur[dd, gev.type] <- 2*out.deoptim$optim$bestval + length(parnames_all[[gev.type]])*log(length(data_set[[dd]]$lsl_max))
  }

  # Naveau (i) model fitting
  for (nav.type in types.of.nav) {
    if(nav.type=='nav3') {auxiliary <- NULL
    } else {auxiliary <- trimmed_forcing(data_set[[dd]]$year_unique, time_forc, temperature_forc)$temperature}
    out.deoptim <- DEoptim(neg_log_like_naveau, lower=bound_lower_set[[nav.type]], upper=bound_upper_set[[nav.type]],
                         DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                         parnames=parnames_all[[nav.type]], data_calib=data_set[[dd]]$lsl_max, auxiliary=auxiliary)
    deoptim.eur[[nav.type]][dd,] <- out.deoptim$optim$bestmem
    colnames(deoptim.eur[[nav.type]]) <- parnames_all[[nav.type]]
    bic.eur[dd, nav.type] <- 2*out.deoptim$optim$bestval + length(parnames_all[[nav.type]])*log(length(data_set[[dd]]$lsl_max))
  }
}

# also include Delfzijl points in this fit/spread
mle.fits <- vector('list', length(types.of.model)); names(mle.fits) <- types.of.model
mle.fits$gev3 <- rbind(deoptim.eur$gev3, deoptim.delfzijl$gev3)
mle.fits$gev4 <- rbind(deoptim.eur$gev4, deoptim.delfzijl$gev4)
mle.fits$gev5 <- rbind(deoptim.eur$gev5, deoptim.delfzijl$gev5)
mle.fits$gev6 <- rbind(deoptim.eur$gev6, deoptim.delfzijl$gev6)
mle.fits$nav3 <- rbind(deoptim.eur$nav3, deoptim.delfzijl$nav3)
mle.fits$nav4 <- rbind(deoptim.eur$nav4, deoptim.delfzijl$nav4)
mle.fits$nav5 <- rbind(deoptim.eur$nav5, deoptim.delfzijl$nav5)
mle.fits$nav6 <- rbind(deoptim.eur$nav6, deoptim.delfzijl$nav6)

print('...done.')

# plot a fit with the empirical survival function
if(FALSE){
dd <- 6
x.eur <- seq(from=0, to=6000, by=10)
parameters <- deoptim.eur$nav3[dd,]; cdf.nav <- naveau_cdf(x=x.eur, kappa=parameters[1], sigma=parameters[2], xi=parameters[3])
parameters <- deoptim.eur$gev3[dd,]; cdf.gev <- pevd(q=x.eur, loc=parameters[1], scale=parameters[2], shape=parameters[3])
esf.levels.eur <- data_set[[dd]]$lsl_max[order(data_set[[dd]]$lsl_max)]
esf.values.eur <- seq(from=length(data_set[[dd]]$lsl_max), to=1, by=-1)/(length(data_set[[dd]]$lsl_max)+1)
plot(esf.levels.eur, log10(esf.values.eur))
lines(x.eur, log10(1-cdf.nav), col='red')
lines(x.eur, log10(1-cdf.gev), col='blue')
}

#
#===============================================================================
# check distributions, fit priors
#===============================================================================
#

print('fitting prior distributions to the MLE parameters...')

# see where Delfzijl parameter fits lie with respect to the rest of them
# (make sure they are not some extreme outlier) (they are last row of each matrix)
if(FALSE){
for (model in types.of.model) {
  x11()
  par(mfrow=c(3,2))
  for (p in 1:length(parnames_all[[model]])) {
    hist(mle.fits[[model]][,p], xlab=parnames_all[[model]][p], main=model)
    lines(c(mle.fits[[model]][length(data_set)+1,p],mle.fits[[model]][length(data_set)+1,p]), c(-1000,1000), type='l', col='red', lwd=2)
  }
}
}

# fit gamma and normal priors
# -> centered at the medians (or Delfzijl MLEs? no, medians, so it is general *prior* belief)
# -> with standard deviation equal to half the max-min range
#    (or do empirical sd? might underestimate though - take wider)

# assign which parameters have which priors
rm(list=c('gamma.priors','normal.priors','uniform.priors'))
gamma.priors <- c('mu','mu0','kappa','kappa0','sigma','sigma0')
normal.priors <- c('xi','mu1','sigma1','xi0','xi1','kappa1')
uniform.priors <- NULL
# test fitting uniform priors
#uniform.priors <- c('mu','mu0','kappa','kappa0','sigma','sigma0','xi','mu1','sigma1','xi0','xi1','kappa1')

priors <- vector('list', length(types.of.model)); names(priors) <- types.of.model
for (model in types.of.model) {
  priors[[model]] <- vector('list', length(parnames_all[[model]])); names(priors[[model]]) <- parnames_all[[model]]
  for (par in parnames_all[[model]]) {
    priors[[model]][[par]] <- vector('list', 3) # type, and 2 distribution parameters
    if(!is.na(match(par, uniform.priors))) {
       names(priors[[model]][[par]]) <- c('type','shape','rate'); priors[[model]][[par]]$type <- 'uniform'
       priors[[model]][[par]]$lower <- bound_lower_set[[model]][match(par,parnames_all[[model]])]
       priors[[model]][[par]]$upper <- bound_upper_set[[model]][match(par,parnames_all[[model]])]
    } else if(!is.na(match(par, gamma.priors))) { # shape=alpha, rate=beta, mean=shape/rate, var=shape/rate^2
      names(priors[[model]][[par]]) <- c('type','shape','rate'); priors[[model]][[par]]$type <- 'gamma'
      priors[[model]][[par]]$rate <- median(mle.fits[[model]][,par]) / (0.5*(max(mle.fits[[model]][,par])-min(mle.fits[[model]][,par])))^2
      priors[[model]][[par]]$shape <- median(mle.fits[[model]][,par]) * priors[[model]][[par]]$rate
    } else if(!is.na(match(par, normal.priors))) {
      names(priors[[model]][[par]]) <- c('type','mean','sd'); priors[[model]][[par]]$type <- 'normal'
      priors[[model]][[par]]$mean <- median(mle.fits[[model]][,par])
      priors[[model]][[par]]$sd   <- 0.5*(max(mle.fits[[model]][,par])-min(mle.fits[[model]][,par]))
      #priors[[model]][[par]]$sd   <- sd(mle.fits[[model]][,par])
    }
  }
}

# redo the kappas for Naveau, sampling from log(kappa) ('lkappa') instead
model <- 'nav3'; par <- 'kappa'
priors[[model]][[par]]$rate <- median(log(mle.fits[[model]][,par])) / var(log(mle.fits[[model]][,par]))
priors[[model]][[par]]$shape <- median(log(mle.fits[[model]][,par])) * priors[[model]][[par]]$rate
names(priors[[model]])[match(par,names(priors[[model]]))] <- 'lkappa'
parnames_all[[model]][match(par, parnames_all[[model]])] <- 'lkappa'

model <- 'nav4'; par <- 'kappa0'
priors[[model]][[par]]$rate <- median(log(mle.fits[[model]][,par])) / var(log(mle.fits[[model]][,par]))
priors[[model]][[par]]$shape <- median(log(mle.fits[[model]][,par])) * priors[[model]][[par]]$rate
names(priors[[model]])[match(par,names(priors[[model]]))] <- 'lkappa0'
parnames_all[[model]][match(par, parnames_all[[model]])] <- 'lkappa0'

model <- 'nav5'; par <- 'kappa0'
priors[[model]][[par]]$rate <- median(log(mle.fits[[model]][,par])) / var(log(mle.fits[[model]][,par]))
priors[[model]][[par]]$shape <- median(log(mle.fits[[model]][,par])) * priors[[model]][[par]]$rate
names(priors[[model]])[match(par,names(priors[[model]]))] <- 'lkappa0'
parnames_all[[model]][match(par, parnames_all[[model]])] <- 'lkappa0'

model <- 'nav6'; par <- 'kappa0'
priors[[model]][[par]]$rate <- median(log(mle.fits[[model]][,par])) / var(log(mle.fits[[model]][,par]))
priors[[model]][[par]]$shape <- median(log(mle.fits[[model]][,par])) * priors[[model]][[par]]$rate
names(priors[[model]])[match(par,names(priors[[model]]))] <- 'lkappa0'
parnames_all[[model]][match(par, parnames_all[[model]])] <- 'lkappa0'

print('...done.')

# plot priors and MLE histograms
if(FALSE){
for (model in types.of.model) {
  x11()
  par(mfrow=c(3,2))
  for (p in 1:length(parnames_all[[model]])) {

# TODO - deal wtih lkappa and lkappa0 here

    range <- max(mle.fits[[model]][,p]) - min(mle.fits[[model]][,p])
    lower <- min(mle.fits[[model]][,p]) - 0.05*range
    upper <- max(mle.fits[[model]][,p]) + 0.05*range
    if(parnames_all[[model]][p] == 'lkappa0' | parnames_all[[model]][p] == 'lkappa') {
      range <- max(log(mle.fits[[model]][,p])) - min(log(mle.fits[[model]][,p]))
      lower <- min(log(mle.fits[[model]][,p])) - 0.05*range
      upper <- max(log(mle.fits[[model]][,p])) + 0.05*range
    }
    x.tmp <- seq(from=lower, to=upper, length.out=1000)
    if(priors[[model]][[p]]$type=='normal') {
      pdf.tmp <- dnorm(x=x.tmp, mean=priors[[model]][[p]]$mean, sd=priors[[model]][[p]]$sd)
    } else if(priors[[model]][[p]]$type=='gamma') {
      pdf.tmp <- dgamma(x=x.tmp, shape=priors[[model]][[p]]$shape, rate=priors[[model]][[p]]$rate)
    }
    if(parnames_all[[model]][p] == 'lkappa0' | parnames_all[[model]][p] == 'lkappa') {
      hist(log(mle.fits[[model]][,p]), xlab=parnames_all[[model]][p], main=model, freq=FALSE)
      lines(log(c(mle.fits[[model]][length(data_set)+1,p],mle.fits[[model]][length(data_set)+1,p])), c(-1000,1000), type='l', col='red', lwd=2)
    } else {
      hist(mle.fits[[model]][,p], xlab=parnames_all[[model]][p], main=model, freq=FALSE)
      lines(c(mle.fits[[model]][length(data_set)+1,p],mle.fits[[model]][length(data_set)+1,p]), c(-1000,1000), type='l', col='red', lwd=2)
    }
    lines(x.tmp, pdf.tmp, lwd=2, col='blue')
  }
}
}

#
#===============================================================================
# save priors and initial values (from Delfzijl DE optim.) and read later when
# calibrating with MCMC
#===============================================================================
#

# rds -> save single object; the only one we need is 'priors'
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.priors <- paste(output.dir,'surge_priors_',appen,'_',today,'.rds', sep='')
filename.mles <- paste(output.dir,'surge_MLEs_',appen,'_',today,'.rds', sep='')
filename.initvals <- paste(output.dir,'surge_initialvalues_',appen,'_',today,'.rds', sep='')
filename.everything <- paste(output.dir,'everything_priors_',appen,'_',today,'.RData', sep='')
filename.temperature <- paste(output.dir,'temperature_forcing_',today,'.csv', sep='')

print(paste('saving priors and initial values as .rds files (',filename.priors,', ',filename.initvals,') to read and use later...',sep=''))

save.image(file=filename.everything)
saveRDS(priors, file=filename.priors)
saveRDS(mle.fits, file=filename.mles)
saveRDS(deoptim.delfzijl, file=filename.initvals)
write.csv(x=cbind(time_forc, temperature_forc), file=filename.temperature, row.names=FALSE)

print('...done.')

#
#===============================================================================
# End
#===============================================================================
#

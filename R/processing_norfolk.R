#===============================================================================
# processing of Norfolk, Virginia, USA, data.
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================

################# TODO modify!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< TODO here now!
################# TODO modify!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< TODO here now!
################# TODO modify!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< TODO here now!


print('starting to process Delfzijl data')

#===
#=== convert Delfzijl tide gauge data to 3-hourly series
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

#
#===============================================================================
# daily POT series, for PP-GPD model
#===============================================================================
#

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
  sl.daily.max[cnt] <- max(sl.3hour.detrended[ind.today], na.rm=TRUE)
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

ind.exceed <- which(sl.daily.max > gpd.threshold)
time.exceed <- time.daily[ind.exceed]
sl.exceed <- sl.daily.max[ind.exceed]
year.exceed <- year.daily.max[ind.exceed]

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

#
#===============================================================================
# subsample smaller sets of the POT/GPD data
#===============================================================================
#

# initialize, and create list elements for these GPD experiments. then later
# remove from each sub-list item the years the experiment will not use
gpd.experiments <- c('gpd30','gpd50','gpd70','gpd90','gpd110')
years.gpd.experiments <- c(30,50,70,90,110); names(years.gpd.experiments) <- gpd.experiments
for (gpd.exp in gpd.experiments) {
  data_calib[[gpd.exp]] <- data_calib$gpd
  ind.experiment <- (length(data_calib$gev_year$year)-years.gpd.experiments[[gpd.exp]]+1):length(data_calib$gev_year$year)
  data_calib[[gpd.exp]]$counts <- data_calib[[gpd.exp]]$counts[ind.experiment]
  data_calib[[gpd.exp]]$time_length <- data_calib[[gpd.exp]]$time_length[ind.experiment]
  data_calib[[gpd.exp]]$excesses <- data_calib[[gpd.exp]]$excesses[ind.experiment]
  data_calib[[gpd.exp]]$year <- data_calib$gev_year$year[ind.experiment]
  data_calib[[gpd.exp]]$counts_all <- NULL
  data_calib[[gpd.exp]]$time_length_all <- NULL
  data_calib[[gpd.exp]]$excesses_all <- NULL
}

#
#===============================================================================
# now do the GEV/Naveau annual block maxima. calculate based on the 3-houlry
# time series (sl.3hour.detrended)
#===============================================================================
#

print('starting to process Delfzijl data to annual block maxima...')

year.unique <- unique(year.3hour)
year.unique <- year.unique[order(year.unique)]
lsl.mean <- rep(NA, length(year.unique))
lsl.max <- rep(NA, length(year.unique))

for (t in 1:length(year.unique)) {
  ind.this.year <- which(year.3hour==year.unique[t])
  lsl.mean[t] <- mean(sl.3hour.detrended[ind.this.year])
  lsl.max[t] <- max(sl.3hour.detrended[ind.this.year])
}

# skip this for now and clip before calibration, like the other Euro. stations
if(FALSE) {
# clip to only the years that overlap with the historical temperature forcing
# (1880-2016) only clip beginning; have projection at other end.
if(year.unique[1] <= time_forc[1]) {
  ind.clip <- which(year.unique < time_forc[1])
  year.unique <- year.unique[-ind.clip]
  lsl.max <- lsl.max[-ind.clip]
} else if(year.unique[1] > time_forc[1]) {
  ind.clip <- which(time_forc < year.unique[1])
  time_forc <- time_forc[-ind.clip]
  temperature_forc <- temperature_forc[-ind.clip]
}
} # end skip

# set up object for passing through aclibration routines
data_calib$gev_year <- vector('list', 2)
names(data_calib$gev_year) <- c('year','lsl_max')
data_calib$gev_year$year <- year.unique
data_calib$gev_year$lsl_max <- lsl.max

print('  ... done.')

#
#===============================================================================
# now do the GEV/Naveau monthyl block maxima
#===============================================================================
#

# TODO

#
#===============================================================================
# save 'data_calib' object as RDS for use later
#===============================================================================
#

output.dir <- '../output/'
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.datacalib <- paste(output.dir,'datacalib_delfzijl',today,'.rds', sep='')
saveRDS(data_calib, file=filename.datacalib)

save.image(file='../output/preprocessing.RData')

#
#===============================================================================
#

print('done processing Delfzijl data')

#
#===============================================================================
# End
#===============================================================================
#

#===============================================================================
# processing of Delfzijl, the Netherlands tide gauge station data
#
# leads to the list object 'data_delfzijl', analogous to the 'data_many' object
# that contains the many tide gauge stations from UHSLC database;
# has all of the necessary information to estimate (MLE) the PP/GPD parameters
# for each of the 4 candidate model structures.
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================

filename.saveprogress <- '../output/processing_delfzijl.RData'

print('starting to process Delfzijl tide gauge data')

#===
#=== read in tide gauge data
#===

#=====================
# function to read the tide gauge data set(s)
read_data <- function(dat.dir, filetype, septype){
  files.tg <- list.files(path=dat.dir, pattern=filetype)
  data <- read.table(paste(dat.dir,files.tg[1],sep=''), header = TRUE, sep=septype)
  if(length(files.tg) > 1) {
    for (ff in 2:length(files.tg)) {
      data <- rbind(data, read.table(paste(dat.dir,files.tg[ff],sep=''), header = TRUE, sep=septype))
    }
  }
  return(data)
}
#=====================

dat.dir <- '~/codes/EVT/data/Delfzijl_Oddo_data/'
data <- read_data(dat.dir=dat.dir, filetype='txt', septype='\t')

# convert sea levels from cm to mm, consistent with the other data
data$sl <- 10* data$sl

data$year   <- as.numeric(substr(as.character(data$date), start=1, stop=4))
data$month  <- as.numeric(substr(as.character(data$date), start=5, stop=6))
data$day    <- as.numeric(substr(as.character(data$date), start=7, stop=8))
data$hour   <- as.numeric(substr(data$time, 1,2))
data$minute <- as.numeric(substr(data$time, 4,5))

# time in days since 01 January 1960
data$time.days <- as.numeric(mdy.date(month=data$month, day=data$day, year=data$year)) + data$hour/24 + data$minute/(24*60)

# create the object to hold the calibration information for the Delfzijl site
# 2 dimensions, one for 'gev_year' experiment information, and one for 'gpd'
# experiment information
data_delfzijl <- vector('list', 2)
names(data_delfzijl) <- c('gpd','gev_year')

#
#===============================================================================
# daily POT series, for PP-GPD model
#===============================================================================
#


tbeg <- proc.time()

# difference between time stamps (in units of days)
time.diff <- diff(data$time.days)

# check that they are in the proper order, ascending
print(paste('Are there any times out of order? ',any(time.diff < 0), sep=''))

# put everything back in order - make sure you do this for the sea levels and
# other fields too, so do as a matrix. also recalculate the time differences,
# which you will need for averaging
data <- data[order(data$time.days),]
time.diff <- diff(data$time.days)

#===
#=== first, need to average up to 3-hourly time series
#===

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
save.image(file=filename.saveprogress)

#===
#=== subtract linear sea-level trend (from fit to monthly means)
#===

# calculate monthly means

dates.new <- date.mdy(time.days.3hour)
date.beg <- date.mdy(min(time.days.3hour))
date.end <- date.mdy(max(time.days.3hour))

# what the years in which we have data?
years.unique <- unique(dates.new$year)

# in each year, what are the months with at least 90% of the data?
months.this.year <- vector('list', length(years.unique))
names(months.this.year) <- years.unique
years.to.remove <- NULL
for (year in years.unique) {
  ind.this.year <- which(dates.new$year == year)
  months.to.keep <- NULL
  for (month in 1:12) {
    ind.this.month <- which(dates.new$month[ind.this.year] == month)
    days.this.month <- monthDays(paste(year,'-',month,'-','1', sep=''))
    hours.this.month <- days.this.month * 24
    # *3 because these are 3-hourly blocks
    perc.data.this.month <- 3*length(ind.this.month)/hours.this.month
    if (perc.data.this.month >= 0.9) {months.to.keep <- c(months.to.keep, month)}
  }
  if(length(months.to.keep)>0) {months.this.year[[year]] <- months.to.keep }
  else                         {years.to.remove <- c(years.to.remove, year)}
}
if(length(years.to.remove)>0) {years.unique <- years.unique[-match(years.to.remove, years.unique)]}

# get the mean time (in days releative to 1 Jan 1960) of the observations of
# each month we are using to fit the trend for SLR
times.month <- rep(NA, length(unlist(months.this.year)))
sl.month    <- rep(NA, length(unlist(months.this.year)))
cnt <- 1
for (year in years.unique) {
  ind.this.year <- which(dates.new$year == year)
  for (month in months.this.year[[year]]) {
    ind.this.month <- which(dates.new$month[ind.this.year] == month)
    times.month[cnt] <- mean(data$time.days[ind.this.year[ind.this.month]])
    sl.month[cnt]    <- mean(data$sl[ind.this.year[ind.this.month]])
    cnt <- cnt + 1
  }
}

# fit trend
slr.trend <- lm(sl.month ~ times.month)
slr.trend.3hour <- slr.trend$coefficients[1] + slr.trend$coefficients[2]*data$time.days

# subtract off from the 1-hourly data
data$sl.detrended <- data$sl - slr.trend.3hour
print('  ... done.')

#===
#=== daily block maxima; calculate 99% quantile as GPD threshold
#===

# how many days in each year have at least 90% of their values?
days.all <- floor(data$time.days)
days.unique <- unique(days.all)
ind.days.to.remove <- NULL
print('... filtering down to do a daily maxima time series of only the days with at least 90% of data ...')
pb <- txtProgressBar(min=min(days.unique),max=max(days.unique),initial=0,style=3)
for (day in days.unique) {
  ind.today <- which(floor(data$time.days) == day)
  # *3 because 3-hourly series instead of 1-hourly
  perc.data.today <- 3*length(ind.today)/24
  if(perc.data.today < 0.9) {ind.days.to.remove <- c(ind.days.to.remove, match(day, days.unique))}
  setTxtProgressBar(pb, day)
}
close(pb)
days.daily.max <- days.unique[-ind.days.to.remove]
n.days <- length(days.daily.max)

# calculate the daily maximum sea levels on the days of 'days.daily.max'
sl.daily.max <- rep(NA, n.days)
years.daily.max <- rep(NA, n.days)
print('... calculating time series of daily maxima ...')
pb <- txtProgressBar(min=0,max=n.days,initial=0,style=3)
for (day in days.daily.max) {
  cnt <- match(day,days.daily.max)
  ind.today <- which(days.all == day)
  sl.daily.max[cnt] <- max(data$sl.detrended[ind.today])
  years.daily.max[cnt] <- data$year[ind.today][1]
  setTxtProgressBar(pb, cnt)
}
close(pb)

#===
#=== find all the excesses, "declustering" = if two are within a day of each
#=== other, take only the maximum of the two (so make sure you save the times
#=== of each excess)
#===

print('... getting threshold excesses ...')

# threshold is 99% quantile of tide gauge's observed values.
# Buchanan et al (2016) use the 99% quantile of the daily maximum time series.
#gpd.threshold <- as.numeric(quantile(data$sl.detrended, .99, na.rm=TRUE))
gpd.threshold <- as.numeric(quantile(sl.daily.max, .99, na.rm=TRUE))

ind.exceed <- which(sl.daily.max > gpd.threshold)
days.exceed <- days.daily.max[ind.exceed]
sl.exceed <- sl.daily.max[ind.exceed]
years.exceed <- years.daily.max[ind.exceed]

declustered.exceed <- decluster_timeseries(time=days.exceed, year=years.exceed, time.series=sl.exceed, min.dt=1)
days.exceed.decl <- declustered.exceed$time
years.exceed.decl <- declustered.exceed$year
sl.exceed.decl <- declustered.exceed$time.series

#===
#=== sub-list object on 'data_delfzijl' to store what is needed to calbirate PP-GPD
#===

data_delfzijl$gpd <- vector('list', 5)
names(data_delfzijl$gpd) <- c('counts','year','time_length','excesses','threshold')
# initialize
data_delfzijl$gpd$threshold <- gpd.threshold
data_delfzijl$gpd$counts <- data_delfzijl$gpd$year <- data_delfzijl$gpd$time_length <- rep(NA, length(years.unique))
data_delfzijl$gpd$excesses <- vector('list', length(years.unique))

for (ind.year in 1:length(years.unique)) {
  ind.hits.this.year <- which(years.exceed.decl == years.unique[ind.year])
  data_delfzijl$gpd$counts[ind.year] <- length(ind.hits.this.year)
  data_delfzijl$gpd$year[ind.year]   <- years.unique[ind.year]
  data_delfzijl$gpd$time_length[ind.year] <- length(which(years.daily.max == years.unique[ind.year]))
  if(length(ind.hits.this.year) > 0) {data_delfzijl$gpd$excesses[[ind.year]] <- sl.exceed.decl[ind.hits.this.year]
  } else                             {data_delfzijl$gpd$excesses[[ind.year]] <- NA}
}

# alternatively, could bin em all together. but this won't allow for potential
# non-stationary behavior in the poisson process
data_delfzijl$gpd$excesses_all <- sl.exceed.decl
data_delfzijl$gpd$counts_all <- length(sl.exceed.decl)
data_delfzijl$gpd$time_length_all <- length(days.daily.max)

# that takes some time, so save the workspace image after each data set
save.image(file=filename.saveprogress)

#
#===============================================================================
# subsample smaller sets of the POT/GPD data
#===============================================================================
#

# initialize, and create list elements for these GPD experiments. then later
# remove from each sub-list item the years the experiment will not use
gpd.experiments <- c('gpd30','gpd50','gpd70','gpd90','gpd110','gpd137')
years.gpd.experiments <- c(30,50,70,90,110,137); names(years.gpd.experiments) <- gpd.experiments
for (gpd.exp in gpd.experiments) {
  data_delfzijl[[gpd.exp]] <- data_delfzijl$gpd
  ind.experiment <- (length(data_delfzijl$gpd$year)-years.gpd.experiments[[gpd.exp]]+1):length(data_delfzijl$gpd$year)
  data_delfzijl[[gpd.exp]]$counts <- data_delfzijl[[gpd.exp]]$counts[ind.experiment]
  data_delfzijl[[gpd.exp]]$time_length <- data_delfzijl[[gpd.exp]]$time_length[ind.experiment]
  data_delfzijl[[gpd.exp]]$excesses <- data_delfzijl[[gpd.exp]]$excesses[ind.experiment]
  data_delfzijl[[gpd.exp]]$year <- data_delfzijl$gpd$year[ind.experiment]
  data_delfzijl[[gpd.exp]]$counts_all <- NULL
  data_delfzijl[[gpd.exp]]$time_length_all <- NULL
  data_delfzijl[[gpd.exp]]$excesses_all <- NULL
}

tend <- proc.time()
print(paste('  ... done. Took ', (tend[3]-tbeg[3])/60, ' minutes.',sep=''))

#
#===============================================================================
# now do the GEV/Naveau annual block maxima. calculate based on the 3-hourly time
# series (data$sl.detrended)
#===============================================================================
#

# give an update to the screen
print('starting to process annual block maxima for Delfzijl data set')

# set up object for passing through calibration routines
data_delfzijl$gev_year <- vector('list', 2)
names(data_delfzijl$gev_year) <- c('year','lsl_max')

data_delfzijl$gev_year$year <- unique(data$year)
data_delfzijl$gev_year$year <- data_delfzijl$gev_year$year[order(data_delfzijl$gev_year$year)]
data_delfzijl$gev_year$lsl_max <- rep(NA, length(data_delfzijl$gev_year$year))
for (t in 1:length(data_delfzijl$gev_year$year)) {
  ind_this_year <- which(data$year==data_delfzijl$gev_year$year[t])
  data_delfzijl$gev_year$lsl_max[t] <- max(data$sl.detrended[ind_this_year])
}

# trim before 1850 (or whenever is time_forc[1]), which is when the temperature
# forcing starts. also make a note of how long each record is
nyear.delfzijl <- NA
if(data_delfzijl$gev_year$year[1] < time_forc[1]) {
  icut <- which(data_delfzijl$gev_year$year < time_forc[1])
  data_delfzijl$gev_year$year <- data_delfzijl$gev_year$year[-icut]
  data_delfzijl$gev_year$lsl_max <- data_delfzijl$gev_year$lsl_max[-icut]
}
nyear.delfzijl <- length(data_delfzijl$gev_year$year)

# that doesn't take as long... so just save it once for good measure
save.image(file=filename.saveprogress)

# save final 'data_delfzijl' object to RDS to use later
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.output <- paste('../data/tidegauge_processed_delfzijl_',today,'.rds', sep='')
saveRDS(data_delfzijl, file=filename.output)

#===============================================================================

print('done processing the Delfzijl tide gauge data set')

#===============================================================================
# End
#===============================================================================

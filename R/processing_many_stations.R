#===============================================================================
# processing of all tide gauge stations' data, with at least 90 years.
#
# leads to the list object 'data_many', which has 28 elements, each corresponding
# to a different tide gauge station from UHSLC database with at least 90 years of
# data (not necessarily continuous, and not necessarily counting for gaps).
# the data_many[[i]] object corresponds to the i-th station (i=1,2,...,28), and
# has all of the necessary information to estimate (MLE) the PP/GPD parameters
# for each of the 4 candidate model structures.
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================

l.installpackages <- FALSE

if(l.installpackages) {
  install.packages('date')
  install.packages('zoo')
  install.packages('Hmisc')
}
library(date)
library(zoo)
library(Hmisc)

print('starting to process many tide gauge stations data')

#===
#=== read in tide gauge data, they're all already hourly series
#===

filetype='csv'
dat.dir <- '~/codes/EVT/data/tide_gauge_long/'
files.tg <- list.files(path=dat.dir, pattern=filetype)

data_set <- vector('list', length(files.tg))
for (dd in 1:length(files.tg)) {
  # print an update of progress to the screen
  print(paste('now reading in data set ',dd,' / ',length(files.tg),sep=''))
  names(data_set)[dd] <- substr(files.tg[dd], start=1, stop=7)
  data.tmp <- read.table(paste(dat.dir,files.tg[dd],sep=''), header=FALSE, sep=',')
  data_set[[dd]] <- vector('list', 5)
  names(data_set[[dd]]) <- c('year','month','day','hour','sl')
  data_set[[dd]]$year <- data.tmp$V1
  data_set[[dd]]$month <- data.tmp$V2
  data_set[[dd]]$day <- data.tmp$V3
  data_set[[dd]]$hour <- data.tmp$V4
  data_set[[dd]]$sl <- data.tmp$V5
  # time in days since 01 January 1960
  data_set[[dd]]$time.days <- as.numeric(mdy.date(month=data_set[[dd]]$month, day=data_set[[dd]]$day,
                                                  year=data_set[[dd]]$year)) + data_set[[dd]]$hour/24
}

# create the object to hold the calibration information for the european tide gauges
data_many <- vector('list', length(data_set))
names(data_many) <- names(data_set)


# note: this is all broken up into different chunks of code and processing so
# the it is digestible and won't die if a single part fo the processing is
# broken.


#
#===============================================================================
# daily POT series, for PP-GPD model
#===============================================================================
#

# difference between time stamps (in units of days)
for (dd in 1:length(files.tg)) {

  print(paste('now processing for pp-gpd data set ',dd,' / ',length(files.tg),sep=''))

  time.diff <- diff(data_set[[dd]]$time.days)

  # check that they are in the proper order, ascending
  print(paste('Are there any times out of order? ',any(time.diff < 0), sep=''))

  # put everything back in order - make sure you do this for the sea levels and
  # other fields too, so do as a matrix. also recalculate the time differences,
  # which you will need for averaging
# Note -  data from PMSLC is all in order
#  data_set[[dd]] <- data_set[[dd]][order(data_set[[dd]]$time.days),]
#  time.diff <- diff(data_set[[dd]]$time.days)

  # what is one hour? in units of days
  one.hours <- 1/24

  # where are there gaps longer than one hour? (+10sec for precision)
  igap <- which(time.diff > (one.hours+10/(24*60*60)))

  #===
  #=== subtract linear sea-level trend (from fit to monthly means)
  #===

  # calcualte monthly means

  dates.new <- date.mdy(data_set[[dd]]$time.days)
  date.beg <- date.mdy(min(data_set[[dd]]$time.days))
  date.end <- date.mdy(max(data_set[[dd]]$time.days))

  # what the years in which we have data?
  years.unique <- unique(dates.new$year)

  # in each year, what are the months with at least 90% of the data?
  months.this.year <- vector('list', length(years.unique))
  names(months.this.year) <- years.unique
  for (year in years.unique) {
    ind.this.year <- which(dates.new$year == year)
    months.to.keep <- NULL
    for (month in 1:12) {
      ind.this.month <- which(dates.new$month[ind.this.year] == month)
      days.this.month <- monthDays(paste(year,'-',month,'-','1', sep=''))
      hours.this.month <- days.this.month * 24
      perc.data.this.month <- length(ind.this.month)/hours.this.month
      if (perc.data.this.month >= 0.9) {months.to.keep <- c(months.to.keep, month)}
    }
    months.this.year[[year]] <- months.to.keep
  }

  # get the mean time (in days releative to 1 Jan 1960) of the observations of
  # each month we are using to fit the trend for SLR
  times.month <- rep(NA, length(unlist(months.this.year)))
  sl.month    <- rep(NA, length(unlist(months.this.year)))
  cnt <- 1
  for (year in years.unique) {
    ind.this.year <- which(dates.new$year == year)
    for (month in months.this.year[[year]]) {
      ind.this.month <- which(dates.new$month[ind.this.year] == month)
      times.month[cnt] <- mean(data_set[[dd]]$time.days[ind.this.year[ind.this.month]])
      sl.month[cnt]    <- mean(data_set[[dd]]$sl[ind.this.year[ind.this.month]])
      cnt <- cnt + 1
    }
  }

  # fit trend
  slr.trend <- lm(sl.month ~ times.month)
  slr.trend.1hour <- slr.trend$coefficients[1] + slr.trend$coefficients[2]*data_set[[dd]]$time.days

  # subtract off from the 1-hourly data
  data_set[[dd]]$sl.detrended <- data_set[[dd]]$sl - slr.trend.1hour
  print('  ... done.')

  ## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< TODO HERE NOW!
  ## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< TODO HERE NOW!
  ##
  ## This needs to be modified such that we use only the information available from each
  ## unique year, and save the length of time from that year (for the PP rate calculation).
  ## Have a "data_calib$time_length" element in "likelihood_ppgpd.R", so use it!
  ##
  ## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< TODO HERE NOW!
  ## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< TODO HERE NOW!

  #===
  #=== daily block maxima; calculate 99% quantile as GPD threshold
  #===

  # how many days in each year have at least 90% of their values?
  days.unique <- unique(floor(data_set[[dd]]$time.days))
  for (day in days.unique) {
    ind.today <- which(floor(data_set[[dd]]$time.days) == day)
  }
  sl.daily.max <- rep(NA, n.days)



  year.daily.max <- rep(NA, n.days) # save which year each came from to organize into blocks
  time.daily <- seq(from=0.5+floor(min(data_set[[dd]]$time.days)), to=0.5+floor(max(data_set[[dd]]$time.days)))
  cnt <- 1
  print('getting time series of daily maxima...')
  pb <- txtProgressBar(min=0,max=n.days,initial=0,style=3)
  for (ttmp in time.daily) {
    ind.today <- which( abs(data_set[[dd]]$time.days-ttmp) < 0.5 )
    if(length(ind.today) > 0) {
      sl.daily.max[cnt] <- max(data_set[[dd]]$sl.detrended[ind.today], na.rm=TRUE)
      year.daily.max[cnt] <- median(data_set[[dd]]$year[ind.today], na.rm=TRUE)
    } else {
      sl.daily.max[cnt] <- NA
      year.daily.max[cnt] <- NA
    }
    setTxtProgressBar(pb, cnt)
    cnt <- cnt+1
  }
  close(pb)
  print('  ... done.')

  #===
  #=== find all the excesses, "declustering" = if two are within a day of each
  #=== other, take only the maximum of the two (so make sure you save the times
  #=== of each excess)
  #===

  print('getting threshold excesses...')

  # threshold is 99% quantile of tide gauge's observed values.
  # Buchanan et al (2016) use the 99% quantile of the daily maximum time series.
  #gpd.threshold <- as.numeric(quantile(data_set[[dd]]$sl.detrended, .99, na.rm=TRUE))
  gpd.threshold <- as.numeric(quantile(sl.daily.max, .99, na.rm=TRUE))

  ind.exceed <- which(sl.daily.max > gpd.threshold)
  time.exceed <- time.daily[ind.exceed]
  sl.exceed <- sl.daily.max[ind.exceed]
  year.exceed <- year.daily.max[ind.exceed]

  declustered.exceed <- decluster_timeseries(time=time.exceed, year=year.exceed, time.series=sl.exceed, min.dt=1)
  time.exceed.decl <- declustered.exceed$time
  year.exceed.decl <- declustered.exceed$year
  sl.exceed.decl <- declustered.exceed$time.series

  #===
  #=== sub-list object on 'data_many' to store what is needed to calbirate PP-GPD
  #===

  year.unique <- unique(data_set[[dd]]$year)
  year.unique <- year.unique[order(year.unique)]
  data_many[[dd]]$gpd <- vector('list', 4)
  names(data_many[[dd]]$gpd) <- c('counts','time_length','excesses','threshold')
  # initialize
  data_many[[dd]]$gpd$threshold <- gpd.threshold
  data_many[[dd]]$gpd$counts <- data_many[[dd]]$gpd$time_length <- rep(NA, length(year.unique))
  data_many[[dd]]$gpd$excesses <- vector('list', length(year.unique))
  for (y in 1:length(year.unique)) {
    hits.this.year <- which(year.exceed.decl==year.unique[y])
    data_many[[dd]]$gpd$counts[y] <- length(hits.this.year)
    data_many[[dd]]$gpd$time_length[y] <- max(time.daily[which(year.daily.max==year.unique[y])]) -
                                     min(time.daily[which(year.daily.max==year.unique[y])]) + 1
    if(length(hits.this.year) > 0) {data_many[[dd]]$gpd$excesses[[y]] <- sl.exceed.decl[hits.this.year]
    } else                         {data_many[[dd]]$gpd$excesses[[y]] <- NA}
  }

  # alternatively, could bin em all together. but this won't allow for potential
  # non-stationary behavior in the poisson process
  data_many[[dd]]$gpd$excesses_all <- sl.exceed.decl
  data_many[[dd]]$gpd$counts_all <- length(sl.exceed.decl)
  data_many[[dd]]$gpd$time_length_all <- max(time.daily) - min(time.daily) + 1

  # that takes some time, so save the workspace image after each data set
  save.image(file='../output/preprocessing.RData')

  print('  ... done.')

}


#
#===============================================================================
# now do the GEV/Naveau annual block maxima. calculate based on the hourly time
# series (data_set[[dd]]$sl.detrended)
#===============================================================================
#

for (dd in 1:length(data_set)) {
  # give an update to the screen
  print(paste('starting to process annual block maxima for European data set ',dd,' / ',length(data_set),sep=''))

  # set up object for passing through aclibration routines
  data_many[[dd]]$gev_year <- vector('list', 2)
  names(data_many[[dd]]$gev_year) <- c('year','lsl_max')

  data_many[[dd]]$gev_year$year <- unique(data_set[[dd]]$year)
  data_many[[dd]]$gev_year$year <- data_many[[dd]]$gev_year$year[order(data_many[[dd]]$gev_year$year)]
  data_many[[dd]]$gev_year$lsl_max <- rep(NA, length(data_many[[dd]]$gev_year$year))
  for (t in 1:length(data_many[[dd]]$gev_year$year)) {
    ind_this_year <- which(data_set[[dd]]$year==data_many[[dd]]$gev_year$year[t])
    data_many[[dd]]$gev_year$lsl_max[t] <- max(data_set[[dd]]$sl.detrended[ind_this_year])
  }
}


# trim before 1850 (or whenever is time_forc[1]), which is when the temperature
# forcing starts. also make a note of how long each record is
nyear <- rep(NA, length(data_many))
for (dd in 1:length(data_many)) {
  if(data_many[[dd]]$gev_year$year[1] < time_forc[1]) {
    icut <- which(data_many[[dd]]$gev_year$year < time_forc[1])
    data_many[[dd]]$gev_year$year_unique <- data_many[[dd]]$gev_year$year_unique[-icut]
    data_many[[dd]]$gev_year$lsl_max <- data_many[[dd]]$gev_year$lsl_max[-icut]
  }
  nyear[dd] <- length(data_many[[dd]]$gev_year$year_unique)
}

# that doesn't take as long... so just save it once for good measure
save.image(file='../output/preprocessing.RData')


#
#===============================================================================
# now do the GEV/Naveau monthyl block maxima
#===============================================================================
#

# TODO


#===============================================================================

print('done processing the European data sets')

#===============================================================================
# End
#===============================================================================

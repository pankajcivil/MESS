#===============================================================================
# processing of many European tide gauge stations data.
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================

print('starting to process many European station data')

#===
#=== read in European tide gauge data, they're all already hourly series
#===

filetype='csv'
dat.dir <- '~/codes/EVT/data/tide_gauge_Europe/'
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
data_europe <- vector('list', length(data_set))
names(data_europe) <- names(data_set)


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
# European data from PMSLC is all in order
#  data_set[[dd]] <- data_set[[dd]][order(data_set[[dd]]$time.days),]
#  time.diff <- diff(data_set[[dd]]$time.days)

  # what is one hour? in units of days
  one.hours <- 1/24

  # where are there gaps longer than one hour? (+10sec for precision)
  igap <- which(time.diff > (one.hours+10/(24*60*60)))

  # where are there gaps shorter than one hours? (-10sec)
  iokay <- which( (time.diff < (one.hours+10/(24*60*60))) & (time.diff > (one.hours-10/(24*60*60))) )
  ishort <- which(time.diff < (one.hours-10/(24*60*60)))

  #===
  #=== subtract linear sea-level trend (from fit to monthly means)
  #===

  # calcualte monthly means

  dates.new <- date.mdy(data_set[[dd]]$time.days)
  date.beg <- date.mdy(min(data_set[[dd]]$time.days))
  date.end <- date.mdy(max(data_set[[dd]]$time.days))

  # how many months are there?
  # this does not account for the fact that there may be missing months, or years
  # this possibility is dealt with below, where we just fill in NA for the
  # missing months
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
    if(length(ind.this.year) > 0) {
      for (mm in min(dates.new$month[ind.this.year]):max(dates.new$month[ind.this.year])) {
        ind.this.month <- which( (dates.new$year == yy) & (dates.new$month == mm) )
        if(length(ind.this.month) > 0) {
          sl.month[cnt] <- mean(data_set[[dd]]$sl[ind.this.month], na.rm=TRUE)
          time.month[cnt] <- mean(data_set[[dd]]$time.days[ind.this.month], na.rm=TRUE)
        } else {
          sl.month[cnt] <- NA
          time.month[cnt] <- NA
        }
        cnt <- cnt+1
        setTxtProgressBar(pb, cnt)
      }
    } else {
      # can use 1:12 here because neither the first nor last year would have
      # no data points; otherwise it wouldn't be on the data set!
      for (mm in 1:12) {
        sl.month[cnt] <- NA
        time.month[cnt] <- NA
        cnt <- cnt+1
        setTxtProgressBar(pb, cnt)
      }
    }
  }
  close(pb)
  print('  ... done.')

  # fit trend
  slr.trend <- lm(sl.month ~ time.month)
  slr.trend.1hour <- slr.trend$coefficients[1] + slr.trend$coefficients[2]*data_set[[dd]]$time.days

  # subtract off from the 1-hourly data
  data_set[[dd]]$sl.detrended <- data_set[[dd]]$sl - slr.trend.1hour
  print('  ... done.')

  #===
  #=== daily block maxima; calculate 99% quantile as GPD threshold
  #===

  n.days <- ceiling(max(data_set[[dd]]$time.days) - min(data_set[[dd]]$time.days))
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
  #=== sub-list object on 'data_europe' to store what is needed to calbirate PP-GPD
  #===

  year.unique <- unique(data_set[[dd]]$year)
  year.unique <- year.unique[order(year.unique)]
  data_europe[[dd]]$gpd <- vector('list', 4)
  names(data_europe[[dd]]$gpd) <- c('counts','time_length','excesses','threshold')
  # initialize
  data_europe[[dd]]$gpd$threshold <- gpd.threshold
  data_europe[[dd]]$gpd$counts <- data_europe[[dd]]$gpd$time_length <- rep(NA, length(year.unique))
  data_europe[[dd]]$gpd$excesses <- vector('list', length(year.unique))
  for (y in 1:length(year.unique)) {
    hits.this.year <- which(year.exceed.decl==year.unique[y])
    data_europe[[dd]]$gpd$counts[y] <- length(hits.this.year)
    data_europe[[dd]]$gpd$time_length[y] <- max(time.daily[which(year.daily.max==year.unique[y])]) -
                                     min(time.daily[which(year.daily.max==year.unique[y])]) + 1
    if(length(hits.this.year) > 0) {data_europe[[dd]]$gpd$excesses[[y]] <- sl.exceed.decl[hits.this.year]
    } else                         {data_europe[[dd]]$gpd$excesses[[y]] <- NA}
  }

  # alternatively, could bin em all together. but this won't allow for potential
  # non-stationary behavior in the poisson process
  data_europe[[dd]]$gpd$excesses_all <- sl.exceed.decl
  data_europe[[dd]]$gpd$counts_all <- length(sl.exceed.decl)
  data_europe[[dd]]$gpd$time_length_all <- max(time.daily) - min(time.daily) + 1

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
  data_europe[[dd]]$gev_year <- vector('list', 2)
  names(data_europe[[dd]]$gev_year) <- c('year','lsl_max')

  data_europe[[dd]]$gev_year$year <- unique(data_set[[dd]]$year)
  data_europe[[dd]]$gev_year$year <- data_europe[[dd]]$gev_year$year[order(data_europe[[dd]]$gev_year$year)]
  data_europe[[dd]]$gev_year$lsl_max <- rep(NA, length(data_europe[[dd]]$gev_year$year))
  for (t in 1:length(data_europe[[dd]]$gev_year$year)) {
    ind_this_year <- which(data_set[[dd]]$year==data_europe[[dd]]$gev_year$year[t])
    data_europe[[dd]]$gev_year$lsl_max[t] <- max(data_set[[dd]]$sl.detrended[ind_this_year])
  }
}


# trim before 1850 (or whenever is time_forc[1]), which is when the temperature
# forcing starts. also make a note of how long each record is
nyear <- rep(NA, length(data_europe))
for (dd in 1:length(data_europe)) {
  if(data_europe[[dd]]$gev_year$year[1] < time_forc[1]) {
    icut <- which(data_europe[[dd]]$gev_year$year < time_forc[1])
    data_europe[[dd]]$gev_year$year_unique <- data_europe[[dd]]$gev_year$year_unique[-icut]
    data_europe[[dd]]$gev_year$lsl_max <- data_europe[[dd]]$gev_year$lsl_max[-icut]
  }
  nyear[dd] <- length(data_europe[[dd]]$gev_year$year_unique)
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

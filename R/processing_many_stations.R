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

filename.saveprogress <- '../output/processing_many_sites.RData'

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

  tbeg <- proc.time()
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
  years.to.remove <- NULL
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

  #===
  #=== daily block maxima; calculate 99% quantile as GPD threshold
  #===

  # how many days in each year have at least 90% of their values?
  days.all <- floor(data_set[[dd]]$time.days)
  days.unique <- unique(days.all)
  ind.days.to.remove <- NULL
  print('... filtering down to do a daily maxima time series of only the days with at least 90% of data ...')
  pb <- txtProgressBar(min=min(days.unique),max=max(days.unique),initial=0,style=3)
  for (day in days.unique) {
    ind.today <- which(floor(data_set[[dd]]$time.days) == day)
    perc.data.today <- length(ind.today)/24
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
    sl.daily.max[cnt] <- max(data_set[[dd]]$sl.detrended[ind.today])
    years.daily.max[cnt] <- data_set[[dd]]$year[ind.today][1]
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
  #gpd.threshold <- as.numeric(quantile(data_set[[dd]]$sl.detrended, .99, na.rm=TRUE))
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
  #=== sub-list object on 'data_many' to store what is needed to calbirate PP-GPD
  #===

  data_many[[dd]]$gpd <- vector('list', 5)
  names(data_many[[dd]]$gpd) <- c('counts','year','time_length','excesses','threshold')
  # initialize
  data_many[[dd]]$gpd$threshold <- gpd.threshold
  data_many[[dd]]$gpd$counts <- data_many[[dd]]$gpd$year <- data_many[[dd]]$gpd$time_length <- rep(NA, length(years.unique))
  data_many[[dd]]$gpd$excesses <- vector('list', length(years.unique))

  for (ind.year in 1:length(years.unique)) {
    ind.hits.this.year <- which(years.exceed.decl == years.unique[ind.year])
    data_many[[dd]]$gpd$counts[ind.year] <- length(ind.hits.this.year)
    data_many[[dd]]$gpd$year[ind.year]   <- years.unique[ind.year]
    data_many[[dd]]$gpd$time_length[ind.year] <- length(which(years.daily.max == years.unique[ind.year]))
    if(length(ind.hits.this.year) > 0) {data_many[[dd]]$gpd$excesses[[ind.year]] <- sl.exceed.decl[ind.hits.this.year]
    } else                             {data_many[[dd]]$gpd$excesses[[ind.year]] <- NA}
  }

  # alternatively, could bin em all together. but this won't allow for potential
  # non-stationary behavior in the poisson process
  data_many[[dd]]$gpd$excesses_all <- sl.exceed.decl
  data_many[[dd]]$gpd$counts_all <- length(sl.exceed.decl)
  data_many[[dd]]$gpd$time_length_all <- length(days.daily.max)

  # that takes some time, so save the workspace image after each data set
  save.image(file=filename.saveprogress)

  tend <- proc.time()
  print(paste('  ... done. Took ', (tend[3]-tbeg[3])/60, ' minutes.',sep=''))

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
    data_many[[dd]]$gev_year$year <- data_many[[dd]]$gev_year$year[-icut]
    data_many[[dd]]$gev_year$lsl_max <- data_many[[dd]]$gev_year$lsl_max[-icut]
  }
  nyear[dd] <- length(data_many[[dd]]$gev_year$year)
}

# that doesn't take as long... so just save it once for good measure
save.image(file=filename.saveprogress)

# save final 'data_many' object to RDS to use later
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.output <- paste('../data/tidegauge_processed_manystations_',today,'.rds', sep='')
saveRDS(data_many, file=filename.output)

#===============================================================================

print('done processing the long tide gauge data sets')

#===============================================================================
# End
#===============================================================================

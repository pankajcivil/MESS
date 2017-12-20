#===============================================================================
# processing of Delfzijl, the Netherlands tide gauge station data
#
# leads to the list object 'data_delfzijl', analogous to the 'data_many' object
# that contains the many tide gauge stations from UHSLC database;
# has all of the necessary information to estimate (MLE) the PP/GPD parameters
# for each of the 4 candidate model structures.
#
# File paths expect you are in the EVT/R directory
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================
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

processing_delfzijl <- function(dt.decluster, detrend.method, pot.threshold) {

  filename.saveprogress <- '../output/processing_delfzijl.RData'

  print('starting to process Delfzijl tide gauge data')

  #===
  #=== read in tide gauge data
  #===

  data <- read.table('../data/id1-DELFZL-187901010000-201701012359.txt', header = TRUE, sep=';', skip=4)
  data <- data[,c(3,4,6)]
  names(data) <- c('date','time','sl')
  data$sl <- data$sl*10 # convert to mm from cm

  # separate date-stamp into year / month / day
  data$year   <- as.numeric(substr(as.character(data$date), start=1, stop=4))
  data$month  <- as.numeric(substr(as.character(data$date), start=6, stop=7))
  data$day    <- as.numeric(substr(as.character(data$date), start=9, stop=10))
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
  #=== first, need to *subsample* up to 3-hourly time series (c.f. email with Rijkwaterstaat, Alexander Bakker)
  #===

  # what is three hours? in units of days
  three.hours <- 3/24

  # where are there gaps longer than three hours? (+10sec for precision)
  igap <- which(time.diff > (three.hours+10/(24*60*60)))

  # where are there gaps shorter than three hours? (-10sec)
  # should find lots of places. need to average these up to 3 hours.
  iokay <- which( (time.diff < (three.hours+10/(24*60*60))) & (time.diff > (three.hours-10/(24*60*60))) )
  ishort <- which(time.diff < (three.hours-10/(24*60*60)))

  # turns out with the Rijkwaterstaat data set, have a lot of different time
  # intervals, but no gaps longer than 3 hours.  upscale the less-than-3-hourly
  # data to 3-hourly

  # in a perfect world, these are the times at which we have obs, and we just
  # need to pluck out the observations at these time points
  time_beg <- min(data$time.days)
  time_end <- max(data$time.days)
  delta_t <- 3/24   # time step of 3 hours, units of days
  time_3hour <- seq(from=time_beg, to=time_end, by=delta_t)

  # go through the entire time.days vector; you know there are no obs gaps longer
  # than 3 hours, so each of the 3-hour intervals has at least one observation in
  # it. grab the one closest to the mark (without going over?)
  # warning: this will take a while (~1 hour, maybe a bit more, on laptop)

  ind_map3hour <- sapply(1:length(time_3hour), function(t) {which.min(abs(time_3hour[t]-data$time.days))})
  sl_3hour <- data$sl[ind_map3hour]
  year_3hour <- data$year[ind_map3hour]

  # that takes a long time, so save the workspace image
  save.image(file=filename.saveprogress)

  #===
  #=== Detrend by either subtracting linear sea-level trend (fit to monthly
  #=== means) or by subtracting annual means (moving 1-year window)
  #===

  print(paste('Detrending using method `',detrend.method,'` ...', sep=''))

  if (detrend.method=='linear') {

    # calculate monthly means

    dates.new <- date.mdy(time_3hour)
    date.beg <- date.mdy(min(time_3hour))
    date.end <- date.mdy(max(time_3hour))

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

  } else if(detrend.method=='annual') {

    # what the years in which we have data?
    dates.new <- date.mdy(time_3hour)
    years.unique <- unique(dates.new$year)

    # get a placeholder -- want to be using the 3-hourly time series
    sl_3hour_detrended <- sl_3hour
    time.days.beg <- min(time_3hour)
    time.days.end <- max(time_3hour)

    pb <- txtProgressBar(min=0,max=length(time_3hour),initial=0,style=3)
    for (tt in 1:length(time_3hour)) {
      # if within half a year of either end of the time series, include either the
      # entire first year or entire last year to get a full year's worth of data in
      # the subtracted mean
      if (time_3hour[tt] - time.days.beg < (365.25*0.5)) {
        ind.close <- which(time_3hour - time.days.beg <= 365.25)
      } else if(time.days.end - time_3hour[tt] < (365.25*0.5)) {
        ind.close <- which(time.days.end - time_3hour <= 365.25)
      } else {
        ind.close <- which(abs(time_3hour-time_3hour[tt]) <= (365.25*0.5) )
      }
      sl_3hour_detrended[tt] <- sl_3hour[tt] - mean(sl_3hour[ind.close])
      setTxtProgressBar(pb, tt)
    }
    close(pb)
  } else {
    print('ERROR: unknown detrend.method value')
  }

  print('  ... done.')

  #===
  #=== daily block maxima; calculate 99% quantile as GPD threshold
  #===

  # how many days in each year have at least 90% of their values?
  days.all <- floor(time_3hour)
  days.unique <- unique(days.all)
  ind.days.to.remove <- NULL
  print('... filtering down to do a daily maxima time series of only the days with at least 90% of data ...')
  pb <- txtProgressBar(min=min(days.unique),max=max(days.unique),initial=0,style=3)
  for (day in days.unique) {
    ind.today <- which(floor(time_3hour) == day)
    # *3 because 3-hourly series instead of 1-hourly
    perc.data.today <- 3*length(ind.today)/24
    if(perc.data.today < 0.9) {ind.days.to.remove <- c(ind.days.to.remove, match(day, days.unique))}
    setTxtProgressBar(pb, day)
  }
  close(pb)
  days.daily.max <- days.unique[-ind.days.to.remove]
  n.days <- length(days.daily.max)


      TODO HERE NOW
      TODO HERE NOW
      TODO HERE NOW




  # calculate the daily maximum sea levels on the days of 'days.daily.max'
  sl.daily.max <- rep(NA, n.days)
  years.daily.max <- rep(NA, n.days)
  print('... calculating time series of daily maxima ...')
  pb <- txtProgressBar(min=0,max=n.days,initial=0,style=3)
  for (day in days.daily.max) {
    cnt <- match(day,days.daily.max)
    ind.today <- which(days.all == day)
    sl.daily.max[cnt] <- max(sl_3hour_detrended[ind.today])
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
  #gpd.threshold <- as.numeric(quantile(sl.daily.max, .99, na.rm=TRUE))
  gpd.threshold <- as.numeric(quantile(sl.daily.max, pot.threshold, na.rm=TRUE))
  data_delfzijl$dt.decluster <- dt.decluster

  ind.exceed <- which(sl.daily.max > gpd.threshold)
  days.exceed <- days.daily.max[ind.exceed]
  sl.exceed <- sl.daily.max[ind.exceed]
  years.exceed <- years.daily.max[ind.exceed]

  declustered.exceed <- decluster_timeseries(time=days.exceed, year=years.exceed, time.series=sl.exceed, min.dt=dt.decluster)
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
  data_delfzijl$gpd$p.threshold <- pot.threshold
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
  filename.output <- paste('../data/tidegauge_processed_deflzijl_decl',data_delfzijl$dt.decluster,'-pot',pot.threshold*100,'-',detrend.method,'_',today,'.rds', sep='')
  saveRDS(data_delfzijl, file=filename.output)

  #===============================================================================

  print('done processing the Delfzijl tide gauge data set')

}

#===============================================================================
# End
#===============================================================================

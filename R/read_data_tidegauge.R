#===============================================================================
# read tide gauge data
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================


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


#=====================
# function to trim temperature forcing to fit TG record unique years
# there might be missing years in TG record, so need to match each year and not
# just plop down an evenly spaced sequence
trimmed_forcing <- function(year_tidegauge, year_temperature, temperature) {
  output <- vector('list', 2); names(output) <- c('time','temperature')
  # check the beginning
  if(year_temperature[1] > year_tidegauge[1]) {print('ERROR - tide gauge record starts before temperature; add support for this situation')}
  # check the end
  if(max(year_temperature) < max(year_tidegauge)) {print('ERROR - tide gauge record ends after temperature; add support for this situation')}
  # match the indices of year_tidegauge within year_temperature
  imatch <- match(year_tidegauge, year_temperature)
  output$time <- year_temperature[imatch]
  output$temperature <- temperature[imatch]
  return(output)
}
#=====================


dat.dir <- '~/codes/EVT/data/tide_gauge_Europe/Hook_of_Holland_Oddo_data/'
data <- read_data(dat.dir=dat.dir, filetype='txt', septype='\t')

year  <- as.numeric(substr(as.character(data$date), start=1, stop=4))
month <- as.numeric(substr(as.character(data$date), start=5, stop=6))
day   <- as.numeric(substr(as.character(data$date), start=7, stop=8))

year.unique <- unique(year)
year.unique <- year.unique[order(year.unique)]
lsl.mean <- rep(NA, length(year.unique))
lsl.max <- rep(NA, length(year.unique))
data$sl.norm <- data$sl

# TODO
# TODO - fix for data gaps or incomplete years (can check first/last points in the year, or maxikappam difference between two points)
# TODO
# will want to convert all dates/times into a standard format (number of days
# after 1880?)

for (t in 1:length(year.unique)) {
  ind.this.year <- which(year==year.unique[t])
  lsl.mean[t] <- mean(data$sl[ind.this.year])
  data$sl.norm[ind.this.year] <- data$sl[ind.this.year] - lsl.mean[t]
  lsl.max[t] <- max(data$sl.norm[ind.this.year])
}

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



#===============================================================================
# End
#===============================================================================

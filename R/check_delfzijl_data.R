#===============================================================================
# check_delfzijl_data.R
#
# Verify that the data set from Alexander, Rijkwaterstaat matches the data set
# from Perry, PSMSL
#===============================================================================

library(date)

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


#===============================================================================
# data set from perry, psmsl

dat.dir <- '../data/Delfzijl_Oddo_data/'
data1 <- read_data(dat.dir=dat.dir, filetype='txt', septype='\t')

# convert sea levels from cm to mm, consistent with the other data
data1$sl <- 10* data1$sl

data1$year   <- as.numeric(substr(as.character(data1$date), start=1, stop=4))
data1$month  <- as.numeric(substr(as.character(data1$date), start=5, stop=6))
data1$day    <- as.numeric(substr(as.character(data1$date), start=7, stop=8))
data1$hour   <- as.numeric(substr(data1$time, 1,2))
data1$minute <- as.numeric(substr(data1$time, 4,5))

# time in days since 01 January 1960
data1$time.days <- as.numeric(mdy.date(month=data1$month, day=data1$day, year=data1$year)) + data1$hour/24 + data1$minute/(24*60)

# difference between time stamps (in units of days)
time.diff <- diff(data1$time.days)

# check that they are in the proper order, ascending
print(paste('Are there any times out of order? ',any(time.diff < 0), sep=''))

# put everything back in order - make sure you do this for the sea levels and
# other fields too, so do as a matrix. also recalculate the time differences,
# which you will need for averaging
data1 <- data1[order(data1$time.days),]
#===============================================================================



#===============================================================================
# data set from alexander, rijkwaterstaat

data2 <- read.table('../data/id1-DELFZL-187901010000-201701012359.txt', header = TRUE, sep=';', skip=4)
names(data2)[6] <- 'sl'
data2$sl <- data2$sl*10 # convert to mm from cm
#===============================================================================



#===============================================================================
sl_diff <- data2$sl[1:length(data1$sl)] - data1$sl
idiff <- which(sl_diff != 0)

# the rijkwaterstaat data is the same, except that the 8-hour gap is filled
# we just leave the missing data out.

#===============================================================================
# end
#===============================================================================


#

#===============================================================================
# processing of Balboa, Panama data.
#
# leads to the list object 'data_balboa', analogous to the 'data_many' object
# that contains the many tide gauge stations from UHSLC database;
# has all of the necessary information to estimate (MLE) the PP/GPD parameters
# for each of the 4 candidate model structures.
#
# Note that this one is in the UHSLC database, so we assume 'processing_many_stations.R'
# has been previously run, and we just pluck off the appropriate station from
# the 'data_many' object that results from that script. Then we save it as its
# own RDS object, for calibration. This is for simplicity and consistency with
# the other two stations for calibration (Delfzijl and Norfolk).
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================

filename.saveprogress <- '../output/processing_balboa.RData'

print('starting to process Balboa tide gauge data')

#===
#=== read in previous processed 'data_many' object, and pluck off Balboa
#===

data_many <- readRDS('../data/tidegauge_processed_manystations_26Jul2017.rds')

data_balboa <- data_many$rqh0302

#
#===============================================================================
# subsample smaller sets of the POT/GPD data
#===============================================================================
#

# initialize, and create list elements for these GPD experiments. then later
# remove from each sub-list item the years the experiment will not use
gpd.experiments <- c('gpd30','gpd50','gpd70','gpd90','gpd107')
years.gpd.experiments <- c(30,50,70,90,107); names(years.gpd.experiments) <- gpd.experiments
for (gpd.exp in gpd.experiments) {
  data_balboa[[gpd.exp]] <- data_balboa$gpd
  ind.experiment <- (length(data_balboa$gpd$year)-years.gpd.experiments[[gpd.exp]]+1):length(data_balboa$gpd$year)
  data_balboa[[gpd.exp]]$counts <- data_balboa[[gpd.exp]]$counts[ind.experiment]
  data_balboa[[gpd.exp]]$time_length <- data_balboa[[gpd.exp]]$time_length[ind.experiment]
  data_balboa[[gpd.exp]]$excesses <- data_balboa[[gpd.exp]]$excesses[ind.experiment]
  data_balboa[[gpd.exp]]$year <- data_balboa$gpd$year[ind.experiment]
  data_balboa[[gpd.exp]]$counts_all <- NULL
  data_balboa[[gpd.exp]]$time_length_all <- NULL
  data_balboa[[gpd.exp]]$excesses_all <- NULL
}

# that doesn't take as long... but lets be consistent I guess and save it.
save.image(file=filename.saveprogress)

# save final 'data_balboa' object to RDS to use later
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.output <- paste('../data/tidegauge_processed_balboa_',today,'.rds', sep='')
saveRDS(data_balboa, file=filename.output)

#===============================================================================

print('done processing the Balboa, Panama tide gauge data set')

#===============================================================================
# End
#===============================================================================

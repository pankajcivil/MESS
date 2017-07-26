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

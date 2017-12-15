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

processing_balboa <- function(dt.decluster, detrend.method, pot.threshold) {

filename.saveprogress <- '../output/processing_balboa.RData'

print('starting to process Balboa tide gauge data')

#===
#=== read in previous processed 'data_many' object, and pluck off Balboa
#===

data_many <- readRDS('../data/tidegauge_processed_manystations_decl3-pot99-annual_10Dec2017.rds')

data_balboa <- data_many$rqh0302
data_balboa$dt.decluster <- dt.decluster

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
filename.output <- paste('../data/tidegauge_processed_balboa_decl',data_balboa$dt.decluster,'-pot',pot.threshold*100,'-',detrend.method,'_',today,'.rds', sep='')
saveRDS(data_balboa, file=filename.output)

#===============================================================================

print('done processing the Balboa, Panama tide gauge data set')

}

#===============================================================================
# End
#===============================================================================

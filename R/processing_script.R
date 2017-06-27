#===============================================================================
# processing_script.R
#
# Process tide gauge data for three experiments:
# 1. GEV and Naveau-(i) (aka Papastathopoulos and Tawn iii) to annual block maxima
# 2. GEV and Naveau-(i) to monthly block maxima (preprocessing needed)
# 3. PP-GPD and Naveau-(i) to POT daily maxima (preprocessing needed for Naveau,
#    and included/not included in two separate experiments for POT/GPD)
#    (POT = Peaks Over Thresholds)
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================

# whoops, extRemes package has a function for this, and it works better for
# the EVT analysis. ah well.
decluster_timeseries <- function(time, year, time.series, min.dt) {
  decluster <- vector('list',3)
  names(decluster) <- c('time','year','time.series')
  tdiff <- diff(time)
  ind.too.close <- which(tdiff <= min.dt)
  if(length(ind.too.close) > 0) {
    # go through the places where there are time.series values too close together,
    # and set to throw the smaller value away.
    # need to account for the possibility that there are more than two values in
    # a 'cluster'.
    # indices where new clusters begin: (tack on first one manually)
    ind.clusters <- c(ind.too.close[1] , ind.too.close[which(diff(ind.too.close) > 1) + 1])
    if(length(ind.clusters) > 1) {
      # case where there are multiple clusters to consider
      # initialize by fixing to remove all of the spots that are too close. then
      # we will remove the indices of each cluster's maximum
      irem <- unique(c(ind.too.close, ind.too.close + 1))
      for (i in 1:length(ind.clusters)) {
        # how long is this cluster?
        if(i < length(ind.clusters)) {
          first <- ind.clusters[i]
          last  <- ind.too.close[which(ind.too.close==ind.clusters[i+1]) -1]+1
        } else {
          first <- ind.clusters[i]
          if (length(which(ind.too.close > first))==0) {last <- first + 1
          } else {last <- first+length(which(ind.too.close > first)) + 1}
        }
        ind.this.cluster <- first:last
        isave <- ind.this.cluster[which.max(time.series[first:last])]
        # remove this cluster's maximum from irem; they might be out of order
        isave <- which(irem==isave)
        irem <- irem[-isave]
      }
    } else {
      # case with only one cluster
      # initialize by fixing to remove all of the cluster. then we will remove
      # the index of the cluster maximum from this
      irem <- unique(c(ind.too.close, ind.too.close + 1))
      irem <- irem[-which.max(time.series[irem])]
    }
    decluster$time <- time[-irem]
    decluster$year <- year[-irem]
    decluster$time.series <- time.series[-irem]
  } else {
    decluster$time <- time
    decluster$year <- year
    decluster$time.series <- time.series
  }
  return(decluster)
}

#
#===============================================================================
# Read tide gauge data for Delfzijl, The Netherlands (site that is the focus of
# the flood risk analysis) and a variety of European tide gauge stations that
# are nearby.
#===============================================================================
#

# Delfzijl
source('processing_delfzijl.R')

# TODO

# Other European stations
source('processing_europe.R')

#
#===============================================================================
# End
#===============================================================================
#

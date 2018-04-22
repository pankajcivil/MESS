#===============================================================================
# This file reads in the output files from the marginal likelihood estimator
# for each model and combines them into a combined .RDS file with bma weights
# and marginal likelihoods.
#
# Requires the RData files from the 'bridge_sample.R' routine as input
#
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

library(Hmisc)

appen <- '_threshold997'
path.ml <- paste('/home/scrim/axw322/codes/EVT/output/bma',appen,sep='')
path.out <- '/home/scrim/axw322/codes/EVT/output'

filename.likelihood <- paste('log_marginal_likelihood',appen,'.rds',sep='')
filename.weights <- paste('bma_weights',appen,'.rds',sep='')

types.of.priors <- 'normalgamma'

site.names <- c('Delfzijl', 'Norfolk')
n_sites <- length(site.names)
bma.weights <- vector('list', n_sites)
names(bma.weights) <- site.names

log.marg.lik <- vector('list', n_sites)
names(log.marg.lik) <- site.names


gpd.models <- c('gpd3','gpd4','gpd5','gpd6')
n_model <- length(gpd.models)


data.lengths <- vector('list', n_sites); names(data.lengths) <- site.names
data.lengths$Norfolk <- c('30','50','70','89')
data.lengths$Delfzijl <- c('30','50','70','90','110','137')

for (site in site.names) {

  bma.weights[[site]] <- vector('list', length(data.lengths[[site]]))
  names(bma.weights[[site]]) <- data.lengths[[site]]

  for (year in data.lengths[[site]]) {
    bma.weights[[site]][[year]] <- rep(NA, n_model)
    names(bma.weights[[site]][[year]]) <- gpd.models
  }
  log.marg.lik[[site]] <- vector('list', length(data.lengths[[site]]))
  names(log.marg.lik[[site]]) <- data.lengths[[site]]
  for (year in data.lengths[[site]]) {
      log.marg.lik[[site]][[year]] <- rep(NA, n_model)
      names(log.marg.lik[[site]][[year]]) <- gpd.models

  }
}

files <- list.files(path=path.ml, full.names=TRUE, recursive=FALSE)

for (file in files) {
  load(file)
  ##site <- capitalize(toString(station))
  station_name <- toString(station)
  site <- paste(toupper(substr(station_name, 1, 1)), substr(station_name, 2, nchar(station_name)), sep="")
  year <- unlist(strsplit(file, split="[_. ]"))[5]
#  data.case <- which.min(abs(as.numeric(levels(data.length)[data.length])-exp.years))]
  log.marg.lik[[site]][[year]][[gpd.model]] <- ml[length(ml)]
}

for (site in site.names) {
  for (year in data.lengths[[site]]) {
    ml <- log.marg.lik[[site]][[year]]
    ml.scale <- ml - max(ml,na.rm=TRUE)
    for (model in gpd.models) {
      if (!is.na(log.marg.lik[[site]][[year]][[model]])) {
        bma.weights[[site]][[year]][[model]] <- exp(ml.scale[model])/sum(exp(ml.scale), na.rm=TRUE)
      }
    }
  }
}

saveRDS(log.marg.lik, paste(path.out,filename.likelihood, sep="/"))
saveRDS(bma.weights, paste(path.out,filename.weights, sep="/"))

#===============================================================================
# End
#===============================================================================

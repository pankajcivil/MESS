###############################################################################
# This file reads in the output files from the marginal likelihood estimator  #
# for each model and combines them into a combined .RDS file with bma weights #
# and marginal likelihoods.                                                   #
###############################################################################
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

path.ml <- '/storage/home/vxs914/work/MESS/ml'
path.out <- '/storage/home/vxs914/work/MESS/output'

filename.out <- 'bma_weights.rds'

types.of.priors <- 'normalgamma'

bma.weights <- vector('list', 3)
names(bma.weights) <- site.names <- c('Delfzijl', 'Balboa', 'Norfolk')

log.marg.lik <- vector('list', 3)
names(log.marg.lik) <- site.names

data.length.names <- c('y30', 'y50', 'y70', 'y90', 'y110', 'y130')

years=c('30','50','70','89','90','107','110','137')

gpd.models <- c('gpd3','gpd4','gpd5','gpd6')

for (site in site.names) {
  
  bma.weights[[site]] <- vector('list', 8)
  names(bma.weights[[site]]) <- years
  for (year in years) {
    bma.weights[[site]][[year]] <- rep(NA, 4)
    names(bma.weights[[site]][[year]]) <- gpd.models
  }
  log.marg.lik[[site]] <- vector('list', 8)
  names(log.marg.lik[[site]]) <- years
  for (year in years) {
      log.marg.lik[[site]][[year]] <- rep(NA, 4)
      names(log.marg.lik[[site]][[year]]) <- gpd.models
  
  }
}

files <- list.files(path=path.ml, full.names=TRUE, recursive=FALSE)

for (file in files) {
  load(file)
  site <- capitalize(station)
  year <- unlist(strsplit(file, split="[_. ]"))[4]
#  data.case <- which.min(abs(as.numeric(levels(data.length)[data.length])-exp.years))]
  log.marg.lik[[site]][[year]][[gpd.model]] <- ml[length(ml)]
  
}

for (site in site.names) {
  for (year in years) {
    ml <- log.marg.lik[[site]][[year]]
    ml.scale <- ml - max(ml,na.rm=TRUE)
    for (model in gpd.models) {
      if (!is.na(log.marg.lik[[site]][[year]][[model]])) {
        bma.weights[[site]][[year]][[model]] <- exp(ml.scale[model])/sum(exp(ml.scale), na.rm=TRUE)
      }
    }
  }
}

saveRDS(log.marg.lik, "output/log_marginal_likelihood.rds")
saveRDS(bma.weights, "output/bma_weights.rds")

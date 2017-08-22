###############################################################################
# This file reads in the output files from the marginal likelihood estimator  #
# for each model and combines them into a combined .RDS file with bma weights #
# and marginal likelihoods.                                                   #
###############################################################################

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

gpd.models <- c('gpd3','gpd4','gpd5','gpd6')

for (site in site.names) {
  
  bma.weights[[site]] <- vector('list', 6)
  names(bma.weights[[site]]) <- data.length.names
  for (data.length in data.length.names) {
    bma.weights[[site]][[data.length]] <- rep(NA, 4)
    names(bma.weights[[site]][[data.length]]) <- gpd.models
  }
  log.marg.lik[[site]] <- vector('list', 6)
  names(log.marg.lik[[site]]) <- data.length.names
  for (data.length in data.length.names) {
      log.marg.lik[[site]][[data.length]] <- rep(NA, 4)
      names(log.marg.lik[[site]][[data.length]]) <- gpd.models
  
  }
}
  
exp.years <- c(30,50,70,90,110,130)

files <- list.files(path=path.ml, full.names=TRUE, recursive=FALSE)

for (file in files) {
  load(file)
  site <- capitalize(station)
  data.case <- data.length.names[which.min(abs(as.numeric(levels(data.length)[data.length])-exp.years))]
  log.marg.lik[[site]][[data.case]][[gpd.model]] <- ml[length(ml)]
}

for (site in site.names) {
  for (data.length in data.length.names) {
    ml <- log.marg.lik[[site]][[data.length]]
    ml.scale <- ml - max(ml,na.rm=TRUE)
    for (model in gpd.models) {
      if (!is.na(log.marg.lik[[site]][[data.length]][[model]])) {
        bma.weights[[site]][[data.length]][[model]] <- exp(ml.scale[model])/sum(exp(ml.scale), na.rm=TRUE)
      }
    }
  }
}

saveRDS(log.marg.lik, "output/log_marginal_likelihood.rds")
saveRDS(bma.weights, "output/bma_weights.rds")
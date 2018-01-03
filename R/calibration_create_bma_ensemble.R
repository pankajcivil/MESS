#===============================================================================
# calibration_create_bma_ensemble.R
#
# Requires RData files from calibration_dayPOT-experiments_driver.R for all
# sites and experiments.
# Draw parameters from the raw parameter results from the MCMC and cook up a
# BMA-weighted ensemble of return levels.
#
# A few things you'll get from the RData files that are nice to have:
#  nnode_mcmc = number of parallel chains
#
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

# mcmc results files
mcmc.balboa <- '../output/everything_mcmc_ppgpd-experiments_balboa_normalgamma_21Dec2017.RData'
mcmc.norfolk <- '../output/everything_mcmc_ppgpd-experiments_norfolk_normalgamma_21Dec2017.RData'
mcmc.delfzijl <- '../output/everything_mcmc_ppgpd-experiments_delfzijl_normalgamma_20Dec2017.RData'


# TODO -- would need to revise the above because only saving the raw MCMC
# output now.  that is, need to do the burn-in


# bma weight results object
bma_weights <- readRDS('../output/bma_weights.rds')

#===============================================================================
# Read and save the raw posterior parameter sets - all of them, not just the
# thinning versions.
#===============================================================================

mcmc.raw.results <- vector('list', nsites); names(mcmc.raw.results) <- site.names
load(mcmc.balboa)
mcmc.raw.results$Balboa <- chains_burned
load(mcmc.norfolk)
mcmc.raw.results$Norfolk <- chains_burned
load(mcmc.delfzijl)
mcmc.raw.results$Delfzijl <- chains_burned

# "unwrap" the list object 'chains_burned' so that you can sample out of it

mcmc.all.parameters <- list.init

for (site in site.names) {
  print(paste('pooling all ',nnode_mcmc,' parallel chains for ',site,'...',sep=''))
  for (dd in 1:all.data[[site]]) {
    print(paste('  experiment ',data.experiment.names[dd],'...',sep=''))
    for (model in types.of.gpd) {
      for (ind.chain in 1:nnode_mcmc) {
        mcmc.all.parameters[[site]][[dd]][[model]] <- rbind(mcmc.all.parameters[[site]][[dd]][[model]], mcmc.raw.results[[site]][[dd]][[model]][[ind.chain]])
      }
    }
  }
}
save.image(file=filename.saveprogress)

#===============================================================================
# Draw 'n.ensemble' sets of parameters out of each of the models for the BMA
# ensemble. Number chosen for consistency with the other ensembles.
#===============================================================================

mcmc.sampled.parameters <- list.init
for (site in site.names) {
  for (dd in 1:all.data[[site]]) {
    for (model in types.of.gpd) {
      ind.bma.sample <- sample(1:nrow(mcmc.all.parameters[[site]][[dd]][[model]]), size=n.ensemble, replace=FALSE)
      mcmc.sampled.parameters[[site]][[dd]][[model]] <- mcmc.all.parameters[[site]][[dd]][[model]][ind.bma.sample,]
    }
  }
}

#===============================================================================
# Calculate return levels for each of the candidate models.
# Then add them up by BMA weights
#===============================================================================

rl100.bmamodels <- list.init
rl100.bma <- list.init
for (site in site.names) {
  for (dd in 1:all.data[[site]]) {
    rl100.bma[[site]][[dd]] <- vector('list', nyears); names(rl100.bma[[site]][[dd]]) <- year.names
    for (model in types.of.gpd) {
      rl100.bmamodels[[site]][[dd]][[model]] <- vector('list', nyears); names(rl100[[site]][[dd]][[model]]) <- year.names
      for (year in year.names) {
        rl100.bmamodels[[site]][[dd]][[model]][[year]] <- rep(NA, n.ensemble)
        print(paste('calculating return levels for BMA ensemble...',site,dd,model,year,sep=' - '))
        tbeg <- proc.time()
        pb <- txtProgressBar(min=0,max=n.ensemble, initial=0,style=3)
        for (sow in 1:n.ensemble) {
          if (length(parnames[[model]])==3) {
            lambda <- mcmc.sampled.parameters[[site]][[dd]][[model]][sow,match('lambda', parnames[[model]])]
            sigma <- mcmc.sampled.parameters[[site]][[dd]][[model]][sow,match('sigma', parnames[[model]])]
            xi <- mcmc.sampled.parameters[[site]][[dd]][[model]][sow,match('xi', parnames[[model]])]
          } else if(length(parnames[[model]])==4) {
            lambda <- mcmc.sampled.parameters[[site]][[dd]][[model]][sow,match('lambda0', parnames[[model]])] +
                      mcmc.sampled.parameters[[site]][[dd]][[model]][sow,match('lambda1', parnames[[model]])]*temperature.years[[year]]
            sigma <- mcmc.sampled.parameters[[site]][[dd]][[model]][sow,match('sigma', parnames[[model]])]
            xi <- mcmc.sampled.parameters[[site]][[dd]][[model]][sow,match('xi', parnames[[model]])]
          } else if(length(parnames[[model]])==5) {
            lambda <- mcmc.sampled.parameters[[site]][[dd]][[model]][sow,match('lambda0', parnames[[model]])] +
                      mcmc.sampled.parameters[[site]][[dd]][[model]][sow,match('lambda1', parnames[[model]])]*temperature.years[[year]]
            sigma <- exp(mcmc.sampled.parameters[[site]][[dd]][[model]][sow,match('sigma0', parnames[[model]])] +
                         mcmc.sampled.parameters[[site]][[dd]][[model]][sow,match('sigma1', parnames[[model]])]*temperature.years[[year]])
            xi <- mcmc.sampled.parameters[[site]][[dd]][[model]][sow,match('xi', parnames[[model]])]
          } else if(length(parnames[[model]])==6) {
            lambda <- mcmc.sampled.parameters[[site]][[dd]][[model]][sow,match('lambda0', parnames[[model]])] +
                      mcmc.sampled.parameters[[site]][[dd]][[model]][sow,match('lambda1', parnames[[model]])]*temperature.years[[year]]
            sigma <- exp(mcmc.sampled.parameters[[site]][[dd]][[model]][sow,match('sigma0', parnames[[model]])] +
                         mcmc.sampled.parameters[[site]][[dd]][[model]][sow,match('sigma1', parnames[[model]])]*temperature.years[[year]])
            xi <- mcmc.sampled.parameters[[site]][[dd]][[model]][sow,match('xi0', parnames[[model]])] +
                  mcmc.sampled.parameters[[site]][[dd]][[model]][sow,match('xi1', parnames[[model]])]*temperature.years[[year]]
          }
          rl100.bmamodels[[site]][[dd]][[model]][[year]][sow] <- rlevd(100, scale=sigma, shape=xi,
                                                                       threshold=data.sites[[site]]$gpd$threshold,
                                                                       type='GP',
                                                                       npy=365.25,
                                                                       rate=lambda)
          setTxtProgressBar(pb, sow)
        }
        close(pb)
      }
    }
    # add up the BMA ensemble from the individual model results calculated above
    # WARNING! - HARD-CODED FOR THE SET OF GPD MODELS USED HERE! - WARNING!
    for (year in year.names) {
      rl100.bma[[site]][[dd]][[year]] <- bma_weights[[site]][[dd]]['gpd3'] * rl100.bmamodels[[site]][[dd]]$gpd3[[year]] +
                                         bma_weights[[site]][[dd]]['gpd4'] * rl100.bmamodels[[site]][[dd]]$gpd4[[year]] +
                                         bma_weights[[site]][[dd]]['gpd5'] * rl100.bmamodels[[site]][[dd]]$gpd5[[year]] +
                                         bma_weights[[site]][[dd]]['gpd6'] * rl100.bmamodels[[site]][[dd]]$gpd6[[year]]
    }
  }
}
save.image(file=filename.saveprogress)

#===============================================================================
# end
#===============================================================================

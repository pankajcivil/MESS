#===============================================================================
# analysis_and_plotting.R
#
# Analysis and plotting for wong, klufas and keller pp/gpd analysis.
# This assumes you've downloaded the Github repo and are running from the
# 'R' directory.
#
# If any this is vague/confusing and you have questions, email me. Seriously.
#
# Questions? Tony Wong (anthony.e.wong@colorado.edu)
#===============================================================================
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

rm(list=ls())

# Some preliminaries

setwd('~/codes/EVT/R')

library(ncdf4)
library(extRemes)
library(Bolstad)

# get some colorblind-friendly colors to plot with
source('colorblindPalette.R')

# set useful directories -- assumes you are in the 'R' directory within the repo
# directory structure
plot.dir <- '../figures/'
output.dir <- '../output/'

# calibrated parameter sets (samples; all should be same size)
# these are the results from 'calibration_dayPOT-experiments_driver.R'
filename.norfolk.normalgamma <-  paste(output.dir,'calibratedParameters_ppgpd-experiments_norfolk_normalgamma_decl3-pot99_31Mar2018.nc', sep='')
filename.delfzijl.normalgamma <- paste(output.dir,'calibratedParameters_ppgpd-experiments_delfzijl_normalgamma_decl3-pot99_31Mar2018.nc',sep='')
filename.norfolk.uniform <-      paste(output.dir,'calibratedParameters_ppgpd-experiments_norfolk_uniform_decl3-pot99_31Mar2018.nc',     sep='')
filename.delfzijl.uniform <-     paste(output.dir,'calibratedParameters_ppgpd-experiments_delfzijl_uniform_decl3-pot99_31Mar2018.nc',    sep='')

filename.norfolk.data <- '../data/tidegauge_processed_norfolk_decl3-pot99-annual_06Dec2017.rds'
filename.delfzijl.data <- '../data/tidegauge_processed_delfzijl_decl3-pot99-annual_20Dec2017.rds'

filename.priors <- '../output/surge_priors_normalgamma_ppgpd_decl3-pot99_28Mar2018.rds'

# file to save progress as you run
filename.saveprogress <- '../output/analysis_inprogress.RData'

#===============================================================================



#===============================================================================
#===============================================================================
# ANALYSIS
#===============================================================================
#===============================================================================




#===============================================================================
# read NAO data for projecting surge levels covarying with NAO index (DJF)
#===============================================================================

# yields 'nao_forc' and 'time_forc', between 1850 and 2100
# both historical and future NAO index standardized relative to 2001-2016 mean/stdev

source('read_data_naoindex.R')

#===============================================================================
# read parameter sets
#===============================================================================

# create list objects to store the parameters
site.names <- c('Delfzijl','Norfolk'); nsites <- length(site.names)
data.lengths <- c(30,50,70,90,110,130); ndata.exp <- length(data.lengths)
data.experiment.names <- rep(NA, length(data.lengths)); for (dd in 1:length(data.lengths)) {data.experiment.names[dd] <- paste('y',data.lengths[dd],sep='')}
types.of.gpd <- c('gpd3','gpd4','gpd5','gpd6'); nmodel <- length(types.of.gpd)

# define a generic list object that will have values for (i) each site,
# (ii) each data length experiment and (iii) each GPD model structure. (in that
# order)
list.init <- vector('list', nsites); names(list.init) <- site.names
for (site in site.names) {
  list.init[[site]] <- vector('list', length(data.lengths)); names(list.init[[site]]) <- data.experiment.names
  for (data.exp in data.experiment.names) {
    list.init[[site]][[data.exp]] <- vector('list', nmodel); names(list.init[[site]][[data.exp]]) <- types.of.gpd
  }
}

# initialize anything you would like to be defined for each site, for each data
# length and for each GPD model structure
gpd.parameters <- list.init
rl100 <- list.init

# index to store within the list level [[data.experiments]] which of the
# experiments corresponds to using all of the data. so you would use:
#    gpd.parameters$delfzijl[[all.data]]$gpd3
# for example to get at the calibrated gpd3 (stationary model) parameters for
# Delfzijl.
all.data <- rep(NA, nsites); names(all.data) <- site.names
data.lengths <- vector('list', nsites); names(data.lengths) <- site.names
data.sites <- vector('list', nsites); names(data.sites) <- site.names

# hard-coding here, even though it is bad practice. the different sites ahve
# different lengths of record, making things a bit complicated.

# parameter nmaes
parnames <- vector('list', nmodel); names(parnames) <- types.of.gpd

for (site in site.names) {
  if(site=='Norfolk') {
    ncdata <- nc_open(filename.norfolk.normalgamma)
    for (model in types.of.gpd) {
      gpd.parameters[[site]]$y30[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd30.',model,sep='')))
      gpd.parameters[[site]]$y50[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd50.',model,sep='')))
      gpd.parameters[[site]]$y70[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd70.',model,sep='')))
      gpd.parameters[[site]]$y90[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd89.',model,sep='')))
      parnames[[model]] <- ncvar_get(ncdata, paste('parnames.',model,sep=''))
    }
    n.ensemble <- nrow(gpd.parameters[[site]]$y30$gpd3)
    data.sites[[site]] <- readRDS(filename.norfolk.data)
    data.lengths[[site]] <- names(data.sites[[site]])[intersect(which(nchar(names(data.sites[[site]]))>3) , grep('gpd', names(data.sites[[site]])))]
    all.data[[site]] <- length(data.lengths[[site]])
  } else if(site=='Delfzijl') {ncdata <- nc_open(filename.delfzijl.normalgamma)
    for (model in types.of.gpd) {
      gpd.parameters[[site]]$y30[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd30.',model,sep='')))
      gpd.parameters[[site]]$y50[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd50.',model,sep='')))
      gpd.parameters[[site]]$y70[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd70.',model,sep='')))
      gpd.parameters[[site]]$y90[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd90.',model,sep='')))
      gpd.parameters[[site]]$y110[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd110.',model,sep='')))
      gpd.parameters[[site]]$y130[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd137.',model,sep='')))
    }
    n.ensemble <- nrow(gpd.parameters[[site]]$y30$gpd3)
    data.sites[[site]] <- readRDS(filename.delfzijl.data)
    data.lengths[[site]] <- names(data.sites[[site]])[intersect(which(nchar(names(data.sites[[site]]))>3) , grep('gpd', names(data.sites[[site]])))]
    all.data[[site]] <- length(data.lengths[[site]])
    if(nrow(gpd.parameters[[site]]$y30$gpd3) != n.ensemble) {print('ERROR - all sites ensembles must be the same size')}
  } else {print('ERROR - unrecognized site name')}
}

#===============================================================================
# calculate return levels
#===============================================================================

# years to grab the return levels
rl.years <- c(2000, 2016, 2065); nyears <- length(rl.years)
year.names <- rep(NA, length(rl.years)); for (year in 1:length(rl.years)) {year.names[year] <- paste('y',rl.years[year],sep='')}

# take the 11-year mean, centered at the year in question
nao.years <- NULL
for (year in rl.years) {
  nao.years <- c(nao.years, mean(nao_forc[which(abs(time_forc-year)<=5)]))
}
names(nao.years) <- year.names

# this is definitely coded like a neanderthal. can come back and code it more
# efficiently..?
for (site in site.names) {
  for (data.len in data.experiment.names[1:all.data[[site]]]) {
    for (model in types.of.gpd) {
      rl100[[site]][[data.len]][[model]] <- vector('list', nyears); names(rl100[[site]][[data.len]][[model]]) <- year.names
      for (year in year.names) {
        rl100[[site]][[data.len]][[model]][[year]] <- rep(NA, n.ensemble)
        print(paste('calculating return levels...',site,data.len,model,year,sep=' - '))
        tbeg <- proc.time()
        pb <- txtProgressBar(min=0,max=n.ensemble, initial=0,style=3)
        for (sow in 1:n.ensemble) {
          if (length(parnames[[model]])==3) {
            lambda <- gpd.parameters[[site]][[data.len]][[model]][sow,match('lambda', parnames[[model]])]
            sigma <- gpd.parameters[[site]][[data.len]][[model]][sow,match('sigma', parnames[[model]])]
            xi <- gpd.parameters[[site]][[data.len]][[model]][sow,match('xi', parnames[[model]])]
          } else if(length(parnames[[model]])==4) {
            lambda <- gpd.parameters[[site]][[data.len]][[model]][sow,match('lambda0', parnames[[model]])] +
                      gpd.parameters[[site]][[data.len]][[model]][sow,match('lambda1', parnames[[model]])]*nao.years[[year]]
            sigma <- gpd.parameters[[site]][[data.len]][[model]][sow,match('sigma', parnames[[model]])]
            xi <- gpd.parameters[[site]][[data.len]][[model]][sow,match('xi', parnames[[model]])]
          } else if(length(parnames[[model]])==5) {
            lambda <- gpd.parameters[[site]][[data.len]][[model]][sow,match('lambda0', parnames[[model]])] +
                      gpd.parameters[[site]][[data.len]][[model]][sow,match('lambda1', parnames[[model]])]*nao.years[[year]]
            sigma <- exp(gpd.parameters[[site]][[data.len]][[model]][sow,match('sigma0', parnames[[model]])] +
                         gpd.parameters[[site]][[data.len]][[model]][sow,match('sigma1', parnames[[model]])]*nao.years[[year]])
            xi <- gpd.parameters[[site]][[data.len]][[model]][sow,match('xi', parnames[[model]])]
          } else if(length(parnames[[model]])==6) {
            lambda <- gpd.parameters[[site]][[data.len]][[model]][sow,match('lambda0', parnames[[model]])] +
                      gpd.parameters[[site]][[data.len]][[model]][sow,match('lambda1', parnames[[model]])]*nao.years[[year]]
            sigma <- exp(gpd.parameters[[site]][[data.len]][[model]][sow,match('sigma0', parnames[[model]])] +
                         gpd.parameters[[site]][[data.len]][[model]][sow,match('sigma1', parnames[[model]])]*nao.years[[year]])
            xi <- gpd.parameters[[site]][[data.len]][[model]][sow,match('xi0', parnames[[model]])] +
                  gpd.parameters[[site]][[data.len]][[model]][sow,match('xi1', parnames[[model]])]*nao.years[[year]]
          }
          rl100[[site]][[data.len]][[model]][[year]][sow] <- rlevd(100, scale=sigma, shape=xi,
                                                                   threshold=data.sites[[site]]$gpd$threshold,
                                                                   type='GP',
                                                                   npy=365.25,
                                                                   rate=lambda)
          setTxtProgressBar(pb, sow)
        }
        close(pb)
      }
    }
  }
}

#===============================================================================
# Bayesian model averaging ensemble and model comparison metrics
#===============================================================================

# need the likelihood functions and the priors
source('likelihood_ppgpd.R')
priors <- readRDS(filename.priors)

llik <- list.init
lpri <- list.init
lpost <- list.init
llik.mod <- list.init
lpost.mod <- list.init

for (site in site.names) {
  for (ind.data in 1:all.data[[site]]) {
    for (model in types.of.gpd) {
      print(paste('calculating log-posts/-likelihoods/-priors ',site,data.experiment.names[ind.data],model, sep=' - '))
      if(model=='gpd3') {auxiliary <- NULL
      } else {auxiliary <- trimmed_forcing(data.sites[[site]][[data.lengths[[site]][ind.data]]]$year, time_forc, nao_forc)$forcing}
      ##lpri[[site]][[ind.data]][[model]] <- sapply(1:n.ensemble, function(sow) {log_prior_ppgpd(parameters=gpd.parameters[[site]][[ind.data]][[model]][sow,], parnames=parnames[[model]], priors=priors, model=model)})
      llik[[site]][[ind.data]][[model]] <- sapply(1:n.ensemble, function(sow) {log_like_ppgpd(parameters=gpd.parameters[[site]][[ind.data]][[model]][sow,], parnames=parnames[[model]], data_calib=data.sites[[site]][[data.lengths[[site]][ind.data]]], auxiliary=auxiliary)})
      ##lpost[[site]][[ind.data]][[model]] <- lpri[[site]][[ind.data]][[model]] + llik[[site]][[ind.data]][[model]]
      llik.mod[[site]][[ind.data]][[model]] <- mean(llik[[site]][[ind.data]][[model]][is.finite(llik[[site]][[ind.data]][[model]])])
      ##lpost.mod[[site]][[ind.data]][[model]] <- mean(lpost[[site]][[ind.data]][[model]][is.finite(lpost[[site]][[ind.data]][[model]])])
    }
  }
}
save.image(file=filename.saveprogress)

#===============================================================================
# Calculate AIC, BIC, DIC for each ensemble
# (of course, this could be wrapped up in the loop structure above, but it is
# nice to compartmentalize the code)
#===============================================================================

# new initialization object for the model comparison metrics, so we have a table
# comparing each of the models instead of a list for each of them
metric.init <- vector('list', nsites); names(metric.init) <- site.names
for (site in site.names) {
  metric.init[[site]] <- mat.or.vec(all.data[[site]], nmodel)
  rownames(metric.init[[site]]) <- data.experiment.names[1:all.data[[site]]]
  colnames(metric.init[[site]]) <- types.of.gpd
}

aic <- metric.init
bic <- metric.init
dic <- metric.init
dev <- list.init      # deviance
par.mean <- list.init # expected value of the model parameters
dev.mean <- list.init # deviance at expected value of the parameters
np.eff <- list.init   # effective number of parameters (p_D)

for (site in site.names) {
  for (ind.data in 1:all.data[[site]]) {
    ndata <- sum(data.sites[[site]][[data.lengths[[site]][ind.data]]]$counts)
    for (model in types.of.gpd) {
      imax <- which.max(llik[[site]][[ind.data]][[model]])
      aic[[site]][ind.data, model] <- -2*llik[[site]][[ind.data]][[model]][imax] + length(parnames[[model]])*2
      bic[[site]][ind.data, model] <- -2*llik[[site]][[ind.data]][[model]][imax] + length(parnames[[model]])*log(ndata)
      dev[[site]][[ind.data]][[model]] <- -2*llik[[site]][[ind.data]][[model]]
      if(model=='gpd3') {auxiliary <- NULL
      } else {auxiliary <- trimmed_forcing(data.sites[[site]][[data.lengths[[site]][ind.data]]]$year, time_forc, nao_forc)$forcing}
      par.mean[[site]][[ind.data]][[model]] <- apply(gpd.parameters[[site]][[ind.data]][[model]], MARGIN=2, FUN=mean)
      dev.mean[[site]][[ind.data]][[model]] <- log_like_ppgpd(parameters=par.mean[[site]][[ind.data]][[model]], parnames=parnames[[model]], data_calib=data.sites[[site]][[data.lengths[[site]][ind.data]]], auxiliary=auxiliary)
      np.eff[[site]][[ind.data]][[model]] <- mean(dev[[site]][[ind.data]][[model]]) - dev.mean[[site]][[ind.data]][[model]]
      dic[[site]][ind.data, model] <- np.eff[[site]][[ind.data]][[model]] + mean(dev[[site]][[ind.data]][[model]])
    }
  }
}

# bma weight results object
bma_weights <- readRDS('../output/bma_weights_threshold99.rds')

# create table 1 for paper
comparison_table1 <- cbind( c(rep('delfzijl', nmodel), rep('norfolk', nmodel)),
                            c(aic$Delfzijl['y130',], aic$Norfolk['y90',]) ,
                            c(bic$Delfzijl['y130',], bic$Norfolk['y90',]) ,
                            c(dic$Delfzijl['y130',], dic$Norfolk['y90',]) ,
                            c(bma_weights$Delfzijl$`137`, bma_weights$Norfolk$`89`))
colnames(comparison_table1) <- c('site','aic','bic','dic','bma')
write.csv(x=comparison_table1, file='../output/model_comparisons_threshold99.csv')

# map from Vivek's BMA weights object to match the data.lengths one here
#bma_weights$Delfzijl <- bma_weights$Delfzijl[-6]   # get rid of `107` experiment
#bma_weights$Delfzijl <- bma_weights$Delfzijl[-4]   # get rid of `89` experiment


# calculate BMA-weighted ensemble using the ensemble members already drawn
# have also checked this making larger (up to 1e6) and separate samples for each
# of the 4 models, and the quantiles do not change by more than a millimeter
if(TRUE) {
  rl100.bma <- list.init
  for (site in site.names) {
    for (dd in 1:length(data.lengths[[site]])) {
      for (year in year.names) {
        rl100.bma[[site]][[dd]][[year]] <- bma_weights[[site]][[dd]]['gpd3']*rl100[[site]][[dd]]$gpd3[[year]] +
                                           bma_weights[[site]][[dd]]['gpd4']*rl100[[site]][[dd]]$gpd4[[year]] +
                                           bma_weights[[site]][[dd]]['gpd5']*rl100[[site]][[dd]]$gpd5[[year]] +
                                           bma_weights[[site]][[dd]]['gpd6']*rl100[[site]][[dd]]$gpd6[[year]]
      }
    }
  }
} else {
  # alternatively, calculate BMA-weighted ensemble using the full MCMC reuslts
  #source('calibration_create_bma_ensemble.R')
}

# save progress
save.image(file=filename.saveprogress)

#===============================================================================




#===============================================================================
#===============================================================================
# PLOTTING
#===============================================================================
#===============================================================================



#
#===============================================================================
# FIGURE 1 - Motivation figure, showing 30-year blocks and the estimated
#            distributions (box-whisker/box-lighter-box) for each block, for
#            each site.

site <- 'Delfzijl'
load(paste('../output/sensitivity_returnlevels_mcmc_',site,'_20Dec2017.RData', sep=''))
returnlevel <- readRDS(paste('../output/sensitivity_returnlevels_',site,'_20Dec2017.rds', sep=''))
nblocks <- length(returnlevel)

## Calculate the quantiles to plot
quantiles.to.grab <- c(.05, .25, .5, .75, .95)
quantile.names <- rep(NULL, length(quantiles.to.grab))
for (qq in 1:length(quantiles.to.grab)) {
  if(quantiles.to.grab[qq] >= .10) {
    quantile.names[qq] <- paste('q',100*quantiles.to.grab[qq], sep='')
  } else if(quantiles.to.grab[qq] < .10 & quantiles.to.grab[qq] >= 0) {
    quantile.names[qq] <- paste('q0',100*quantiles.to.grab[qq], sep='')
  }
}
returnlevel.quantiles <- mat.or.vec(nblocks, length(quantiles.to.grab))
colnames(returnlevel.quantiles) <- quantile.names
for (bb in 1:nblocks) {
  # the /1000 is to convert to m from mm
  returnlevel.quantiles[bb,] <- quantile(returnlevel[[bb]], quantiles.to.grab)/1000
}

# useful analysis:
widest <- max( returnlevel.quantiles[,'q95'] - returnlevel.quantiles[,'q05'] )
narrowest <- min( returnlevel.quantiles[,'q95'] - returnlevel.quantiles[,'q05'] )
highest <- max( returnlevel.quantiles[,'q50'])
lowest <- min( returnlevel.quantiles[,'q50'])

## Useful for plotting - centers of the time blocks used in the experiments
block.years.center <- apply(X=block.years, MARGIN=1, FUN=median)

block.colors <- colorRampPalette(c("darkslateblue","royalblue","turquoise1"),space="Lab")(max(nblocks))
block.colors.lighter <- paste(block.colors, "70", sep="")


## The actual figure

pdf(paste(plot.dir,'stormsurge_sensitivity_boxwhisker_',site,'.pdf',sep=''),width=4,height=3,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(.8,.7,.15,.2))
halfwidth <- 2 # half the width of the boxes, in years
# put the first median bar down, to get hte plot started
plot(c(block.years.center[1]-halfwidth, block.years.center[1]+halfwidth), rep(returnlevel.quantiles[1,'q50'],2),
     type='l', lwd=3, col='black', xlim=c(1900,2004), ylim=c(0,8), xlab='', ylab='', las=1)
# now add the darker 25-75% range polygon before the median bars, ...
for (bb in 1:nblocks) {
    times.beg.end <- c(block.years.center[bb]-halfwidth, block.years.center[bb]+halfwidth)
    polygon(c(times.beg.end,rev(times.beg.end)), c(returnlevel.quantiles[bb,c('q25','q25')],rev(returnlevel.quantiles[bb,c('q75','q75')])),
            col=block.colors[1], border=NA)
}
# ... and add the lighter 5-95% range polygon before the median bars too...
for (bb in 1:nblocks) {
    times.beg.end <- c(block.years.center[bb]-halfwidth, block.years.center[bb]+halfwidth)
    polygon(c(times.beg.end,rev(times.beg.end)), c(returnlevel.quantiles[bb,c('q05','q05')],rev(returnlevel.quantiles[bb,c('q95','q95')])),
            col=block.colors.lighter[1], border=NA)
}
# ... so the bars are on top
for (bb in 1:nblocks) {lines(c(block.years.center[bb]-halfwidth, block.years.center[bb]+halfwidth),
                                   rep(returnlevel.quantiles[bb,'q50'],2), lwd=3, col='black')}
text(1895, 0.5, site, pos=4)
mtext('Year', side=1, line=2.4, cex=1);
mtext('100-year return level [m]', side=2, line=2.2, cex=1);
dev.off()

#===============================================================================
#



#
#===============================================================================
# FIGURE 2 - Top row: current surge levels; bottom row: projected 2065 surge
#            levels relative to present. Columns: different sites.
#            Projected distributions of 100-year surge level by BMA, relative to
#            each of the individual model structures.

# either do the convolution by fitting kernel density estimates, or sample all
# (or a large number of) possible combinations of BMA and other model ensemble
# (n.ensemble * n.ensemble possibilities) SOW to construct the distribution of
# [BMA] - [specific model]
# (old commits have the convolution codes)

# all index combinations between the BMA ensemble and each specific model ensemble
# don't even set up 'all.indices', because it is probably going to be HUGE (and
# in R, therefore slow)
#all.indices <- expand.grid(bma=1:n.ensemble, model=1:n.ensemble)
n.sample <- 1e6 # 1e6 chosen because quantiles seem to stabilize; 1e4 not enough,
                # but 1e5 and 1e6 yield similar results across all 3 sites
sample.indices <- expand.grid(bma=1:n.ensemble, model=1:n.ensemble)[sample(x=(1:(n.ensemble*n.ensemble)), size=n.sample, replace=FALSE),]

# * rl100.bma_mod is to sample from the distribution of [BMA - model] for 2016
#   100-year return levels and 2065 return levels relative to 2016.
# * f.rl100.bma_mod is to fit KDEs for plotting nicely
# * rl100.x are x coordinates for the KDEs, in meters, because we will start using meters after the density estimates
rl100.bma_mod        <- f.rl100.bma_mod        <- vector('list', nsites)
names(rl100.bma_mod) <- names(f.rl100.bma_mod) <- site.names
rl100.dx <- 0.01 # 1 cm sea level is the minimum difference we care about
rl100.x <- seq(-30, 30, by=rl100.dx)

for (site in site.names) {
  rl100.bma_mod[[site]]        <- f.rl100.bma_mod[[site]]        <- vector('list', nmodel)
  names(rl100.bma_mod[[site]]) <- names(f.rl100.bma_mod[[site]]) <- types.of.gpd
  for (model in types.of.gpd) {
    rl100.bma_mod[[site]][[model]]        <- f.rl100.bma_mod[[site]][[model]]        <- vector('list', length(year.names))
    names(rl100.bma_mod[[site]][[model]]) <- names(f.rl100.bma_mod[[site]][[model]]) <- year.names
    rl100.bma_mod[[site]][[model]]$y2016 <- rl100.bma[[site]][[all.data[[site]]]]$y2016[sample.indices[,1]] - rl100[[site]][[all.data[[site]]]][[model]]$y2016[sample.indices[,2]]
    rl100.bma_mod[[site]][[model]]$y2065 <- rl100.bma[[site]][[all.data[[site]]]]$y2065[sample.indices[,1]] - rl100[[site]][[all.data[[site]]]][[model]]$y2065[sample.indices[,2]]
    f.rl100.bma_mod[[site]][[model]]$y2016 <- density(x=rl100.bma_mod[[site]][[model]]$y2016, from=1000*min(rl100.x), to=1000*max(rl100.x), n=length(rl100.x), kernel='gaussian')
    f.rl100.bma_mod[[site]][[model]]$y2065 <- density(x=rl100.bma_mod[[site]][[model]]$y2065, from=1000*min(rl100.x), to=1000*max(rl100.x), n=length(rl100.x), kernel='gaussian')
    # normalize
    f.rl100.bma_mod[[site]][[model]]$y2016 <- f.rl100.bma_mod[[site]][[model]]$y2016$y/sintegral(x=rl100.x, fx=f.rl100.bma_mod[[site]][[model]]$y2016$y)$value
    f.rl100.bma_mod[[site]][[model]]$y2065 <- f.rl100.bma_mod[[site]][[model]]$y2065$y/sintegral(x=rl100.x, fx=f.rl100.bma_mod[[site]][[model]]$y2065$y)$value
  }
}

# box-whiskers for the plot

## Calculate the quantiles to plot
quantiles.to.grab <- c(.05, .25, .5, .75, .95)
quantile.names <- rep(NULL, length(quantiles.to.grab))
for (qq in 1:length(quantiles.to.grab)) {
  if(quantiles.to.grab[qq] >= .10) {
    quantile.names[qq] <- paste('q',100*quantiles.to.grab[qq], sep='')
  } else if(quantiles.to.grab[qq] < .10 & quantiles.to.grab[qq] >= 0) {
    quantile.names[qq] <- paste('q0',100*quantiles.to.grab[qq], sep='')
  }
}
# initialize like this, and then within each of the 'year.names' (last element)
# replace with the quantiles and rename 'quantile.names'
q.rl100.bma_mod <- rl100.bma_mod
for (site in site.names) {
  for (model in types.of.gpd) {
    for (year in year.names[2:3]) {
      q.rl100.bma_mod[[site]][[model]][[year]] <- quantile(rl100.bma_mod[[site]][[model]][[year]], quantiles.to.grab)
      names(q.rl100.bma_mod[[site]][[model]][[year]]) <- quantile.names
    }
  }
}

medians <- vector('list', nsites); names(medians) <- site.names
widths.50 <- vector('list', nsites); names(widths.50) <- site.names
widths.90 <- vector('list', nsites); names(widths.90) <- site.names
for (site in site.names) {
  medians[[site]] <- widths.50[[site]] <- widths.90[[site]] <- mat.or.vec(2, nmodel)
  rownames(medians[[site]]) <- rownames(widths.50[[site]]) <- rownames(widths.90[[site]]) <- year.names[2:3]
  colnames(medians[[site]]) <- colnames(widths.50[[site]]) <- colnames(widths.90[[site]]) <- types.of.gpd
  for (model in types.of.gpd) {
    for (year in year.names[2:3]) {
      # 0.001 * is to convert mm to m
      widths.50[[site]][year,model] <- 0.001 * (q.rl100.bma_mod[[site]][[model]][[year]]['q75'] - q.rl100.bma_mod[[site]][[model]][[year]]['q25'])
      widths.90[[site]][year,model] <- 0.001 * (q.rl100.bma_mod[[site]][[model]][[year]]['q95'] - q.rl100.bma_mod[[site]][[model]][[year]]['q05'])
      medians[[site]][year,model] <- 0.001 * q.rl100.bma_mod[[site]][[model]][[year]]['q50']
    }
  }
}

# Calculate the survival functions (empirical)

# empirical cumulative density (distribution) and survival function
# (the probabilities)
ecdf.vals.bma <- seq(from=0, to=1, length.out=n.sample)
esf.vals.bma  <- 1 - ecdf.vals.bma
ecdf.vals <- seq(from=0, to=1, length.out=n.ensemble)
esf.vals  <- 1 - ecdf.vals
x.lims <- c(0, 100, 10000) # lower lim, upper lim, number of nodes

# values corresponding to those probabilities
sf.rl100.bma_mod <- vector('list', nsites); names(sf.rl100.bma_mod) <- site.names
sf.rl100 <- vector('list', nsites); names(sf.rl100) <- site.names
pdf.rl100 <- vector('list', nsites); names(pdf.rl100) <- site.names

for (site in site.names) {
  sf.rl100.bma_mod[[site]] <- vector('list', nmodel); names(sf.rl100.bma_mod[[site]]) <- types.of.gpd
  sf.rl100[[site]] <- vector('list', nmodel+1); names(sf.rl100[[site]]) <- c(types.of.gpd, 'bma')
  pdf.rl100[[site]] <- vector('list', nmodel+1); names(pdf.rl100[[site]]) <- c(types.of.gpd, 'bma')
  for (model in types.of.gpd) {
    sf.rl100.bma_mod[[site]][[model]] <- vector('list', nyears); names(sf.rl100.bma_mod[[site]][[model]]) <- year.names
    sf.rl100[[site]][[model]] <- vector('list', nyears); names(sf.rl100[[site]][[model]]) <- year.names
    pdf.rl100[[site]][[model]] <- vector('list', nyears); names(pdf.rl100[[site]][[model]]) <- year.names
    for (year in year.names[2:3]) {
      sf.rl100.bma_mod[[site]][[model]][[year]] <- 0.001*rl100.bma_mod[[site]][[model]][[year]][order(rl100.bma_mod[[site]][[model]][[year]])]
      sf.rl100[[site]][[model]][[year]] <- 0.001*rl100[[site]][[all.data[[site]]]][[model]][[year]][order(rl100[[site]][[all.data[[site]]]][[model]][[year]])]
      pdf.rl100[[site]][[model]][[year]] <- density(x=0.001*rl100[[site]][[all.data[[site]]]][[model]][[year]], from=x.lims[1], to=x.lims[2], n=x.lims[3], kernel='gaussian')
    }
  }
  for (year in year.names[2:3]) {
    sf.rl100[[site]]$bma[[year]] <- 0.001*rl100.bma[[site]][[all.data[[site]]]][[year]][order(rl100.bma[[site]][[all.data[[site]]]][[year]])]
    pdf.rl100[[site]]$bma[[year]] <- density(x=0.001*rl100.bma[[site]][[all.data[[site]]]][[year]], from=x.lims[1], to=x.lims[2], n=x.lims[3], kernel='gaussian')
  }
}

#
# cap the pdfs that are just spikes and give no useful information
# otherwise they'll just be hard to plot with good, god-fearing pdfs
#
capping <- function(fx, softcap, hardcap) {
  ind_high <- which(fx > softcap)
  if(length(ind_high)==0) {
    return(fx)
  } else {
    scl <- ((fx[ind_high] - softcap)/(max(fx)-softcap)) * 0.5 * (hardcap-softcap)
    fx_cap <- fx
    fx_cap[ind_high] <- scl + softcap
    return(fx_cap)
  }
}

saveRDS(rl100.bma_mod, file='bma_minus_model_anomalies.rds')

#
# the actual figure -- WITH pdfs and survival function of ACTUAL SURGE LEVELS
#

pdf(paste(plot.dir,'returnlevels_pdf_sf.pdf',sep=''),width=7,height=9.5,colormodel='rgb')

par(mfrow=c(4,2), mai=c(.67,.60,.25,.01))
#####################
# pdfs in year 2016 #
#####################
year <- 'y2016'
site <- 'Delfzijl'
plot(pdf.rl100[[site]]$bma[[year]]$x, pdf.rl100[[site]]$bma[[year]]$y, type='l', xlim=c(3.75, 6.25), ylim=c(0, 3),
     lwd=2.5, lty=1, xlab='', ylab='', xaxs='i', yaxs='i', yaxt='n', axes=FALSE, col='black')
lines(pdf.rl100[[site]]$gpd4[[year]]$x, pdf.rl100[[site]]$gpd4[[year]]$y, col='darkorange3', lwd=2, lty=2)
lines(pdf.rl100[[site]]$gpd5[[year]]$x, pdf.rl100[[site]]$gpd5[[year]]$y, col='mediumslateblue', lwd=2, lty=2)
lines(pdf.rl100[[site]]$gpd6[[year]]$x, pdf.rl100[[site]]$gpd6[[year]]$y, col='mediumvioletred', lwd=2, lty=2)
lines(pdf.rl100[[site]]$gpd3[[year]]$x, pdf.rl100[[site]]$gpd3[[year]]$y, col='seagreen', lwd=2, lty=4)
u <- par("usr")
arrows(u[1], u[3], u[1], .99*u[4], code = 2, length=.15, xpd = TRUE)
mtext('Probability density', side=2, line=1, cex=.75);
mtext('Delfzijl', side=3, line=.6, cex=0.75)
mtext(side=3, text=expression(bold(' a.')), line=.6, cex=0.75, adj=0)
mtext('100-year return level in 2016 [m]', side=1, line=3.4, cex=0.75);
axis(1, at=seq(3.75, 6.25, .25), labels=c('','4','','4.5','','5','','5.5','','6',''), cex.axis=1.1)
legend(5.2, 3, c('ST','NS1','NS2','NS3','BMA'), lty=c(4,2,2,2,1), cex=1.1, bty='n', lwd=2,
       col=c('seagreen','darkorange3','mediumslateblue','mediumvioletred', 'black'))
#
par(mai=c(.67,.42,.25,.19))
site <- 'Norfolk'
plot(pdf.rl100[[site]]$bma[[year]]$x, pdf.rl100[[site]]$bma[[year]]$y, type='l', xlim=c(1.5, 4), ylim=c(0, 3),
     lwd=2.5, lty=1, xlab='', ylab='', xaxs='i', yaxs='i', yaxt='n', axes=FALSE, col='black')
lines(pdf.rl100[[site]]$gpd4[[year]]$x, pdf.rl100[[site]]$gpd4[[year]]$y, col='darkorange3', lwd=2, lty=2)
lines(pdf.rl100[[site]]$gpd5[[year]]$x, pdf.rl100[[site]]$gpd5[[year]]$y, col='mediumslateblue', lwd=2, lty=2)
lines(pdf.rl100[[site]]$gpd6[[year]]$x, pdf.rl100[[site]]$gpd6[[year]]$y, col='mediumvioletred', lwd=2, lty=2)
lines(pdf.rl100[[site]]$gpd3[[year]]$x, pdf.rl100[[site]]$gpd3[[year]]$y, col='seagreen', lwd=2, lty=4)
#lines(c(0,0), c(-100,100), col='black', lwd=1.5)
lines(c(-100,100), c(0,0), col='black', lwd=1.5) # replace the x-axis wonkily
u <- par("usr")
arrows(u[1], u[3], u[1], .99*u[4], code = 2, length=.15, xpd = TRUE)
mtext('Probability density', side=2, line=1, cex=.75);
mtext('Norfolk', side=3, line=.6, cex=0.75)
mtext(side=3, text=expression(bold(' b.')), line=.6, cex=0.75, adj=0)
mtext('100-year return level in 2016 [m]', side=1, line=3.4, cex=0.75);
axis(1, at=seq(1.5, 4, .25), labels=c('1.5','','2','','2.5','','3','','3.5','','4'), cex.axis=1.1)

#####################
#   sf in year 2016 #
#####################

par(mai=c(.67,.60,.25,.01))
site <- 'Delfzijl'
plot(c(0,sf.rl100[[site]]$bma[[year]]), c(0,log10(esf.vals)), type='l', xlim=c(3.75, 6.25), ylim=c(-2.5, 0.2),
     lwd=2.5, lty=1, xlab='', ylab='', xaxs='i', yaxs='i', yaxt='n', axes=FALSE, col='black')
lines(sf.rl100[[site]]$gpd4[[year]], log10(esf.vals), col='darkorange3', lwd=2, lty=2)
lines(sf.rl100[[site]]$gpd5[[year]], log10(esf.vals), col='mediumslateblue', lwd=2, lty=2)
lines(sf.rl100[[site]]$gpd6[[year]], log10(esf.vals), col='mediumvioletred', lwd=2, lty=2)
lines(sf.rl100[[site]]$gpd3[[year]], log10(esf.vals), col='seagreen', lwd=2, lty=4)
mtext('Survival function', side=2, line=3, cex=.75);
mtext(side=3, text=expression(bold(' c.')), line=.1, cex=0.75, adj=0)
mtext('100-year return level in 2016 [m]', side=1, line=3.4, cex=0.75);
axis(1, at=seq(3.75, 6.25, .25), labels=c('','4','','4.5','','5','','5.5','','6',''), cex.axis=1.1)
axis(2, at=seq(-4,0), label=parse(text=paste("10^", seq(-4,0), sep="")), las=1, cex.axis=1.1, mgp=c(3,.7,0))
#
par(mai=c(.67,.42,.25,.19))
site <- 'Norfolk'
plot(c(0,sf.rl100[[site]]$bma[[year]]), c(0,log10(esf.vals)), type='l', xlim=c(1.5, 4), ylim=c(-2.5, 0.2),
     lwd=2.5, lty=1, xlab='', ylab='', xaxs='i', yaxs='i', yaxt='n', axes=FALSE, col='black')
lines(sf.rl100[[site]]$gpd4[[year]], log10(esf.vals), col='darkorange3', lwd=2, lty=2)
lines(sf.rl100[[site]]$gpd5[[year]], log10(esf.vals), col='mediumslateblue', lwd=2, lty=2)
lines(sf.rl100[[site]]$gpd6[[year]], log10(esf.vals), col='mediumvioletred', lwd=2, lty=2)
lines(sf.rl100[[site]]$gpd3[[year]], log10(esf.vals), col='seagreen', lwd=2, lty=4)
mtext(side=3, text=expression(bold(' d.')), line=.1, cex=0.75, adj=0)
mtext('100-year return level in 2016 [m]', side=1, line=3.4, cex=0.75);
axis(1, at=seq(1.5, 4, .25), labels=c('1.5','','2','','2.5','','3','','3.5','','4'), cex.axis=1.1)
axis(2, at=seq(-4,0), label=parse(text=paste("10^", seq(-4,0), sep="")), las=1, cex.axis=1.1, mgp=c(3,.7,0))

#####################
# pdfs in year 2065 #
#####################
year <- 'y2065'
#
par(mai=c(.67,.60,.25,.01))
site <- 'Delfzijl'
plot(pdf.rl100[[site]]$bma[[year]]$x, pdf.rl100[[site]]$bma[[year]]$y, type='l', xlim=c(3.75, 6.25), ylim=c(0, 3),
     lwd=2.5, lty=1, xlab='', ylab='', xaxs='i', yaxs='i', yaxt='n', axes=FALSE, col='black')
lines(pdf.rl100[[site]]$gpd4[[year]]$x, pdf.rl100[[site]]$gpd4[[year]]$y, col='darkorange3', lwd=2, lty=2)
lines(pdf.rl100[[site]]$gpd5[[year]]$x, pdf.rl100[[site]]$gpd5[[year]]$y, col='mediumslateblue', lwd=2, lty=2)
lines(pdf.rl100[[site]]$gpd6[[year]]$x, pdf.rl100[[site]]$gpd6[[year]]$y, col='mediumvioletred', lwd=2, lty=2)
lines(pdf.rl100[[site]]$gpd3[[year]]$x, pdf.rl100[[site]]$gpd3[[year]]$y, col='seagreen', lwd=2, lty=4)
u <- par("usr")
arrows(u[1], u[3], u[1], .99*u[4], code = 2, length=.15, xpd = TRUE)
mtext('Probability density', side=2, line=1, cex=.75);
mtext(side=3, text=expression(bold(' e.')), line=.6, cex=0.75, adj=0)
mtext('100-year return level in 2065 [m]', side=1, line=3.4, cex=0.75);
axis(1, at=seq(3.75, 6.25, .25), labels=c('','4','','4.5','','5','','5.5','','6',''), cex.axis=1.1)
#
par(mai=c(.67,.42,.25,.19))
site <- 'Norfolk'
plot(pdf.rl100[[site]]$bma[[year]]$x, pdf.rl100[[site]]$bma[[year]]$y, type='l', xlim=c(1.5, 4), ylim=c(0, 3),
     lwd=2.5, lty=1, xlab='', ylab='', xaxs='i', yaxs='i', yaxt='n', axes=FALSE, col='black')
lines(pdf.rl100[[site]]$gpd4[[year]]$x, pdf.rl100[[site]]$gpd4[[year]]$y, col='darkorange3', lwd=2, lty=2)
lines(pdf.rl100[[site]]$gpd5[[year]]$x, pdf.rl100[[site]]$gpd5[[year]]$y, col='mediumslateblue', lwd=2, lty=2)
lines(pdf.rl100[[site]]$gpd6[[year]]$x, pdf.rl100[[site]]$gpd6[[year]]$y, col='mediumvioletred', lwd=2, lty=2)
lines(pdf.rl100[[site]]$gpd3[[year]]$x, pdf.rl100[[site]]$gpd3[[year]]$y, col='seagreen', lwd=2, lty=4)
u <- par("usr")
arrows(u[1], u[3], u[1], .99*u[4], code = 2, length=.15, xpd = TRUE)
mtext(side=3, text=expression(bold(' f.')), line=.6, cex=0.75, adj=0)
mtext('100-year return level in 2065 [m]', side=1, line=3.4, cex=0.75);
axis(1, at=seq(1.5, 4, .25), labels=c('1.5','','2','','2.5','','3','','3.5','','4'), cex.axis=1.1)

#####################
#   sf in year 2065 #
#####################

par(mai=c(.67,.60,.25,.01))
site <- 'Delfzijl'
plot(c(0,sf.rl100[[site]]$bma[[year]]), c(0,log10(esf.vals)), type='l', xlim=c(3.75, 6.25), ylim=c(-2.5, 0.2),
     lwd=2.5, lty=1, xlab='', ylab='', xaxs='i', yaxs='i', yaxt='n', axes=FALSE, col='black')
lines(sf.rl100[[site]]$gpd4[[year]], log10(esf.vals), col='darkorange3', lwd=2, lty=2)
lines(sf.rl100[[site]]$gpd5[[year]], log10(esf.vals), col='mediumslateblue', lwd=2, lty=2)
lines(sf.rl100[[site]]$gpd6[[year]], log10(esf.vals), col='mediumvioletred', lwd=2, lty=2)
lines(sf.rl100[[site]]$gpd3[[year]], log10(esf.vals), col='seagreen', lwd=2, lty=4)
mtext('Survival function', side=2, line=3, cex=.75);
mtext(side=3, text=expression(bold(' g.')), line=.1, cex=0.75, adj=0)
mtext('100-year return level in 2065 [m]', side=1, line=3.4, cex=0.75);
axis(1, at=seq(3.75, 6.25, .25), labels=c('','4','','4.5','','5','','5.5','','6',''), cex.axis=1.1)
axis(2, at=seq(-4,0), label=parse(text=paste("10^", seq(-4,0), sep="")), las=1, cex.axis=1.1, mgp=c(3,.7,0))
#
par(mai=c(.67,.42,.25,.19))
site <- 'Norfolk'
plot(c(0,sf.rl100[[site]]$bma[[year]]), c(0,log10(esf.vals)), type='l', xlim=c(1.5, 4), ylim=c(-2.5, 0.2),
     lwd=2.5, lty=1, xlab='', ylab='', xaxs='i', yaxs='i', yaxt='n', axes=FALSE, col='black')
lines(sf.rl100[[site]]$gpd4[[year]], log10(esf.vals), col='darkorange3', lwd=2, lty=2)
lines(sf.rl100[[site]]$gpd5[[year]], log10(esf.vals), col='mediumslateblue', lwd=2, lty=2)
lines(sf.rl100[[site]]$gpd6[[year]], log10(esf.vals), col='mediumvioletred', lwd=2, lty=2)
lines(sf.rl100[[site]]$gpd3[[year]], log10(esf.vals), col='seagreen', lwd=2, lty=4)
mtext(side=3, text=expression(bold(' h.')), line=.1, cex=0.75, adj=0)
mtext('100-year return level in 2065 [m]', side=1, line=3.4, cex=0.75);
axis(1, at=seq(1.5, 4, .25), labels=c('1.5','','2','','2.5','','3','','3.5','','4'), cex.axis=1.1)
axis(2, at=seq(-4,0), label=parse(text=paste("10^", seq(-4,0), sep="")), las=1, cex.axis=1.1, mgp=c(3,.7,0))

dev.off()

#===============================================================================
#


#===============================================================================

#
# the actual FIGURE 2 -- pdfs and bar-plots of ACTUAL SURGE LEVELS
#

pdf(paste(plot.dir,'returnlevels_pdf_bar.pdf',sep=''),width=7,height=5.5,colormodel='rgb')

par(mfrow=c(2,2), mai=c(.55,.50,.35,.02))
#####################
# pdfs in year 2016 #
#####################
year <- 'y2016'
site <- 'Delfzijl'
offset = 0.5
plot(pdf.rl100[[site]]$bma[[year]]$x, offset+pdf.rl100[[site]]$bma[[year]]$y, type='l', xlim=c(3.75, 5.75), ylim=c(0, 3+offset),
     lwd=2.5, lty=1, xlab='', ylab='', xaxs='i', yaxs='i', yaxt='n', axes=FALSE, col='black')
lines(pdf.rl100[[site]]$gpd4[[year]]$x, offset+pdf.rl100[[site]]$gpd4[[year]]$y, col='darkorange3', lwd=2, lty=2)
lines(pdf.rl100[[site]]$gpd5[[year]]$x, offset+pdf.rl100[[site]]$gpd5[[year]]$y, col='mediumslateblue', lwd=2, lty=2)
lines(pdf.rl100[[site]]$gpd6[[year]]$x, offset+pdf.rl100[[site]]$gpd6[[year]]$y, col='mediumvioletred', lwd=2, lty=2)
lines(pdf.rl100[[site]]$gpd3[[year]]$x, offset+pdf.rl100[[site]]$gpd3[[year]]$y, col='seagreen', lwd=2, lty=4)
y1 = c(.29, .45)
polygon(c(q.rl100.mod[[site]]['gpd3',c(2,2)], q.rl100.mod[[site]]['gpd3',c(6,6)]), c(y1,rev(y1)), col='palegreen1', border=NA)
polygon(c(q.rl100.mod[[site]]['gpd3',c(3,3)], q.rl100.mod[[site]]['gpd3',c(5,5)]), c(y1,rev(y1)), col='seagreen3', border=NA)
lines(c(q.rl100.mod[[site]]['gpd3',4], q.rl100.mod[[site]]['gpd3',4]), y1, lwd=2.5, col='olivedrab')
y2 = c(.06, .23)
polygon(c(q.rl100.mod[[site]]['bma',c(2,2)], q.rl100.mod[[site]]['bma',c(6,6)]), c(y2,rev(y2)), col='gray60', border=NA)
polygon(c(q.rl100.mod[[site]]['bma',c(3,3)], q.rl100.mod[[site]]['bma',c(5,5)]), c(y2,rev(y2)), col='gray30', border=NA)
lines(c(q.rl100.mod[[site]]['bma',4], q.rl100.mod[[site]]['bma',4]), y2, lwd=2.5, col='black')
u <- par("usr")
arrows(u[1], u[3], u[1], .99*u[4], code = 2, length=.15, xpd = TRUE)
mtext('Probability density', side=2, line=1, cex=1);
mtext('Delfzijl', side=3, line=.6, cex=1)
mtext(side=3, text=expression(bold(' a.')), line=.6, cex=1, adj=0)
mtext('100-year return level in 2016 [m]', side=1, line=2.4, cex=0.9);
axis(1, at=seq(3.75, 6.25, .25), labels=c('','4','','4.5','','5','','5.5','','6',''), cex.axis=1.05)
legend(5, 3.6, c('ST','NS1','NS2','NS3','BMA'), lty=c(4,2,2,2,1), cex=1, bty='n', lwd=2,
       col=c('seagreen','darkorange3','mediumslateblue','mediumvioletred', 'black'))
#
par(mai=c(.55,.26,.35,.26))
site <- 'Norfolk'
offset = 0.5
plot(pdf.rl100[[site]]$bma[[year]]$x, offset+pdf.rl100[[site]]$bma[[year]]$y, type='l', xlim=c(1.5, 3.5), ylim=c(0, 3.05+offset),
     lwd=2.5, lty=1, xlab='', ylab='', xaxs='i', yaxs='i', yaxt='n', axes=FALSE, col='black')
lines(pdf.rl100[[site]]$gpd4[[year]]$x, offset+pdf.rl100[[site]]$gpd4[[year]]$y, col='darkorange3', lwd=2, lty=2)
lines(pdf.rl100[[site]]$gpd5[[year]]$x, offset+pdf.rl100[[site]]$gpd5[[year]]$y, col='mediumslateblue', lwd=2, lty=2)
lines(pdf.rl100[[site]]$gpd6[[year]]$x, offset+pdf.rl100[[site]]$gpd6[[year]]$y, col='mediumvioletred', lwd=2, lty=2)
lines(pdf.rl100[[site]]$gpd3[[year]]$x, offset+pdf.rl100[[site]]$gpd3[[year]]$y, col='seagreen', lwd=2, lty=4)
y1 = c(.29, .45)
polygon(c(q.rl100.mod[[site]]['gpd3',c(9,9)], q.rl100.mod[[site]]['gpd3',c('q95','q95')]), c(y1,rev(y1)), col='palegreen1', border=NA)
polygon(c(q.rl100.mod[[site]]['gpd3',c(10,10)], q.rl100.mod[[site]]['gpd3',c('q75','q75')]), c(y1,rev(y1)), col='seagreen3', border=NA)
lines(c(q.rl100.mod[[site]]['gpd3',11], q.rl100.mod[[site]]['gpd3','q50']), y1, lwd=2.5, col='olivedrab')
y2 = c(.06, .23)
polygon(c(q.rl100.mod[[site]]['bma',c(9,9)], q.rl100.mod[[site]]['bma',c(13,13)]), c(y2,rev(y2)), col='gray60', border=NA)
polygon(c(q.rl100.mod[[site]]['bma',c(10,10)], q.rl100.mod[[site]]['bma',c(12,12)]), c(y2,rev(y2)), col='gray30', border=NA)
lines(c(q.rl100.mod[[site]]['bma',11], q.rl100.mod[[site]]['bma',11]), y2, lwd=2.5, col='black')
#lines(c(0,0), c(-100,100), col='black', lwd=1.5)
lines(c(-100,100), c(0,0), col='black', lwd=1.5) # replace the x-axis wonkily
u <- par("usr")
arrows(u[1], u[3], u[1], .99*u[4], code = 2, length=.15, xpd = TRUE)
mtext('Norfolk', side=3, line=.6, cex=1)
mtext(side=3, text=expression(bold(' b.')), line=.6, cex=1, adj=0)
mtext('100-year return level in 2016 [m]', side=1, line=2.4, cex=.9);
axis(1, at=seq(1.5, 4, .25), labels=c('1.5','','2','','2.5','','3','','3.5','','4'), cex.axis=1.05)

#####################
# pdfs in year 2065 #
#####################
year <- 'y2065'
#
par(mai=c(.65,.50,.3,.02))
site <- 'Delfzijl'
offset = 0.5
plot(pdf.rl100[[site]]$bma[[year]]$x, offset+pdf.rl100[[site]]$bma[[year]]$y, type='l', xlim=c(3.75, 5.75), ylim=c(0, 3+offset),
     lwd=2.5, lty=1, xlab='', ylab='', xaxs='i', yaxs='i', yaxt='n', axes=FALSE, col='black')
lines(pdf.rl100[[site]]$gpd4[[year]]$x, offset+pdf.rl100[[site]]$gpd4[[year]]$y, col='darkorange3', lwd=2, lty=2)
lines(pdf.rl100[[site]]$gpd5[[year]]$x, offset+pdf.rl100[[site]]$gpd5[[year]]$y, col='mediumslateblue', lwd=2, lty=2)
lines(pdf.rl100[[site]]$gpd6[[year]]$x, offset+pdf.rl100[[site]]$gpd6[[year]]$y, col='mediumvioletred', lwd=2, lty=2)
lines(pdf.rl100[[site]]$gpd3[[year]]$x, offset+pdf.rl100[[site]]$gpd3[[year]]$y, col='seagreen', lwd=2, lty=4)
y1 = c(.29, .45)
polygon(c(q.rl100.mod[[site]]['gpd3',c(9,9)], q.rl100.mod[[site]]['gpd3',c('q95','q95')]), c(y1,rev(y1)), col='palegreen1', border=NA)
polygon(c(q.rl100.mod[[site]]['gpd3',c(10,10)], q.rl100.mod[[site]]['gpd3',c('q75','q75')]), c(y1,rev(y1)), col='seagreen3', border=NA)
lines(c(q.rl100.mod[[site]]['gpd3',11], q.rl100.mod[[site]]['gpd3','q50']), y1, lwd=2.5, col='olivedrab')
y2 = c(.06, .23)
polygon(c(q.rl100.mod[[site]]['bma',c(9,9)], q.rl100.mod[[site]]['bma',c(13,13)]), c(y2,rev(y2)), col='gray60', border=NA)
polygon(c(q.rl100.mod[[site]]['bma',c(10,10)], q.rl100.mod[[site]]['bma',c(12,12)]), c(y2,rev(y2)), col='gray30', border=NA)
lines(c(q.rl100.mod[[site]]['bma',11], q.rl100.mod[[site]]['bma',11]), y2, lwd=2.5, col='black')
u <- par("usr")
arrows(u[1], u[3], u[1], .99*u[4], code = 2, length=.15, xpd = TRUE)
mtext('Probability density', side=2, line=1, cex=.9);
mtext(side=3, text=expression(bold(' c.')), line=.6, cex=1, adj=0)
mtext('100-year return level in 2065 [m]', side=1, line=2.4, cex=0.9);
axis(1, at=seq(3.75, 6.25, .25), labels=c('','4','','4.5','','5','','5.5','','6',''), cex.axis=1.05)
#
par(mai=c(.65,.26,.3,.26))
site <- 'Norfolk'
offset = 0.5
plot(pdf.rl100[[site]]$bma[[year]]$x, offset+pdf.rl100[[site]]$bma[[year]]$y, type='l', xlim=c(1.5, 3.5), ylim=c(0, 3.05+offset),
     lwd=2.5, lty=1, xlab='', ylab='', xaxs='i', yaxs='i', yaxt='n', axes=FALSE, col='black')
lines(pdf.rl100[[site]]$gpd4[[year]]$x, offset+pdf.rl100[[site]]$gpd4[[year]]$y, col='darkorange3', lwd=2, lty=2)
lines(pdf.rl100[[site]]$gpd5[[year]]$x, offset+pdf.rl100[[site]]$gpd5[[year]]$y, col='mediumslateblue', lwd=2, lty=2)
lines(pdf.rl100[[site]]$gpd6[[year]]$x, offset+pdf.rl100[[site]]$gpd6[[year]]$y, col='mediumvioletred', lwd=2, lty=2)
lines(pdf.rl100[[site]]$gpd3[[year]]$x, offset+pdf.rl100[[site]]$gpd3[[year]]$y, col='seagreen', lwd=2, lty=4)
y1 = c(.29, .45)
polygon(c(q.rl100.mod[[site]]['gpd3',c(9,9)], q.rl100.mod[[site]]['gpd3',c('q95','q95')]), c(y1,rev(y1)), col='palegreen1', border=NA)
polygon(c(q.rl100.mod[[site]]['gpd3',c(10,10)], q.rl100.mod[[site]]['gpd3',c('q75','q75')]), c(y1,rev(y1)), col='seagreen3', border=NA)
lines(c(q.rl100.mod[[site]]['gpd3',11], q.rl100.mod[[site]]['gpd3','q50']), y1, lwd=2.5, col='olivedrab')
y2 = c(.06, .23)
polygon(c(q.rl100.mod[[site]]['bma',c(9,9)], q.rl100.mod[[site]]['bma',c(13,13)]), c(y2,rev(y2)), col='gray60', border=NA)
polygon(c(q.rl100.mod[[site]]['bma',c(10,10)], q.rl100.mod[[site]]['bma',c(12,12)]), c(y2,rev(y2)), col='gray30', border=NA)
lines(c(q.rl100.mod[[site]]['bma',11], q.rl100.mod[[site]]['bma',11]), y2, lwd=2.5, col='black')
u <- par("usr")
arrows(u[1], u[3], u[1], .99*u[4], code = 2, length=.15, xpd = TRUE)
mtext(side=3, text=expression(bold(' d.')), line=.6, cex=1, adj=0)
mtext('100-year return level in 2065 [m]', side=1, line=2.4, cex=0.9);
axis(1, at=seq(1.5, 4, .25), labels=c('1.5','','2','','2.5','','3','','3.5','','4'), cex.axis=1.05)

dev.off()

#===============================================================================
#




#
#===============================================================================
# Further analysis:
#==================

# By how much do the distribution medians differ between BMA and the ST model?

print('BMA - ST medians in 2016:')
print(paste('  Delfzijl: ',median(rl100.bma$Delfzijl$y130$y2016)-median(rl100$Delfzijl$y130$gpd3$y2016), sep=''))
print(paste('  Norfolk: ',median(rl100.bma$Norfolk$y90$y2016)-median(rl100$Norfolk$y90$gpd3$y2016), sep=''))

print('BMA - ST medians in 2065:')
print(paste('  Delfzijl:',median(rl100.bma$Delfzijl$y130$y2065)-median(rl100$Delfzijl$y130$gpd3$y2065), sep=''))
print(paste('  Norfolk:',median(rl100.bma$Norfolk$y90$y2065)-median(rl100$Norfolk$y90$gpd3$y2065), sep=''))

#===============================================================================
#




#
#===============================================================================
# FIGURE 3 - Box-whisker (or box-lighter-box) distributions (horizontally) of
#            100-year return level for varying lengths of data employed for the
#            BMA ensemble. Different sites are 3 horizontally oriented panels.

## Calculate the quantiles to plot
quantiles.to.grab <- c(.05, .25, .5, .75, .95)
quantile.names <- rep(NULL, length(quantiles.to.grab))
for (qq in 1:length(quantiles.to.grab)) {
  if(quantiles.to.grab[qq] >= .10) {
    quantile.names[qq] <- paste('q',100*quantiles.to.grab[qq], sep='')
  } else if(quantiles.to.grab[qq] < .10 & quantiles.to.grab[qq] >= 0) {
    quantile.names[qq] <- paste('q0',100*quantiles.to.grab[qq], sep='')
  }
}
year.for.quantiles <- year.names[2]
returnlevel.quantiles <- vector('list', nsites); names(returnlevel.quantiles) <- site.names
for (site in site.names) {
  returnlevel.quantiles[[site]] <- mat.or.vec(length(data.lengths[[site]]), length(quantiles.to.grab))
  rownames(returnlevel.quantiles[[site]]) <- data.lengths[[site]]
  colnames(returnlevel.quantiles[[site]]) <- quantile.names
  for (dd in 1:length(data.lengths[[site]])) {
    for (year in year.names) {
      # the /1000 is to convert to m from mm
      #returnlevel.quantiles[[site]][dd,] <- quantile(rl100[[site]][[dd]]$bma[[year.for.quantiles]], quantiles.to.grab)/1000
      returnlevel.quantiles[[site]][dd,] <- quantile(rl100.bma[[site]][[dd]][[year.for.quantiles]], quantiles.to.grab)/1000
    }
  }
}

## Useful for plotting - centers of the time blocks used in the experiments
y.datalengths <- seq(from=30, to=137, by=20)

#
# make the figure
#
block.colors <- rev(colorRampPalette(c("darkslateblue","royalblue","turquoise1"),space="Lab")(max(all.data)))
block.colors.lighter <- paste(block.colors, "70", sep="")

pdf(paste(plot.dir,'datalengths_boxwhisker.pdf',sep=''),width=5.5,height=4,colormodel='cmyk')
par(mfrow=c(1,2), mai=c(.8, .7, .4, .1))
halfwidth <- 2 # half the width of the boxes, in years
# put the first median bar down, to get hte plot started
site <- 'Delfzijl'
plot(rep(returnlevel.quantiles[[site]][1,'q50'],2), c(y.datalengths[1]-halfwidth, y.datalengths[1]+halfwidth),
     type='l', lwd=1.5, col='black', xlim=c(4, 7), ylim=c(130,30), xlab='', ylab='', las=1, xaxt='n', yaxt='n', cex.axis=1.15)
# ... and add the 25-75% range polygon...
for (dd in 1:length(data.lengths[[site]])) {
    times.beg.end <- c(y.datalengths[dd]-halfwidth, y.datalengths[dd]+halfwidth)
    polygon(c(returnlevel.quantiles[[site]][dd,c('q25','q25')],rev(returnlevel.quantiles[[site]][dd,c('q75','q75')])), c(times.beg.end,rev(times.beg.end)),
            col=block.colors[1], border=NA)
}
# ... and add the lighter 5-95% range polygon before the median bars too...
for (dd in 1:length(data.lengths[[site]])) {
    times.beg.end <- c(y.datalengths[dd]-halfwidth, y.datalengths[dd]+halfwidth)
    polygon(c(returnlevel.quantiles[[site]][dd,c('q05','q05')],rev(returnlevel.quantiles[[site]][dd,c('q95','q95')])), c(times.beg.end,rev(times.beg.end)),
            col=block.colors.lighter[1], border=NA)
}
# ... so the bars are on top
for (dd in 1:length(data.lengths[[site]])) {lines(rep(returnlevel.quantiles[[site]][dd,'q50'],2), c(y.datalengths[dd]-halfwidth, y.datalengths[dd]+halfwidth), lwd=1.5, col='black')}
mtext('Delfzijl', side=3, line=.6, cex=1)
mtext(side=3, text=expression(bold(' a.')), line=.6, cex=1, adj=0)
mtext('Years of data', side=2, line=2.4, cex=1);
mtext('100-year return level [m]', side=1, line=2.5, cex=1);
axis(1, at=seq(4, 7, by=0.5), labels=c('4','','5','','6','','7'), cex.axis=1)
axis(2, at=y.datalengths, labels=y.datalengths, cex.axis=1)

# put the first median bar down, to get hte plot started
site <- 'Norfolk'
plot(rep(returnlevel.quantiles[[site]][1,'q50'],2), c(y.datalengths[1]-halfwidth, y.datalengths[1]+halfwidth),
     type='l', lwd=1.5, col='black', xlim=c(1.5, 3.25), ylim=c(130,30), xlab='', ylab='', las=1, xaxt='n', yaxt='n', cex.axis=1.15)
# ... and add the 25-75% range polygon...
for (dd in 1:length(data.lengths[[site]])) {
    times.beg.end <- c(y.datalengths[dd]-halfwidth, y.datalengths[dd]+halfwidth)
    polygon(c(returnlevel.quantiles[[site]][dd,c('q25','q25')],rev(returnlevel.quantiles[[site]][dd,c('q75','q75')])), c(times.beg.end,rev(times.beg.end)),
            col=block.colors[1], border=NA)
}
# ... and add the lighter 5-95% range polygon before the median bars too...
for (dd in 1:length(data.lengths[[site]])) {
    times.beg.end <- c(y.datalengths[dd]-halfwidth, y.datalengths[dd]+halfwidth)
    polygon(c(returnlevel.quantiles[[site]][dd,c('q05','q05')],rev(returnlevel.quantiles[[site]][dd,c('q95','q95')])), c(times.beg.end,rev(times.beg.end)),
            col=block.colors.lighter[1], border=NA)
}
# ... so the bars are on top
for (dd in 1:length(data.lengths[[site]])) {lines(rep(returnlevel.quantiles[[site]][dd,'q50'],2), c(y.datalengths[dd]-halfwidth, y.datalengths[dd]+halfwidth), lwd=1.5, col='black')}
mtext('Norfolk', side=3, line=.6, cex=1)
mtext(side=3, text=expression(bold(' b.')), line=.6, cex=1, adj=0)
mtext('100-year return level [m]', side=1, line=2.5, cex=1);
axis(1, at=seq(1,3.5,by=0.25), labels=c('1','','1.5','','2','','2.5','','3','','3.5'), cex.axis=1)
axis(2, at=y.datalengths, labels=y.datalengths, cex.axis=1)

dev.off()

#===============================================================================
#


#
#===============================================================================
# FIGURE 4 - Bayesian model averaging weights (equation (XX)) for the four
#            candidate models, using (a) 30 years of tide gauge data from
#            Delfzijl, (b) 50 years of data, (c) 70 years of data, (d) 90 years
#            of data, (e) 110 years of data and (f) 137 years of data. Higher
#            values imply better model-data match. Different sites are different
#            columns.

library(plyr)
library(reshape2)
library(ggplot2)

bma.weights <- readRDS('../output/bma_weights_threshold99.rds')
bma.dfs <- lapply(bma.weights, data.frame, stringsAsFactors = FALSE)
bma.all <- ldply(bma.dfs, .id = 'Site')
bma.all$Model <- rep(c('ST','NS1','NS2','NS3'), 2)
bma.all$Model <- ordered(bma.all$Model, levels=c('ST', 'NS1', 'NS2', 'NS3'))
colnames(bma.all)[2:9] <- gsub("X", "", colnames(bma.all)[2:9], fixed=TRUE)
bma.all$ModelType <- ifelse(bma.all$Model == 'ST', 'Stationary', 'Nonstationary')
bma.all$ModelType <- ordered(bma.all$ModelType, levels=c('Stationary', 'Nonstationary'))

bma.melt <- melt(bma.all, id.vars=c('Site','Model', 'ModelType'),variable.name="DataLength", value.name= "ModelWeight")
bma.melt$DataLength <- as.numeric(levels(bma.melt$DataLength)[bma.melt$DataLength])

pdf(paste(plot.dir,'bma_weights_stationarytype.pdf',sep=''), height=3, width=5)
p <- ggplot(bma.melt[!is.na(bma.melt$ModelWeight),]) +
    geom_line(aes(x=DataLength, y=ModelWeight, linetype=ModelType, color=Model)) +
    geom_point(aes(x=DataLength, y=ModelWeight, color=Model), shape=20) +
    facet_grid(Site~.) + theme_bw(base_size=9) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_color_brewer('Model',type='qual',palette='Set2') +
    scale_linetype_discrete('Model Class') +
    scale_y_continuous('BMA Weight [dimensionless]', breaks=seq(0,0.6,0.1)) +
    scale_x_continuous('Data Length [years]', breaks=seq(30,130,20), expand=c(.01,.01)) +
    theme(strip.text.y = element_text(angle = 90))
print(p)
dev.off()
#===============================================================================
#



#
#===============================================================================
# SOM FIGURE S1 -- show the histograms of the MLE parameters and superimpose the
#                  prior distribution (normal or gamma)

deoptim.all <- readRDS('../output/surge_MLEs_ppgpd_decl3-pot99_28Mar2018.rds')
priors_normalgamma <- readRDS('../output/surge_priors_normalgamma_ppgpd_decl3-pot99_28Mar2018.rds')

nbins <- 12 # note that there are only 30 sites...
frac.ran <- 0.35 # extend axes above/below max/min range

# get standard limits for each of the 6 parameters
lims <- mat.or.vec(6,2)
pp <- 1 # lambda0
  parameters.pooled <- c(deoptim.all$gpd3[,1], deoptim.all$gpd4[,1], deoptim.all$gpd5[,1], deoptim.all$gpd6[,1])
  ran <- diff(quantile(parameters.pooled, c(0,1)))
  lims[pp,] <- c(min(parameters.pooled) - frac.ran*ran , max(parameters.pooled) + frac.ran*ran)
pp <- 2 # lambda1
  parameters.pooled <- c(deoptim.all$gpd4[,2], deoptim.all$gpd5[,2], deoptim.all$gpd6[,2])
  ran <- diff(quantile(parameters.pooled, c(0,1)))
  lims[pp,] <- c(min(parameters.pooled) - frac.ran*ran , max(parameters.pooled) + frac.ran*ran)
pp <- 3 # sigma0
  parameters.pooled <- c(log(deoptim.all$gpd3[,2]), log(deoptim.all$gpd4[,3]), deoptim.all$gpd5[,3], deoptim.all$gpd6[,3])
  ran <- diff(quantile(parameters.pooled, c(0,1)))
  lims[pp,] <- c(min(parameters.pooled) - frac.ran*ran , max(parameters.pooled) + frac.ran*ran)
pp <- 4 # sigma1
  parameters.pooled <- c(deoptim.all$gpd5[,4], deoptim.all$gpd6[,4])
  ran <- diff(quantile(parameters.pooled, c(0,1)))
  lims[pp,] <- c(min(parameters.pooled) - frac.ran*ran , max(parameters.pooled) + frac.ran*ran)
pp <- 5 # xi0
  parameters.pooled <- c(deoptim.all$gpd5[,5], deoptim.all$gpd6[,5])
  ran <- diff(quantile(parameters.pooled, c(0,1)))
  lims[pp,] <- c(min(parameters.pooled) - frac.ran*ran , max(parameters.pooled) + frac.ran*ran)
pp <- 6 # xi1
  parameters.pooled <- c(deoptim.all$gpd6[,6])
  ran <- diff(quantile(parameters.pooled, c(0,1)))
  lims[pp,] <- c(min(parameters.pooled) - frac.ran*ran , max(parameters.pooled) + frac.ran*ran)



pdf(paste(plot.dir,'priors_normalgamma_SOM.pdf',sep=''), height=7, width=10, colormodel='cmyk')
par(mfrow=c(4,6))
##=============================
model <- 'gpd3'
pp <- 1; par(mai=c(.25,.59,.25,.01))
box.width <- diff(lims[1,])/nbins; box.edges <- seq(from=lims[1,1], to=lims[1,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[1,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[1,1], to=lims[1,2], length.out=1000); lines(x.tmp, dgamma(x=x.tmp, shape=priors_normalgamma[[model]]$lambda$shape, rate=priors_normalgamma[[model]]$lambda$rate), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
mtext('Probability density', side=2, line=1.2, cex=1)
mtext('ST', side=2, line=3, cex=1)
plot.new()
pp <- 2; par(mai=c(.25,.3,.25,.3))
box.width <- diff(lims[3,])/nbins; box.edges <- seq(from=lims[3,1], to=lims[3,2], by=box.width)
hist(log(deoptim.all[[model]][,pp]), xlim=lims[3,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
rate.tmp <- median(log(deoptim.all[[model]][,2])) / (0.5*(max(log(deoptim.all[[model]][,2]))-min(log(deoptim.all[[model]][,2]))))^2
shape.tmp <- median(log(deoptim.all[[model]][,2])) * rate.tmp
x.tmp <- seq(from=lims[3,1], to=lims[3,2], length.out=1000); lines(x.tmp, dgamma(x=x.tmp, shape=shape.tmp, rate=rate.tmp), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
plot.new()
pp <- 3; par(mai=c(.25,.3,.25,.3))
box.width <- diff(lims[5,])/nbins; box.edges <- seq(from=lims[5,1], to=lims[5,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[5,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[5,1], to=lims[5,2], length.out=1000); lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[model]]$xi$mean, sd=priors_normalgamma[[model]]$xi$sd), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
plot.new()
##=============================
model <- 'gpd4'
pp <- 1; par(mai=c(.3,.59,.2,.01))
box.width <- diff(lims[1,])/nbins; box.edges <- seq(from=lims[1,1], to=lims[1,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[1,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[1,1], to=lims[1,2], length.out=1000); lines(x.tmp, dgamma(x=x.tmp, shape=priors_normalgamma[[model]]$lambda0$shape, rate=priors_normalgamma[[model]]$lambda0$rate), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
mtext('Probability density', side=2, line=1.2, cex=1)
mtext('NS1', side=2, line=3, cex=1)
pp <- 2; par(mai=c(.3,.35,.2,.25))
box.width <- diff(lims[2,])/nbins; box.edges <- seq(from=lims[2,1], to=lims[2,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[2,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[2,1], to=lims[2,2], length.out=1000); lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[model]]$lambda1$mean, sd=priors_normalgamma[[model]]$lambda1$sd), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
pp <- 3; par(mai=c(.3,.3,.2,.3))
box.width <- diff(lims[3,])/nbins; box.edges <- seq(from=lims[3,1], to=lims[3,2], by=box.width)
hist(log(deoptim.all[[model]][,pp]), xlim=lims[3,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
rate.tmp <- median(log(deoptim.all[[model]][,3])) / (0.5*(max(log(deoptim.all[[model]][,3]))-min(log(deoptim.all[[model]][,3]))))^2
shape.tmp <- median(log(deoptim.all[[model]][,3])) * rate.tmp
x.tmp <- seq(from=lims[3,1], to=lims[3,2], length.out=1000); lines(x.tmp, dgamma(x=x.tmp, shape=shape.tmp, rate=rate.tmp), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
plot.new()
pp <- 4; par(mai=c(.3,.3,.2,.3))
box.width <- diff(lims[5,])/nbins; box.edges <- seq(from=lims[5,1], to=lims[5,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[5,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[5,1], to=lims[5,2], length.out=1000); lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[model]]$xi$mean, sd=priors_normalgamma[[model]]$xi$sd), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
plot.new()
##=============================
model <- 'gpd5'
pp <- 1; par(mai=c(.35,.59,.15,.01))
box.width <- diff(lims[1,])/nbins; box.edges <- seq(from=lims[1,1], to=lims[1,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[1,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[1,1], to=lims[1,2], length.out=1000); lines(x.tmp, dgamma(x=x.tmp, shape=priors_normalgamma[[model]]$lambda0$shape, rate=priors_normalgamma[[model]]$lambda0$rate), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
mtext('Probability density', side=2, line=1.2, cex=1)
mtext('NS2', side=2, line=3, cex=1)
pp <- 2; par(mai=c(.35,.35,.15,.25))
box.width <- diff(lims[2,])/nbins; box.edges <- seq(from=lims[2,1], to=lims[2,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[2,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[2,1], to=lims[2,2], length.out=1000); lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[model]]$lambda1$mean, sd=priors_normalgamma[[model]]$lambda1$sd), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
pp <- 3; par(mai=c(.35,.3,.15,.3))
box.width <- diff(lims[3,])/nbins; box.edges <- seq(from=lims[3,1], to=lims[3,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[3,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[3,1], to=lims[3,2], length.out=1000); lines(x.tmp, dgamma(x=x.tmp, shape=priors_normalgamma[[model]]$sigma0$shape, rate=priors_normalgamma[[model]]$sigma0$rate), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
pp <- 4; par(mai=c(.35,.3,.15,.3))
box.width <- diff(lims[4,])/nbins; box.edges <- seq(from=lims[4,1], to=lims[4,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[4,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[4,1], to=lims[4,2], length.out=1000); lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[model]]$sigma1$mean, sd=priors_normalgamma[[model]]$sigma1$sd), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
pp <- 5; par(mai=c(.35,.3,.15,.3))
box.width <- diff(lims[5,])/nbins; box.edges <- seq(from=lims[5,1], to=lims[5,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[5,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, xaxt='n', yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[5,1], to=lims[5,2], length.out=1000); lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[model]]$xi$mean, sd=priors_normalgamma[[model]]$xi$sd), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
plot.new()
##=============================
model <- 'gpd6'
pp <- 1; par(mai=c(.5,.59,.01,.01))
box.width <- diff(lims[1,])/nbins; box.edges <- seq(from=lims[1,1], to=lims[1,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[1,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[1,1], to=lims[1,2], length.out=1000); lines(x.tmp, dgamma(x=x.tmp, shape=priors_normalgamma[[model]]$lambda0$shape, rate=priors_normalgamma[[model]]$lambda0$rate), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
mtext('Probability density', side=2, line=1.2, cex=1)
mtext(expression(lambda[0]), side=1, line=2.7, cex=1)
mtext('NS3', side=2, line=3, cex=1)
pp <- 2; par(mai=c(.5,.35,.01,.25))
box.width <- diff(lims[2,])/nbins; box.edges <- seq(from=lims[2,1], to=lims[2,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[2,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[2,1], to=lims[2,2], length.out=1000); lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[model]]$lambda1$mean, sd=priors_normalgamma[[model]]$lambda1$sd), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
mtext(expression(lambda[1]), side=1, line=2.7, cex=1)
pp <- 3; par(mai=c(.5,.3,.01,.3))
box.width <- diff(lims[3,])/nbins; box.edges <- seq(from=lims[3,1], to=lims[3,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[3,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[3,1], to=lims[3,2], length.out=1000); lines(x.tmp, dgamma(x=x.tmp, shape=priors_normalgamma[[model]]$sigma0$shape, rate=priors_normalgamma[[model]]$sigma0$rate), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
mtext(expression(sigma[0]), side=1, line=2.7, cex=1)
pp <- 4; par(mai=c(.5,.3,.01,.3))
box.width <- diff(lims[4,])/nbins; box.edges <- seq(from=lims[4,1], to=lims[4,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[4,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[4,1], to=lims[4,2], length.out=1000); lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[model]]$sigma1$mean, sd=priors_normalgamma[[model]]$sigma1$sd), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
mtext(expression(sigma[1]), side=1, line=2.7, cex=1)
pp <- 5; par(mai=c(.5,.3,.01,.3))
box.width <- diff(lims[5,])/nbins; box.edges <- seq(from=lims[5,1], to=lims[5,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[5,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[5,1], to=lims[5,2], length.out=1000); lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[model]]$xi0$mean, sd=priors_normalgamma[[model]]$xi0$sd), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
mtext(expression(xi[0]), side=1, line=2.7, cex=1)
pp <- 6; par(mai=c(.5,.3,.01,.3))
box.width <- diff(lims[6,])/nbins; box.edges <- seq(from=lims[6,1], to=lims[6,2], by=box.width)
hist(deoptim.all[[model]][,pp], xlim=lims[6,], freq=FALSE, main='', xlab='', ylab='',
     breaks=box.edges, yaxt='n', yaxs='i')
x.tmp <- seq(from=lims[6,1], to=lims[6,2], length.out=1000); lines(x.tmp, dnorm(x=x.tmp, mean=priors_normalgamma[[model]]$xi1$mean, sd=priors_normalgamma[[model]]$xi1$sd), col='red', lwd=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
mtext(expression(xi[1]), side=1, line=2.7, cex=1)
##=============================
dev.off()

#===============================================================================
#


#
#===============================================================================
# SOM FIGURE S2 -- show the calibrated distributions of parameters using the
# (and S3, S4)     normal and gamma priors, with uniform results superimposed
#                  as a dashed curve. Do this for all three sites (S2-S4 figs)

gpd.parameters.uniform <- list.init
for (site in site.names) {
  if(site=='Norfolk') {
    ncdata <- nc_open(filename.norfolk.uniform)
    for (model in types.of.gpd) {gpd.parameters.uniform[[site]]$y90[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd89.',model,sep='')))}
  } else if(site=='Balboa') {ncdata <- nc_open(filename.balboa.uniform)
    for (model in types.of.gpd) {gpd.parameters.uniform[[site]]$y110[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd107.',model,sep='')))}
  } else if(site=='Delfzijl') {ncdata <- nc_open(filename.delfzijl.uniform)
    for (model in types.of.gpd) {gpd.parameters.uniform[[site]]$y130[[model]] <- t(ncvar_get(ncdata, paste('parameters.gpd137.',model,sep='')))}
  } else {print('ERROR - unrecognized site name')}
}

kde.uniform <- kde.normalgamma <- list.init
for (site in site.names) {
  for (model in types.of.gpd) {
    kde.normalgamma[[site]][[all.data[[site]]]][[model]] <- vector('list', length(parnames[[model]])); names(kde.normalgamma[[site]][[all.data[[site]]]][[model]]) <- parnames[[model]]
    for (pp in 1:length(parnames[[model]])) {
      kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]] <- density(gpd.parameters[[site]][[all.data[[site]]]][[model]][,pp])
      kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]] <- density(gpd.parameters.uniform[[site]][[all.data[[site]]]][[model]][,pp])
    }
  }
}

# get limits that bound all of the parameters across both priors and all sites
lims.p <- lims
pp <- 1 # lambda0
lims.p[pp,1] <- min(kde.normalgamma[[1]][[all.data[[1]]]]$gpd3$lambda$x[which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd3$lambda$y > 0)[1]],
                    kde.normalgamma[[1]][[all.data[[1]]]]$gpd4$lambda0$x[which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd4$lambda0$y > 0)[1]],
                    kde.normalgamma[[1]][[all.data[[1]]]]$gpd5$lambda0$x[which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd5$lambda0$y > 0)[1]],
                    kde.normalgamma[[1]][[all.data[[1]]]]$gpd6$lambda0$x[which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd6$lambda0$y > 0)[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd3$lambda$x[which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd3$lambda$y > 0)[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd4$lambda0$x[which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd4$lambda0$y > 0)[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd5$lambda0$x[which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd5$lambda0$y > 0)[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd6$lambda0$x[which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd6$lambda0$y > 0)[1]])
lims.p[pp,2] <- max(kde.normalgamma[[1]][[all.data[[1]]]]$gpd3$lambda$x[rev(which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd3$lambda$y > 0))[1]],
                    kde.normalgamma[[1]][[all.data[[1]]]]$gpd4$lambda0$x[rev(which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd4$lambda0$y > 0))[1]],
                    kde.normalgamma[[1]][[all.data[[1]]]]$gpd5$lambda0$x[rev(which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd5$lambda0$y > 0))[1]],
                    kde.normalgamma[[1]][[all.data[[1]]]]$gpd6$lambda0$x[rev(which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd6$lambda0$y > 0))[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd3$lambda$x[rev(which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd3$lambda$y > 0))[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd4$lambda0$x[rev(which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd4$lambda0$y > 0))[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd5$lambda0$x[rev(which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd5$lambda0$y > 0))[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd6$lambda0$x[rev(which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd6$lambda0$y > 0))[1]])
pp <- 2 # lambda1
lims.p[pp,1] <- min(kde.normalgamma[[1]][[all.data[[1]]]]$gpd4$lambda1$x[which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd4$lambda1$y > 0)[1]],
                    kde.normalgamma[[1]][[all.data[[1]]]]$gpd5$lambda1$x[which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd5$lambda1$y > 0)[1]],
                    kde.normalgamma[[1]][[all.data[[1]]]]$gpd6$lambda1$x[which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd6$lambda1$y > 0)[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd4$lambda1$x[which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd4$lambda1$y > 0)[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd5$lambda1$x[which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd5$lambda1$y > 0)[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd6$lambda1$x[which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd6$lambda1$y > 0)[1]])
lims.p[pp,2] <- max(kde.normalgamma[[1]][[all.data[[1]]]]$gpd4$lambda1$x[rev(which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd4$lambda1$y > 0))[1]],
                    kde.normalgamma[[1]][[all.data[[1]]]]$gpd5$lambda1$x[rev(which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd5$lambda1$y > 0))[1]],
                    kde.normalgamma[[1]][[all.data[[1]]]]$gpd6$lambda1$x[rev(which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd6$lambda1$y > 0))[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd4$lambda1$x[rev(which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd4$lambda1$y > 0))[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd5$lambda1$x[rev(which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd5$lambda1$y > 0))[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd6$lambda1$x[rev(which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd6$lambda1$y > 0))[1]])
pp <- 3 # sigma0
lims.p[pp,1] <- min(log(kde.normalgamma[[1]][[all.data[[1]]]]$gpd3$sigma$x[which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd3$sigma$y > 0)[1]]),
                    log(kde.normalgamma[[1]][[all.data[[1]]]]$gpd4$sigma$x[which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd4$sigma$y > 0)[1]]),
                    kde.normalgamma[[1]][[all.data[[1]]]]$gpd5$sigma0$x[which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd5$sigma0$y > 0)[1]],
                    kde.normalgamma[[1]][[all.data[[1]]]]$gpd6$sigma0$x[which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd6$sigma0$y > 0)[1]],
                    log(kde.normalgamma[[2]][[all.data[[2]]]]$gpd3$sigma$x[which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd3$sigma$y > 0)[1]]),
                    log(kde.normalgamma[[2]][[all.data[[2]]]]$gpd4$sigma$x[which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd4$sigma$y > 0)[1]]),
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd5$sigma0$x[which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd5$sigma0$y > 0)[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd6$sigma0$x[which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd6$sigma0$y > 0)[1]])
lims.p[pp,2] <- max(log(kde.normalgamma[[1]][[all.data[[1]]]]$gpd3$sigma$x[rev(which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd3$sigma$y > 0))[1]]),
                    log(kde.normalgamma[[1]][[all.data[[1]]]]$gpd4$sigma$x[rev(which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd4$sigma$y > 0))[1]]),
                    kde.normalgamma[[1]][[all.data[[1]]]]$gpd5$sigma0$x[rev(which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd5$sigma0$y > 0))[1]],
                    kde.normalgamma[[1]][[all.data[[1]]]]$gpd6$sigma0$x[rev(which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd6$sigma0$y > 0))[1]],
                    log(kde.normalgamma[[2]][[all.data[[2]]]]$gpd3$sigma$x[rev(which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd3$sigma$y > 0))[1]]),
                    log(kde.normalgamma[[2]][[all.data[[2]]]]$gpd4$sigma$x[rev(which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd4$sigma$y > 0))[1]]),
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd5$sigma0$x[rev(which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd5$sigma0$y > 0))[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd6$sigma0$x[rev(which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd6$sigma0$y > 0))[1]])
pp <- 4 # sigma1
lims.p[pp,1] <- min(kde.normalgamma[[1]][[all.data[[1]]]]$gpd5$sigma1$x[which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd5$sigma1$y > 0)[1]],
                    kde.normalgamma[[1]][[all.data[[1]]]]$gpd6$sigma1$x[which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd6$sigma1$y > 0)[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd5$sigma1$x[which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd5$sigma1$y > 0)[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd6$sigma1$x[which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd6$sigma1$y > 0)[1]])
lims.p[pp,2] <- max(kde.normalgamma[[1]][[all.data[[1]]]]$gpd5$sigma1$x[rev(which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd5$sigma1$y > 0))[1]],
                    kde.normalgamma[[1]][[all.data[[1]]]]$gpd6$sigma1$x[rev(which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd6$sigma1$y > 0))[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd5$sigma1$x[rev(which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd5$sigma1$y > 0))[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd6$sigma1$x[rev(which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd6$sigma1$y > 0))[1]])
pp <- 5 # xi0
lims.p[pp,1] <- min(kde.normalgamma[[1]][[all.data[[1]]]]$gpd3$xi$x[which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd3$xi$y > 0)[1]],
                    kde.normalgamma[[1]][[all.data[[1]]]]$gpd4$xi0$x[which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd4$xi0$y > 0)[1]],
                    kde.normalgamma[[1]][[all.data[[1]]]]$gpd5$xi0$x[which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd5$xi0$y > 0)[1]],
                    kde.normalgamma[[1]][[all.data[[1]]]]$gpd6$xi0$x[which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd6$xi0$y > 0)[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd3$xi$x[which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd3$xi$y > 0)[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd4$xi0$x[which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd4$xi0$y > 0)[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd5$xi0$x[which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd5$xi0$y > 0)[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd6$xi0$x[which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd6$xi0$y > 0)[1]])
lims.p[pp,2] <- max(kde.normalgamma[[1]][[all.data[[1]]]]$gpd3$xi$x[rev(which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd3$xi$y > 0))[1]],
                    kde.normalgamma[[1]][[all.data[[1]]]]$gpd4$xi0$x[rev(which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd4$xi0$y > 0))[1]],
                    kde.normalgamma[[1]][[all.data[[1]]]]$gpd5$xi0$x[rev(which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd5$xi0$y > 0))[1]],
                    kde.normalgamma[[1]][[all.data[[1]]]]$gpd6$xi0$x[rev(which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd6$xi0$y > 0))[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd3$xi$x[rev(which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd3$xi$y > 0))[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd4$xi0$x[rev(which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd4$xi0$y > 0))[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd5$xi0$x[rev(which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd5$xi0$y > 0))[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd6$xi0$x[rev(which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd6$xi0$y > 0))[1]])
pp <- 6 # xi1
lims.p[pp,1] <- min(kde.normalgamma[[1]][[all.data[[1]]]]$gpd6$xi1$x[which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd6$xi1$y > 0)[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd6$xi1$x[which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd6$xi1$y > 0)[1]])
lims.p[pp,2] <- max(kde.normalgamma[[1]][[all.data[[1]]]]$gpd6$xi1$x[rev(which(kde.normalgamma[[1]][[all.data[[1]]]]$gpd6$xi1$y > 0))[1]],
                    kde.normalgamma[[2]][[all.data[[2]]]]$gpd6$xi1$x[rev(which(kde.normalgamma[[2]][[all.data[[2]]]]$gpd6$xi1$y > 0))[1]])



#
# Delfzijl posterior results
#

pdf(paste(plot.dir,'posteriors_normalgamma_delfzijl_SOM.pdf',sep=''), height=7, width=10, colormodel='cmyk')
par(mfrow=c(4,6))
##=============================
site <- 'Delfzijl'
model <- 'gpd3'
pp <- 1; par(mai=c(.25,.59,.25,.01))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', axes=FALSE, xlim=lims.p[1,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
mtext('Probability density', side=2, line=1.2, cex=1)
mtext('ST', side=2, line=3, cex=1)

# add a legend:
par(mai=c(.01,.01,.01,.01))
plot.new()
text(0.2, .77, 'Priors:', cex=1.2)
legend(0,0.75, c('normal/gamma','uniform'), lty=c(1,2), cex=1.2, bty='n', lwd=c(2,2))

pp <- 2; par(mai=c(.25,.3,.25,.3))
plot(log(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x), kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', axes=FALSE, xlim=lims.p[3,])
lines(log(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x), kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
plot.new()
pp <- 3; par(mai=c(.25,.3,.25,.3))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', axes=FALSE, xlim=lims.p[5,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
plot.new()
##=============================
model <- 'gpd4'
pp <- 1; par(mai=c(.25,.59,.25,.01))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', axes=FALSE, xlim=lims.p[1,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
mtext('Probability density', side=2, line=1.2, cex=1)
mtext('NS1', side=2, line=3, cex=1)
pp <- 2; par(mai=c(.25,.35,.25,.25))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', axes=FALSE, xlim=lims.p[2,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
pp <- 3; par(mai=c(.25,.3,.25,.3))
plot(log(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x), kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', axes=FALSE, xlim=lims.p[3,])
lines(log(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x), kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
plot.new()
pp <- 4; par(mai=c(.25,.3,.25,.3))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', axes=FALSE, xlim=lims.p[5,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
plot.new()
##=============================
model <- 'gpd5'
pp <- 1; par(mai=c(.25,.59,.25,.01))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', axes=FALSE, xlim=lims.p[1,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
mtext('Probability density', side=2, line=1.2, cex=1)
mtext('NS2', side=2, line=3, cex=1)
pp <- 2; par(mai=c(.25,.35,.25,.25))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', axes=FALSE, xlim=lims.p[2,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
pp <- 3; par(mai=c(.25,.3,.25,.3))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', axes=FALSE, xlim=lims.p[3,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
pp <- 4; par(mai=c(.25,.3,.25,.3))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', axes=FALSE, xlim=lims.p[4,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
pp <- 5; par(mai=c(.25,.3,.25,.3))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', axes=FALSE, xlim=lims.p[5,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
plot.new()
##=============================
model <- 'gpd6'
pp <- 1; par(mai=c(.5,.59,.01,.01))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', yaxt='n', yaxs='i', bty='n', xlim=lims.p[1,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
mtext('Probability density', side=2, line=1.2, cex=1)
mtext(expression(lambda[0]), side=1, line=2.7, cex=1)
mtext('NS3', side=2, line=3, cex=1)
pp <- 2; par(mai=c(.5,.35,.01,.25))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', yaxt='n', yaxs='i', bty='n', xlim=lims.p[2,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
mtext(expression(lambda[1]), side=1, line=2.7, cex=1)
pp <- 3; par(mai=c(.5,.3,.01,.3))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', yaxt='n', yaxs='i', bty='n', xlim=lims.p[3,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
mtext(expression(sigma[0]), side=1, line=2.7, cex=1)
pp <- 4; par(mai=c(.5,.3,.01,.3))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', yaxt='n', yaxs='i', bty='n', xlim=lims.p[4,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
mtext(expression(sigma[1]), side=1, line=2.7, cex=1)
pp <- 5; par(mai=c(.5,.3,.01,.3))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', yaxt='n', yaxs='i', bty='n', xlim=lims.p[5,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
mtext(expression(xi[0]), side=1, line=2.7, cex=1)
pp <- 6; par(mai=c(.5,.3,.01,.3))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', yaxt='n', yaxs='i', bty='n', xlim=lims.p[6,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
mtext(expression(xi[1]), side=1, line=2.7, cex=1)
##=============================
dev.off()


#
# Norfolk posterior results
#

pdf(paste(plot.dir,'posteriors_normalgamma_norfolk_SOM.pdf',sep=''), height=7, width=10, colormodel='cmyk')
par(mfrow=c(4,6))
##=============================
site <- 'Norfolk'
model <- 'gpd3'
pp <- 1; par(mai=c(.25,.59,.25,.01))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', axes=FALSE, xlim=lims.p[1,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
mtext('Probability density', side=2, line=1.2, cex=1)
mtext('ST', side=2, line=3, cex=1)

# add a legend:
par(mai=c(.01,.01,.01,.01))
plot.new()
text(0.2, .77, 'Priors:', cex=1.2)
legend(0,0.75, c('normal/gamma','uniform'), lty=c(1,2), cex=1.2, bty='n', lwd=c(2,2))

pp <- 2; par(mai=c(.25,.3,.25,.3))
plot(log(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x), kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', axes=FALSE, xlim=lims.p[3,])
lines(log(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x), kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
plot.new()
pp <- 3; par(mai=c(.25,.3,.25,.3))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', axes=FALSE, xlim=lims.p[5,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
plot.new()
##=============================
model <- 'gpd4'
pp <- 1; par(mai=c(.25,.59,.25,.01))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', axes=FALSE, xlim=lims.p[1,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
mtext('Probability density', side=2, line=1.2, cex=1)
mtext('NS1', side=2, line=3, cex=1)
pp <- 2; par(mai=c(.25,.35,.25,.25))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', axes=FALSE, xlim=lims.p[2,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
pp <- 3; par(mai=c(.25,.3,.25,.3))
plot(log(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x), kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', axes=FALSE, xlim=lims.p[3,])
lines(log(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x), kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
plot.new()
pp <- 4; par(mai=c(.25,.3,.25,.3))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', axes=FALSE, xlim=lims.p[5,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
plot.new()
##=============================
model <- 'gpd5'
pp <- 1; par(mai=c(.25,.59,.25,.01))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', axes=FALSE, xlim=lims.p[1,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
mtext('Probability density', side=2, line=1.2, cex=1)
mtext('NS2', side=2, line=3, cex=1)
pp <- 2; par(mai=c(.25,.35,.25,.25))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', axes=FALSE, xlim=lims.p[2,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
pp <- 3; par(mai=c(.25,.3,.25,.3))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', axes=FALSE, xlim=lims.p[3,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
pp <- 4; par(mai=c(.25,.3,.25,.3))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', axes=FALSE, xlim=lims.p[4,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
pp <- 5; par(mai=c(.25,.3,.25,.3))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', axes=FALSE, xlim=lims.p[5,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
plot.new()
##=============================
model <- 'gpd6'
pp <- 1; par(mai=c(.5,.59,.01,.01))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', yaxt='n', yaxs='i', bty='n', xlim=lims.p[1,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
mtext('Probability density', side=2, line=1.2, cex=1)
mtext(expression(lambda[0]), side=1, line=2.7, cex=1)
mtext('NS3', side=2, line=3, cex=1)
pp <- 2; par(mai=c(.5,.35,.01,.25))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', yaxt='n', yaxs='i', bty='n', xlim=lims.p[2,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
mtext(expression(lambda[1]), side=1, line=2.7, cex=1)
pp <- 3; par(mai=c(.5,.3,.01,.3))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', yaxt='n', yaxs='i', bty='n', xlim=lims.p[3,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
mtext(expression(sigma[0]), side=1, line=2.7, cex=1)
pp <- 4; par(mai=c(.5,.3,.01,.3))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', yaxt='n', yaxs='i', bty='n', xlim=lims.p[4,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
mtext(expression(sigma[1]), side=1, line=2.7, cex=1)
pp <- 5; par(mai=c(.5,.3,.01,.3))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', yaxt='n', yaxs='i', bty='n', xlim=lims.p[5,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
mtext(expression(xi[0]), side=1, line=2.7, cex=1)
pp <- 6; par(mai=c(.5,.3,.01,.3))
plot(kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.normalgamma[[site]][[all.data[[site]]]][[model]][[pp]]$y,
     type='l', lwd=2, lty=1, main='', xlab='', ylab='', yaxt='n', yaxs='i', bty='n', xlim=lims.p[6,])
lines(kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$x, kde.uniform[[site]][[all.data[[site]]]][[model]][[pp]]$y,
      type='l', lwd=2, lty=2)
u <- par('usr'); arrows(u[1], u[3], u[1], .95*u[4], code=2, length=.15, xpd=TRUE)
lines(c(-1e4,1e4),c(0,0), lty=1, col='black', lwd=2)
mtext(expression(xi[1]), side=1, line=2.7, cex=1)
##=============================
dev.off()


#===============================================================================
#




#
#===============================================================================
# SOM FIGURE 5 - Top row: current surge levels; bottom row: projected 2065 surge
#                levels. Columns: different sites.
#                Projected distributions of 100-year surge level by BMA, and
#                each of the individual model structures.

library(Bolstad)

# all index combinations between the BMA ensemble and each specific model ensemble
# don't even set up 'all.indices', because it is going to be HUGE (and in R,
# therefore slow)
#all.indices <- expand.grid(bma=1:n.ensemble, model=1:n.ensemble)
n.sample <- 1e6 # 1e6 chosen because quantiles seem to stabilize; 1e4 not enough,
                # but 1e5 and 1e6 yield similar results across all 3 sites
sample.indices <- expand.grid(bma=1:n.ensemble, model=1:n.ensemble)[sample(x=(1:(n.ensemble*n.ensemble)), size=n.sample, replace=FALSE),]

# * rl100.mod is to sample from the distribution of [model] (including [bma])
#   for 2016 100-year return levels and 2065 return levels relative to 2016.
# * f.rl100.mod is to fit KDEs for plotting nicely
# * rl100.x are x coordinates for the KDEs, in meters, because we will start using meters after the density estimates
rl100.mod        <- f.rl100.mod        <- vector('list', nsites)
names(rl100.mod) <- names(f.rl100.mod) <- site.names
rl100.dx <- 0.02 # 2 cm sea level is the minimum difference we care about
rl100.x <- seq(-30, 30, by=rl100.dx)

for (site in site.names) {
  rl100.mod[[site]]        <- f.rl100.mod[[site]]        <- vector('list', nmodel)
  names(rl100.mod[[site]]) <- names(f.rl100.mod[[site]]) <- types.of.gpd
  for (model in types.of.gpd) {
    rl100.mod[[site]][[model]]        <- f.rl100.mod[[site]][[model]]        <- vector('list', length(year.names))
    names(rl100.mod[[site]][[model]]) <- names(f.rl100.mod[[site]][[model]]) <- year.names
    rl100.mod[[site]][[model]]$y2016 <- rl100[[site]][[all.data[[site]]]][[model]]$y2016[sample.indices[,2]]
    rl100.mod[[site]][[model]]$y2065 <- rl100[[site]][[all.data[[site]]]][[model]]$y2065[sample.indices[,2]]
    f.rl100.mod[[site]][[model]]$y2016 <- density(x=rl100.mod[[site]][[model]]$y2016, from=1000*min(rl100.x), to=1000*max(rl100.x), n=length(rl100.x), kernel='gaussian')
    f.rl100.mod[[site]][[model]]$y2065 <- density(x=rl100.mod[[site]][[model]]$y2065, from=1000*min(rl100.x), to=1000*max(rl100.x), n=length(rl100.x), kernel='gaussian')
    # normalize
    f.rl100.mod[[site]][[model]]$y2016 <- f.rl100.mod[[site]][[model]]$y2016$y/sintegral(x=rl100.x, fx=f.rl100.mod[[site]][[model]]$y2016$y)$value
    f.rl100.mod[[site]][[model]]$y2065 <- f.rl100.mod[[site]][[model]]$y2065$y/sintegral(x=rl100.x, fx=f.rl100.mod[[site]][[model]]$y2065$y)$value
  }
  rl100.mod[[site]]$bma        <- f.rl100.mod[[site]]$bma        <- vector('list', length(year.names))
  names(rl100.mod[[site]]$bma) <- names(f.rl100.mod[[site]]$bma) <- year.names
  rl100.mod[[site]]$bma$y2016 <- rl100.bma[[site]][[all.data[[site]]]]$y2016[sample.indices[,1]]
  rl100.mod[[site]]$bma$y2065 <- rl100.bma[[site]][[all.data[[site]]]]$y2065[sample.indices[,1]]
  f.rl100.mod[[site]]$bma$y2016 <- density(x=rl100.mod[[site]]$bma$y2016, from=1000*min(rl100.x), to=1000*max(rl100.x), n=length(rl100.x), kernel='gaussian')
  f.rl100.mod[[site]]$bma$y2065 <- density(x=rl100.mod[[site]]$bma$y2065, from=1000*min(rl100.x), to=1000*max(rl100.x), n=length(rl100.x), kernel='gaussian')
  # normalize
  f.rl100.mod[[site]]$bma$y2016 <- f.rl100.mod[[site]]$bma$y2016$y/sintegral(x=rl100.x, fx=f.rl100.mod[[site]]$bma$y2016$y)$value
  f.rl100.mod[[site]]$bma$y2065 <- f.rl100.mod[[site]]$bma$y2065$y/sintegral(x=rl100.x, fx=f.rl100.mod[[site]]$bma$y2065$y)$value
}

# box-whiskers for the plot

## Calculate the quantiles to plot
quantiles.to.grab <- c(0, .05, .25, .5, .75, .95, 1)
quantile.names <- rep(NULL, length(quantiles.to.grab))
for (qq in 1:length(quantiles.to.grab)) {
  if(quantiles.to.grab[qq] >= .10) {
    quantile.names[qq] <- paste('q',100*quantiles.to.grab[qq], sep='')
  } else if(quantiles.to.grab[qq] < .10 & quantiles.to.grab[qq] >= 0) {
    quantile.names[qq] <- paste('q0',100*quantiles.to.grab[qq], sep='')
  }
}

# create matrices for each of the sites
q.rl100.mod <- vector('list', nsites); names(q.rl100.mod) <- site.names
for (site in site.names) {
  q.rl100.mod[[site]] <- mat.or.vec(5, 2*length(quantile.names))
  rownames(q.rl100.mod[[site]]) <- c(types.of.gpd, 'bma')
  colnames(q.rl100.mod[[site]]) <- c(quantile.names, quantile.names)
  for (model in c(types.of.gpd, 'bma')) {
    q.rl100.mod[[site]][model,1:length(quantile.names)] <- 0.001 * quantile(rl100.mod[[site]][[model]]$y2016, quantiles.to.grab)
    q.rl100.mod[[site]][model,(length(quantile.names)+1):(2*length(quantile.names))] <- 0.001 * quantile(rl100.mod[[site]][[model]]$y2065, quantiles.to.grab)
  }
}

## Here we write the CSV file that is the basis for supplemental tables S1-S3
write.csv(x=rbind(q.rl100.mod$Delfzijl, q.rl100.mod$Norfolk), file='../output/returnlevels_threshold99.csv')

#===============================================================================
#



#
#===============================================================================
# SOM FIGURE - BMA weights for POT threshold of 95th and 99.7th quantiles
#              and 1-day declustering time-scale


library(plyr)
library(reshape2)
library(ggplot2)

bma.weights <- readRDS('../output/bma_weights_threshold95.rds')
bma.dfs <- lapply(bma.weights, data.frame, stringsAsFactors = FALSE)
bma.all <- ldply(bma.dfs, .id = 'Site')
bma.all$Model <- rep(c('ST','NS1','NS2','NS3'),2)
bma.all$Model <- ordered(bma.all$Model, levels=c('ST', 'NS1', 'NS2', 'NS3'))
colnames(bma.all)[2:9] <- gsub("X", "", colnames(bma.all)[2:9], fixed=TRUE)
bma.all$ModelType <- ifelse(bma.all$Model == 'ST', 'Stationary', 'Nonstationary')
bma.all$ModelType <- ordered(bma.all$ModelType, levels=c('Stationary', 'Nonstationary'))

bma.melt <- melt(bma.all, id.vars=c('Site','Model', 'ModelType'),variable.name="DataLength", value.name= "ModelWeight")
bma.melt$DataLength <- as.numeric(levels(bma.melt$DataLength)[bma.melt$DataLength])

pdf(paste(plot.dir,'bma_weights_threshold95_SOM.pdf',sep=''), height=3, width=5)
p <- ggplot(bma.melt[!is.na(bma.melt$ModelWeight),]) +
  geom_line(aes(x=DataLength, y=ModelWeight, linetype=ModelType, color=Model)) +
  geom_point(aes(x=DataLength, y=ModelWeight, color=Model), shape=20) +
  facet_grid(Site~.) + theme_bw(base_size=9) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_brewer('Model',type='qual',palette='Set2') +
  scale_linetype_discrete('Model Class') +
  scale_y_continuous('BMA Weight [dimensionless]', breaks=seq(0,0.6,0.1)) +
  scale_x_continuous('Data Length [years]', breaks=seq(30,130,20), expand=c(.01,.01))
print(p)
dev.off()


bma.weights <- readRDS('../output/bma_weights_threshold997.rds')
bma.dfs <- lapply(bma.weights, data.frame, stringsAsFactors = FALSE)
bma.all <- ldply(bma.dfs, .id = 'Site')
bma.all$Model <- rep(c('ST','NS1','NS2','NS3'),2)
bma.all$Model <- ordered(bma.all$Model, levels=c('ST', 'NS1', 'NS2', 'NS3'))
colnames(bma.all)[2:9] <- gsub("X", "", colnames(bma.all)[2:9], fixed=TRUE)
bma.all$ModelType <- ifelse(bma.all$Model == 'ST', 'Stationary', 'Nonstationary')
bma.all$ModelType <- ordered(bma.all$ModelType, levels=c('Stationary', 'Nonstationary'))

bma.melt <- melt(bma.all, id.vars=c('Site','Model', 'ModelType'),variable.name="DataLength", value.name= "ModelWeight")
bma.melt$DataLength <- as.numeric(levels(bma.melt$DataLength)[bma.melt$DataLength])

pdf(paste(plot.dir,'bma_weights_threshold997_SOM.pdf',sep=''), height=3, width=5)
p <- ggplot(bma.melt[!is.na(bma.melt$ModelWeight),]) +
  geom_line(aes(x=DataLength, y=ModelWeight, linetype=ModelType, color=Model)) +
  geom_point(aes(x=DataLength, y=ModelWeight, color=Model), shape=20) +
  facet_grid(Site~.) + theme_bw(base_size=9) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_brewer('Model',type='qual',palette='Set2') +
  scale_linetype_discrete('Model Class') +
  scale_y_continuous('BMA Weight [dimensionless]', breaks=seq(0,0.6,0.1)) +
  scale_x_continuous('Data Length [years]', breaks=seq(30,130,20), expand=c(.01,.01))
print(p)
dev.off()



bma.weights <- readRDS('../output/bma_weights_decl1.rds')
bma.dfs <- lapply(bma.weights, data.frame, stringsAsFactors = FALSE)
bma.all <- ldply(bma.dfs, .id = 'Site')
bma.all$Model <- rep(c('ST','NS1','NS2','NS3'),2)
bma.all$Model <- ordered(bma.all$Model, levels=c('ST', 'NS1', 'NS2', 'NS3'))
colnames(bma.all)[2:9] <- gsub("X", "", colnames(bma.all)[2:9], fixed=TRUE)
bma.all$ModelType <- ifelse(bma.all$Model == 'ST', 'Stationary', 'Nonstationary')
bma.all$ModelType <- ordered(bma.all$ModelType, levels=c('Stationary', 'Nonstationary'))

bma.melt <- melt(bma.all, id.vars=c('Site','Model', 'ModelType'),variable.name="DataLength", value.name= "ModelWeight")
bma.melt$DataLength <- as.numeric(levels(bma.melt$DataLength)[bma.melt$DataLength])

pdf(paste(plot.dir,'bma_weights_decl1_SOM.pdf',sep=''), height=3, width=5)
p <- ggplot(bma.melt[!is.na(bma.melt$ModelWeight),]) +
  geom_line(aes(x=DataLength, y=ModelWeight, linetype=ModelType, color=Model)) +
  geom_point(aes(x=DataLength, y=ModelWeight, color=Model), shape=20) +
  facet_grid(Site~.) + theme_bw(base_size=9) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_brewer('Model',type='qual',palette='Set2') +
  scale_linetype_discrete('Model Class') +
  scale_y_continuous('BMA Weight [dimensionless]', breaks=seq(0,0.6,0.1)) +
  scale_x_continuous('Data Length [years]', breaks=seq(30,130,20), expand=c(.01,.01))
print(p)
dev.off()

#===============================================================================
#



#===============================================================================
# Return period analysis
#===============================================================================

load('returnperiod_analysis_18Apr2018.RData')


lreturn_periods <- log10(return_periods)


# dim(returnlevels.q$Delfzijl$gpd3$y2016) = 100   3


#############
# year 2016 #
year <- 'y2016'
#############
pdf(paste(plot.dir,'return_periods_levels_2016.pdf',sep=''),width=6.5,height=10,colormodel='rgb')

par(mfrow=c(4,2), mai=c(.5,.60,.4,.2))

###  GPD3  ###
model <- 'gpd3'
site <- 'Delfzijl'
plot(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', lwd=2,
     col='seagreen',  ylim=c(3000, 6200), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]][[model]][[year]][,1],rev(returnlevels.q[[site]][[model]][[year]][,3])),
        col=rgb(mycol[3,1],mycol[3,2],mycol[3,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', col='seagreen', lwd=2, lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]]$bma[[year]][,1],rev(returnlevels.q[[site]]$bma[[year]][,3])),
        col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]]$bma[[year]][,2], type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2)
mtext('Return level [m]', side=2, line=2.3, cex=.75);
mtext('Delfzijl', side=3, line=1.5, cex=0.75)
mtext(side=3, text=expression(bold(' a.')), line=.6, cex=0.75, adj=0)
mtext('Return period [years]', side=1, line=2.4, cex=0.75);
xlabs <- c(2, 20, 100, 500)
axis(1, at=log10(xlabs), labels=xlabs, cex.axis=1.1)
axis(2, at=seq(3000, 6500, 500), labels=c('3','','4','','5','','6',''))
legend(log10(2.5), 6000, c('ST', 'BMA'), bty='n',
       fill=c(rgb(mycol[3,1],mycol[3,2],mycol[3,3],.4), rgb(mycol[6,1],mycol[6,2],mycol[6,3],.4)))

site <- 'Norfolk'
plot(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', lwd=2,
     col='seagreen',  ylim=c(1200, 3900), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]][[model]][[year]][,1],rev(returnlevels.q[[site]][[model]][[year]][,3])),
        col=rgb(mycol[3,1],mycol[3,2],mycol[3,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', col='seagreen', lwd=2, lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]]$bma[[year]][,1],rev(returnlevels.q[[site]]$bma[[year]][,3])),
        col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]]$bma[[year]][,2], type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2)
mtext('Return level [m]', side=2, line=2.3, cex=.75);
mtext('Norfolk', side=3, line=1.5, cex=0.75)
mtext(side=3, text=expression(bold(' b.')), line=.6, cex=0.75, adj=0)
mtext('Return period [years]', side=1, line=2.4, cex=0.75);
xlabs <- c(2, 20, 100, 500)
axis(1, at=log10(xlabs), labels=xlabs, cex.axis=1.1)
axis(2, at=seq(1200, 4000, 500), labels=c('','2','','3','','4'))

###  GPD4  ###
model <- 'gpd4'
site <- 'Delfzijl'
plot(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', lwd=2,
     col='seagreen',  ylim=c(3000, 6200), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]][[model]][[year]][,1],rev(returnlevels.q[[site]][[model]][[year]][,3])),
        col=rgb(mycol[3,1],mycol[3,2],mycol[3,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', col='seagreen', lwd=2, lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]]$bma[[year]][,1],rev(returnlevels.q[[site]]$bma[[year]][,3])),
        col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]]$bma[[year]][,2], type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2)
mtext('Return level [m]', side=2, line=2.3, cex=.75);
mtext(side=3, text=expression(bold(' c.')), line=.6, cex=0.75, adj=0)
mtext('Return period [years]', side=1, line=2.4, cex=0.75);
xlabs <- c(2, 20, 100, 500)
axis(1, at=log10(xlabs), labels=xlabs, cex.axis=1.1)
axis(2, at=seq(3000, 6500, 500), labels=c('3','','4','','5','','6',''))
legend(log10(2.5), 6000, c('NS1', 'BMA'), bty='n',
       fill=c(rgb(mycol[3,1],mycol[3,2],mycol[3,3],.4), rgb(mycol[6,1],mycol[6,2],mycol[6,3],.4)))

site <- 'Norfolk'
plot(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', lwd=2,
     col='seagreen',  ylim=c(1200, 3900), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]][[model]][[year]][,1],rev(returnlevels.q[[site]][[model]][[year]][,3])),
        col=rgb(mycol[3,1],mycol[3,2],mycol[3,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', col='seagreen', lwd=2, lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]]$bma[[year]][,1],rev(returnlevels.q[[site]]$bma[[year]][,3])),
        col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]]$bma[[year]][,2], type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2)
mtext('Return level [m]', side=2, line=2.3, cex=.75);
mtext(side=3, text=expression(bold(' d.')), line=.6, cex=0.75, adj=0)
mtext('Return period [years]', side=1, line=2.4, cex=0.75);
xlabs <- c(2, 20, 100, 500)
axis(1, at=log10(xlabs), labels=xlabs, cex.axis=1.1)
axis(2, at=seq(1200, 4000, 500), labels=c('','2','','3','','4'))

###  GPD5  ###
model <- 'gpd5'
site <- 'Delfzijl'
plot(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', lwd=2,
     col='seagreen',  ylim=c(3000, 6200), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]][[model]][[year]][,1],rev(returnlevels.q[[site]][[model]][[year]][,3])),
        col=rgb(mycol[3,1],mycol[3,2],mycol[3,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', col='seagreen', lwd=2, lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]]$bma[[year]][,1],rev(returnlevels.q[[site]]$bma[[year]][,3])),
        col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]]$bma[[year]][,2], type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2)
mtext('Return level [m]', side=2, line=2.3, cex=.75);
mtext(side=3, text=expression(bold(' e.')), line=.6, cex=0.75, adj=0)
mtext('Return period [years]', side=1, line=2.4, cex=0.75);
xlabs <- c(2, 20, 100, 500)
axis(1, at=log10(xlabs), labels=xlabs, cex.axis=1.1)
axis(2, at=seq(3000, 6500, 500), labels=c('3','','4','','5','','6',''))
legend(log10(2.5), 6000, c('NS2', 'BMA'), bty='n',
       fill=c(rgb(mycol[3,1],mycol[3,2],mycol[3,3],.4), rgb(mycol[6,1],mycol[6,2],mycol[6,3],.4)))

site <- 'Norfolk'
plot(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', lwd=2,
     col='seagreen',  ylim=c(1200, 3900), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]][[model]][[year]][,1],rev(returnlevels.q[[site]][[model]][[year]][,3])),
        col=rgb(mycol[3,1],mycol[3,2],mycol[3,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', col='seagreen', lwd=2, lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]]$bma[[year]][,1],rev(returnlevels.q[[site]]$bma[[year]][,3])),
        col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]]$bma[[year]][,2], type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2)
mtext('Return level [m]', side=2, line=2.3, cex=.75);
mtext(side=3, text=expression(bold(' f.')), line=.6, cex=0.75, adj=0)
mtext('Return period [years]', side=1, line=2.4, cex=0.75);
xlabs <- c(2, 20, 100, 500)
axis(1, at=log10(xlabs), labels=xlabs, cex.axis=1.1)
axis(2, at=seq(1200, 4000, 500), labels=c('','2','','3','','4'))

###  GPD6  ###
model <- 'gpd6'
site <- 'Delfzijl'
plot(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', lwd=2,
     col='seagreen',  ylim=c(3000, 6200), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]][[model]][[year]][,1],rev(returnlevels.q[[site]][[model]][[year]][,3])),
        col=rgb(mycol[3,1],mycol[3,2],mycol[3,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', col='seagreen', lwd=2, lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]]$bma[[year]][,1],rev(returnlevels.q[[site]]$bma[[year]][,3])),
        col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]]$bma[[year]][,2], type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2)
mtext('Return level [m]', side=2, line=2.3, cex=.75);
mtext(side=3, text=expression(bold(' g.')), line=.6, cex=0.75, adj=0)
mtext('Return period [years]', side=1, line=2.4, cex=0.75);
xlabs <- c(2, 20, 100, 500)
axis(1, at=log10(xlabs), labels=xlabs, cex.axis=1.1)
axis(2, at=seq(3000, 6500, 500), labels=c('3','','4','','5','','6',''))
legend(log10(2.5), 6000, c('NS3', 'BMA'), bty='n',
       fill=c(rgb(mycol[3,1],mycol[3,2],mycol[3,3],.4), rgb(mycol[6,1],mycol[6,2],mycol[6,3],.4)))

site <- 'Norfolk'
plot(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', lwd=2,
     col='seagreen',  ylim=c(1200, 3900), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]][[model]][[year]][,1],rev(returnlevels.q[[site]][[model]][[year]][,3])),
        col=rgb(mycol[3,1],mycol[3,2],mycol[3,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', col='seagreen', lwd=2, lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]]$bma[[year]][,1],rev(returnlevels.q[[site]]$bma[[year]][,3])),
        col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]]$bma[[year]][,2], type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2)
mtext('Return level [m]', side=2, line=2.3, cex=.75);
mtext(side=3, text=expression(bold(' h.')), line=.6, cex=0.75, adj=0)
mtext('Return period [years]', side=1, line=2.4, cex=0.75);
xlabs <- c(2, 20, 100, 500)
axis(1, at=log10(xlabs), labels=xlabs, cex.axis=1.1)
axis(2, at=seq(1200, 4000, 500), labels=c('','2','','3','','4'))

dev.off()


#############
# year 2065 #
year <- 'y2065'
#############
pdf(paste(plot.dir,'return_periods_levels_2065.pdf',sep=''),width=6.5,height=10,colormodel='rgb')

par(mfrow=c(4,2), mai=c(.5,.60,.4,.2))

###  GPD3  ###
model <- 'gpd3'
site <- 'Delfzijl'
plot(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', lwd=2,
     col='seagreen',  ylim=c(3000, 6200), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]][[model]][[year]][,1],rev(returnlevels.q[[site]][[model]][[year]][,3])),
        col=rgb(mycol[3,1],mycol[3,2],mycol[3,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', col='seagreen', lwd=2, lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]]$bma[[year]][,1],rev(returnlevels.q[[site]]$bma[[year]][,3])),
        col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]]$bma[[year]][,2], type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2)
mtext('Return level [m]', side=2, line=2.3, cex=.75);
mtext('Delfzijl', side=3, line=1.5, cex=0.75)
mtext(side=3, text=expression(bold(' a.')), line=.6, cex=0.75, adj=0)
mtext('Return period [years]', side=1, line=2.4, cex=0.75);
xlabs <- c(2, 20, 100, 500)
axis(1, at=log10(xlabs), labels=xlabs, cex.axis=1.1)
axis(2, at=seq(3000, 6500, 500), labels=c('3','','4','','5','','6',''))
legend(log10(2.5), 6000, c('ST', 'BMA'), bty='n',
       fill=c(rgb(mycol[3,1],mycol[3,2],mycol[3,3],.4), rgb(mycol[6,1],mycol[6,2],mycol[6,3],.4)))

site <- 'Norfolk'
plot(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', lwd=2,
     col='seagreen',  ylim=c(1200, 3900), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]][[model]][[year]][,1],rev(returnlevels.q[[site]][[model]][[year]][,3])),
        col=rgb(mycol[3,1],mycol[3,2],mycol[3,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', col='seagreen', lwd=2, lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]]$bma[[year]][,1],rev(returnlevels.q[[site]]$bma[[year]][,3])),
        col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]]$bma[[year]][,2], type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2)
mtext('Return level [m]', side=2, line=2.3, cex=.75);
mtext('Norfolk', side=3, line=1.5, cex=0.75)
mtext(side=3, text=expression(bold(' b.')), line=.6, cex=0.75, adj=0)
mtext('Return period [years]', side=1, line=2.4, cex=0.75);
xlabs <- c(2, 20, 100, 500)
axis(1, at=log10(xlabs), labels=xlabs, cex.axis=1.1)
axis(2, at=seq(1200, 4000, 500), labels=c('','2','','3','','4'))

###  GPD4  ###
model <- 'gpd4'
site <- 'Delfzijl'
plot(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', lwd=2,
     col='seagreen',  ylim=c(3000, 6200), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]][[model]][[year]][,1],rev(returnlevels.q[[site]][[model]][[year]][,3])),
        col=rgb(mycol[3,1],mycol[3,2],mycol[3,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', col='seagreen', lwd=2, lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]]$bma[[year]][,1],rev(returnlevels.q[[site]]$bma[[year]][,3])),
        col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]]$bma[[year]][,2], type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2)
mtext('Return level [m]', side=2, line=2.3, cex=.75);
mtext(side=3, text=expression(bold(' c.')), line=.6, cex=0.75, adj=0)
mtext('Return period [years]', side=1, line=2.4, cex=0.75);
xlabs <- c(2, 20, 100, 500)
axis(1, at=log10(xlabs), labels=xlabs, cex.axis=1.1)
axis(2, at=seq(3000, 6500, 500), labels=c('3','','4','','5','','6',''))
legend(log10(2.5), 6000, c('NS1', 'BMA'), bty='n',
       fill=c(rgb(mycol[3,1],mycol[3,2],mycol[3,3],.4), rgb(mycol[6,1],mycol[6,2],mycol[6,3],.4)))

site <- 'Norfolk'
plot(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', lwd=2,
     col='seagreen',  ylim=c(1200, 3900), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]][[model]][[year]][,1],rev(returnlevels.q[[site]][[model]][[year]][,3])),
        col=rgb(mycol[3,1],mycol[3,2],mycol[3,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', col='seagreen', lwd=2, lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]]$bma[[year]][,1],rev(returnlevels.q[[site]]$bma[[year]][,3])),
        col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]]$bma[[year]][,2], type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2)
mtext('Return level [m]', side=2, line=2.3, cex=.75);
mtext(side=3, text=expression(bold(' d.')), line=.6, cex=0.75, adj=0)
mtext('Return period [years]', side=1, line=2.4, cex=0.75);
xlabs <- c(2, 20, 100, 500)
axis(1, at=log10(xlabs), labels=xlabs, cex.axis=1.1)
axis(2, at=seq(1200, 4000, 500), labels=c('','2','','3','','4'))

###  GPD5  ###
model <- 'gpd5'
site <- 'Delfzijl'
plot(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', lwd=2,
     col='seagreen',  ylim=c(3000, 6200), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]][[model]][[year]][,1],rev(returnlevels.q[[site]][[model]][[year]][,3])),
        col=rgb(mycol[3,1],mycol[3,2],mycol[3,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', col='seagreen', lwd=2, lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]]$bma[[year]][,1],rev(returnlevels.q[[site]]$bma[[year]][,3])),
        col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]]$bma[[year]][,2], type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2)
mtext('Return level [m]', side=2, line=2.3, cex=.75);
mtext(side=3, text=expression(bold(' e.')), line=.6, cex=0.75, adj=0)
mtext('Return period [years]', side=1, line=2.4, cex=0.75);
xlabs <- c(2, 20, 100, 500)
axis(1, at=log10(xlabs), labels=xlabs, cex.axis=1.1)
axis(2, at=seq(3000, 6500, 500), labels=c('3','','4','','5','','6',''))
legend(log10(2.5), 6000, c('NS2', 'BMA'), bty='n',
       fill=c(rgb(mycol[3,1],mycol[3,2],mycol[3,3],.4), rgb(mycol[6,1],mycol[6,2],mycol[6,3],.4)))

site <- 'Norfolk'
plot(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', lwd=2,
     col='seagreen',  ylim=c(1200, 3900), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]][[model]][[year]][,1],rev(returnlevels.q[[site]][[model]][[year]][,3])),
        col=rgb(mycol[3,1],mycol[3,2],mycol[3,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', col='seagreen', lwd=2, lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]]$bma[[year]][,1],rev(returnlevels.q[[site]]$bma[[year]][,3])),
        col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]]$bma[[year]][,2], type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2)
mtext('Return level [m]', side=2, line=2.3, cex=.75);
mtext(side=3, text=expression(bold(' f.')), line=.6, cex=0.75, adj=0)
mtext('Return period [years]', side=1, line=2.4, cex=0.75);
xlabs <- c(2, 20, 100, 500)
axis(1, at=log10(xlabs), labels=xlabs, cex.axis=1.1)
axis(2, at=seq(1200, 4000, 500), labels=c('','2','','3','','4'))

###  GPD6  ###
model <- 'gpd6'
site <- 'Delfzijl'
plot(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', lwd=2,
     col='seagreen',  ylim=c(3000, 6200), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]][[model]][[year]][,1],rev(returnlevels.q[[site]][[model]][[year]][,3])),
        col=rgb(mycol[3,1],mycol[3,2],mycol[3,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', col='seagreen', lwd=2, lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]]$bma[[year]][,1],rev(returnlevels.q[[site]]$bma[[year]][,3])),
        col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]]$bma[[year]][,2], type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2)
mtext('Return level [m]', side=2, line=2.3, cex=.75);
mtext(side=3, text=expression(bold(' g.')), line=.6, cex=0.75, adj=0)
mtext('Return period [years]', side=1, line=2.4, cex=0.75);
xlabs <- c(2, 20, 100, 500)
axis(1, at=log10(xlabs), labels=xlabs, cex.axis=1.1)
axis(2, at=seq(3000, 6500, 500), labels=c('3','','4','','5','','6',''))
legend(log10(2.5), 6000, c('NS3', 'BMA'), bty='n',
       fill=c(rgb(mycol[3,1],mycol[3,2],mycol[3,3],.4), rgb(mycol[6,1],mycol[6,2],mycol[6,3],.4)))

site <- 'Norfolk'
plot(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', lwd=2,
     col='seagreen',  ylim=c(1200, 3900), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]][[model]][[year]][,1],rev(returnlevels.q[[site]][[model]][[year]][,3])),
        col=rgb(mycol[3,1],mycol[3,2],mycol[3,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]][[model]][[year]][,2], type='l', col='seagreen', lwd=2, lty=3)
polygon(c(lreturn_periods, rev(lreturn_periods)),
        c(returnlevels.q[[site]]$bma[[year]][,1],rev(returnlevels.q[[site]]$bma[[year]][,3])),
        col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.3), border=NA)
lines(lreturn_periods, returnlevels.q[[site]]$bma[[year]][,2], type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2)
mtext('Return level [m]', side=2, line=2.3, cex=.75);
mtext(side=3, text=expression(bold(' h.')), line=.6, cex=0.75, adj=0)
mtext('Return period [years]', side=1, line=2.4, cex=0.75);
xlabs <- c(2, 20, 100, 500)
axis(1, at=log10(xlabs), labels=xlabs, cex.axis=1.1)
axis(2, at=seq(1200, 4000, 500), labels=c('','2','','3','','4'))

dev.off()
#===============================================================================




#===============================================================================
# End
#===============================================================================

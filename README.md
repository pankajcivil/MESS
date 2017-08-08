# MESS: Models for Extreme Storm Surges

Questions? Tony Wong <twong@psu.edu>

## Synopsis

This repository is sandbox for calibrating extreme value models to statistically model the distributions of extreme storm surge events. These models include:
* generalized extreme value distributions (GEV)
* Poisson process/generalized Pareto distributions (PP/GPD)
* [Naveau et al. (2016)](http://onlinelibrary.wiley.com/doi/10.1002/2015WR018552/abstract) / [Papastathopoulos and Tawn (2013)](http://www.sciencedirect.com/science/article/pii/S0378375812002388) models (currently, just type (i) of Naveau et al.)
* others?

## Directory structure

./
   * MESS "home" directory

./calibration/
   * R scripts related to the calibration of the statistical models, including reading data, likelihood functions

./data/
   * tide gauge data for calibration, and temperature data for forcing the non-stationary models
   * also contains fitted RData and RDS objects, for prior distributions  and maximum likelihood parameter estimates (MLEs)

./output/
   * statistical model output (i.e., posterior parameter values, thinned samples from calibrated ensembles, ...)

./R/
   * calibration, plotting and analysis codes

## Motivation

The motivation for these codes is to provide a sandbox for calibrating models using large sets of tide gauge data, and really diving under the hood of some commonly used methods for assessing current storm surge risk as well as projecting future surge levels.

## Installation

To obtain the model codes:
~~~~
git clone https://github.com/tonyewong/MESS.git
~~~~

## Workflow to reproduce the results of Wong et al. (arXiv link once ready)

TODO


## Another Hopefully Useful Code Example: Calibrating GEV models using a long tide gauge record

Suppose you would like to use GEV distributions under a variety of non-stationary model structures (a la [Grinsted et al., 2013](http://www.pnas.org/content/110/14/5369.abstract)) to assess storm surge risk.

1. Checkout the model codes and data.

All of the codes are written assuming that you operate in the R directory. Open an R session.
~~~~
git clone https://github.com/tonyewong/MESS.git
cd MESS/R
R
~~~~

For the paper ([Wong et al., 2017](insert an arXiv link once it is ready)), we create priors by taking all of the tide gauge stations from the [UHSLC database](http://uhslc.soest.hawaii.edu/data/?rq) with at least 90 years of data available. One could also use a spreadsheet tool (provided in this repository as `tg_finder.xls`) to determine all of the tide gauge stations within a center latitude/longitude region around your station of interest, with at least a user-defined number of years of data available. We used this tool to find the ~30 tide gauge stations for the Wong et al. paper, lat/lon limits of 1000 degrees (i.e., anywhere...) and a lower limit of 90 years of data. (Note that as of this writing, Norfolk, Virginia has about 89 years of usable data, because we throw some out due to incomplete coverage.)

Let us proceed under the assumption that you are interested to fit annual block maxima to the Delfzijl, the Netherlands (pronounced approximately as "Delf-shale", with my apologies to my Dutch friends :-)) tide gauge data set. It is a rather unique set of data in that it is long (about 137 years) and continuous (longest gap is 8 hours, and there are only a few gaps).

2. Find and obtain the tide gauge data you would like to use to fit prior distributions.

Preliminary step, before fitting priors to all tide gauge sites from the UHSLC database with at least 30 years of data and within +/-30 degrees lat/lon of Delfzijl, the Netherlands -- First, we interactively use the `tg_finder.xls` tool to find all of the tide gauge sites that fit these criteria. This requires the "extra date functions" add-in for Excel. If you do not have access to Excel (and cannot figure out how to do it in OpenOffice Calc, which admittedly I have put no effort into trying) and/or cannot figure out **how** to do it, feel free to contact Tony at <anthony.e.wong@colorado.edu>. (Seriously, I will be happy to help.)

You can also skip ahead to the fun part because in the `MESS/output/old` directory, there are a few files that are the result of fitting MLEs for four different GEV-style models and four different Naveau-i-style models. They vary from completely stationary ('gev3' and 'nav3' in the codes) to completely non-stationary (all parameters non-stationary), covarying with **global mean temperature normalized to have mean zero from 1901-2000 inclusive** ('gev6' and 'nav6' in the codes). The relevant output files (which should resemble what you would find by running the above MLE calibration scripts) are:
~~~
../output/old/surge_priors_gev_nav_19Jun2017.rds
../output/old/surge_initialvalues_gev_nav_19Jun2017.rds
../output/old/surge_MLEs_gev_nav_19Jun2017.rds
~~~

3. Process the data into a form that can be input to the likelihood functions.

If you did not elect to skip ahead and just use the rds files given above, then the next step is to calibrate (via maximum likelihood) parameter estimates for each candidate model, for each of the tide gauge stations that meets whatever your criteria were from step 2.

If you have already done the processing, you will need to run the `processing_delfzijl.R` and `processing_europe_yearBM.R` routines and then point the `filename.processing` variable in `fit_priors_yearBM_driver.R` to the RData file output from them. If you ahve already run then, then set `filename.processing` in `fit_priors_yearBM_driver.R` to the RData file that resulted from them. Note that these R routines are probably in the `MESS/R/old` directory, because the annual block maxima experiments were not presented in the main paper accompanying this repository.

For these experiments using annual block maxima, the processing amounts to detrending and taking annual block maxima, and potentially throwing out a few years because of poor data coverage. This also constructs a `data_europe` list object, where the top-level list elements correspond to all of the different available tide gauge sites for the fitting. Note that these processing scripts do take a bit of time (more than 5 minutes, less than an hour) because they process for the peaks-over-thresholds (POT)/GPD experiments as well. If you are only interested in the annual block maxima, then you might also be interested to crack open the processing scripts and get rid of all the unnecessary POT stuff.
~~~~
source('read_data_temperature.R')   # need temperature to force the non-stationary models
source('processing_delfzijl.R')     # this should be in the main R directory within the MESS
source('./old/processing_europe.R') # this should be in the R/old directory
~~~~

4. Finally! We are ready to fit some prior distributions.

This routine uses a differential optimization algorithm in order to estimate the maximum likelihood parameters for each type of model, for each of the tide gauge sites. Then, a normal or gamma distribution is fit to the resulting set of MLE parameters (one parameter per tide gauge site). Normal distributions are used for the parameters which should have infinite support, and gamma distributions are used for the parameters with half-infinite support (non-negative, specifically, for example, the GEV scale parameter). This will yield an RDS file (a single R object stored to a file, which you can conveniently read later) that contains the fitted `priors` object. The elements of the `priors` list object are different models; the elements at the model level are the parameter names; the elements at the parameters level are `type` and two other elements that define the distribution. For normally-distributed parameters, these elements are `mean` and `sd`; for gamma-distributed parameters, these elements are `shape` and `rate`.
~~~~
source('./old/fit_priors_yearBM_driver.R')  # this should be in the R/old directory
~~~~

5. Run the MCMC calibration.

Before you run, make sure that you set the relevant input file names (e.g., the prior distribution fits and processed data objects for calibration). To do this, open up the `./old/calibration_yearBM_driver.R` script and modify the following:
~~~~
filename.priors   <-   # put the file you created, or would like to use, here
filename.initvals <-   # put the file you created, or would like to use, here # file holding the 'deoptim.delfzijl' object
filename.mles <-       # put the file you created, or would like to use, here # file holding the 'mle.fits' object
filename.datacalib <-  # put the file you created, or would like to use, here # file holding the 'data_calib' object (calibration data)
~~~~

Also modify to fit your local directory structure in the following line. Make sure we run from your local MESS/R directory.
~~~~
setwd('/home/scrim/axw322/codes/EVT/R')
~~~~

You may also want to modify the settings for the MCMC calibration. These are below. The preliminary chains are done to get a good estimate for the transition covariance matrix for the production calibration, so that the initial estimate is better and the algorithm converges faster. Note that gamma_mcmc represents a trade-off between adapting quickly (and converging faster) and adapting too quickly and getting stuck in a local posterior mode. For the simple models here, multi-modal posteriors are not a problem, at least in my experience. Leave gamma_mcmc at 0.5 unless you have a good reason not to. Also, if you want to do Gelman and Rubin diagnostics for convergence, you should use `nnode_mcmc_prod000` > 1 (the more the merrier!). "Good" values of `niter_mcmc_prod000` are in the 100s of thousands (will probably hack off about 10-50 thousand for burn-in).
~~~~
niter_mcmc_prelim000 <-       # number of MCMC iterations (PRELIMINARY chains)
nnode_mcmc_prelim000 <-       # number of CPUs to use (PRELIMINARY chains)
niter_mcmc_prod000 <-         # number of MCMC iterations (PRODUCTION chains)
nnode_mcmc_prod000 <-         # number of CPUs to use (PRODUCTION chains)
gamma_mcmc000 <-              # speed of adaptation (0.5=faster, 1=slowest)
~~~~

If you do not know how many cores you can run with (`nnode_mcmc_prod000`), you can check the number of *physical* and *logical* cores in your machine, from R, as shown below. Note that the number of physical cores is the number of actual pieces of hardware are in there that are distinct processing units. The number of logical cores accounts for each of the physical cores potentially being able to *hyperthread*, which means a single physical core can have two (or more) *threads*, each of which behaves like a distinct core. But they share computational resources, so a two core, two thread machine will not perform as well as a "true" four-core machine (with four physical cores). Still, even on a 2-by-2 machine allows for 4 parallel MCMC chains, which in turns yields more robust ensemble statistics (i.e., Gelman and Rubin diagnostics).
~~~~
library(parallel)
detectCores(logical = TRUE)      # this is the number of logical cores your machine has
detectCores(logical = FALSE)     # this is the number of physical cores your machine has
~~~~

It is also worth noting here that the use of an adaptive MCMC sampler means that if the transition covariance is changing too much, then you will not technically be sampling from the stationary (posterior) distribution. To be rigorous, one could run the adaptive sampler until convergence (using, say, the Gelman and Rubin or Heidelberger and Welch diagnostics to assess convergence), then restart sampling without adaptation but using the final value from the earlier adaptive sampling as the initial value for the non-adaptive sampling, and using the final transition covariance matrix. But the sampler employed here diminishes the adaptation as the iteration progresses.

Okay. Now that that is out of the way, let's have some fun.

This routine will calibrate each of four types of GEV model and each of four types of Naveau-i model, all using the settings you set above.
~~~~
source('./old/calibration_yearBM_driver.R')
~~~~


## Contributors

Please enjoy the code and offer me any suggestions. I am always happy to chat about how the models could be improved, or the design of new numerical experiments especially.

Questions? Tony Wong (twong@psu.edu)

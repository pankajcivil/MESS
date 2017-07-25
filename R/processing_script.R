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

#
#===============================================================================
# Read tide gauge data for Delfzijl, The Netherlands (site that is the focus of
# the flood risk analysis) and a variety of European tide gauge stations that
# are nearby.
#===============================================================================
#

# Delfzijl
source('processing_delfzijl.R')

# save your work! (done periodically within each script, but you would rather not re-do all that, right?)
save.image(file='../output/preprocessing.RData')

# Other European stations
source('processing_europe.R')

# save your work! (done periodically within each script, but you would rather not re-do all that, right?)
save.image(file='../output/preprocessing.RData')

#
#===============================================================================
# End
#===============================================================================
#

#===============================================================================
# read tide gauge data
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================

read_data <- function(dat.dir, filetype, septype){

  files.tg <- list.files(path=dat.dir, pattern=filetype)

  data <- read.table(paste(dat.dir,files.tg[1],sep=''), header = TRUE, sep=septype)
  if(length(files.tg) > 1) {
    for (ff in 2:length(files.tg)) {
      data <- rbind(data, read.table(paste(dat.dir,files.tg[ff],sep=''), header = TRUE, sep=septype))
    }
  }

  return(data)
}

#===============================================================================
# End
#===============================================================================

#===============================================================================
# fingerprint BRICK SLR to local, normalizing at time.beg, projecting until
# time.end
#
# valid values for rcp: {'RCP26','RCP45','RCP85'}
#
# questions? Tony Wong (twong@psu.edu)
#===============================================================================

source('BRICK_LSL.R')

get_lsl <- function(filename.sealevelrise, rcp, lat, lon, time.beg, time.end, dt) {

  ncdata <- nc_open(filename.sealevelrise)
  slr_lws  <- ncvar_get(ncdata, paste('LWS_',rcp,sep=''))
  slr_te   <- ncvar_get(ncdata, paste('TE_',rcp,sep=''))
  slr_ais  <- ncvar_get(ncdata, paste('AIS_',rcp,sep=''))
  slr_gis  <- ncvar_get(ncdata, paste('GIS_',rcp,sep=''))
  slr_gsic <- ncvar_get(ncdata, paste('GSIC_',rcp,sep=''))
  time_proj <- ncvar_get(ncdata, 'time_proj')
  ens <- ncvar_get(ncdata, 'ens')
  nc_close(ncdata)

  lsl.tmp <- brick_lsl(lat.in=lat, lon.in=lon, slr_gis=slr_gis, slr_gsic=slr_gsic,
                       slr_ais=slr_ais, slr_te=slr_te, slr_lws=slr_lws,
                       n.time=length(time_proj))

  # use just the median?
  lsl <- apply(X=lsl.tmp, MARGIN=1, FUN=median)

  time <- seq(from=time.beg, to=time.end, by=dt)
  lsl <- lsl[which(time_proj==time.beg):which(time_proj==time.end)]

  # normalize to the beginning of the design period
  lsl <- lsl - lsl[1]

  # interpolate linearly to get sub-annual times, if necessary

  # TODO
    # TODO
      # TODO
  if (dt != 1) {}

  lsl.out <- data.frame(cbind(time,lsl))

  return(lsl.out)
}

#===============================================================================
# end
#===============================================================================

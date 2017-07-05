# example plot

library(extRemes)

x.gev <- seq(from=0, to=5000, by=10)
pdf.gev <- devd(x=x.gev, loc=1500, scale=400, shape=0.05, type='GEV')

x.tail <- x.gev[x.gev > 2600]
y.tail <- pdf.gev[x.gev > 2600]

pdf('~/Downloads/gev_probability_density_example.pdf', width=6,height=4,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(.8,.4,.1,.2))

plot(x.gev/1000, pdf.gev, type='l', axes=FALSE, xlim=c(0,4), ylim=c(0,1e-3),
     xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', col='black', lwd=2)
u <- par("usr")
arrows(0, u[3],0, .95*u[4], code = 2, length=.15, xpd = TRUE)
mtext('Probability density', side=2, line=0.1, cex=1);
mtext('Surge level [m]', side=1, line=2.3, cex=1);
axis(1,seq(0,10), cex.axis=1.3)
polygon(c(x.tail, rev(x.tail))/1000, c(y.tail, rep(0, length(y.tail))), col=rgb(1,0,0,.6), border=NA)

arrows(3.2, 4e-4, 2.8, 0.7*pdf.gev[which(x.gev==2800)], length=0.15, angle=30, lty=1, lwd=3, code=2)
text(3.45,5e-4,"flood\        \nprobability",cex=1.1,srt=0)

lines(c(x.tail[1]/1000, x.tail[1]/1000), c(7.05e-4, 0), lty=1, lwd=3)
#arrows(x.tail[1]/1000, 6e-4, x.tail[1]/1000, pdf.gev[which(x.gev==x.tail[1])], length=0.15, angle=30, lty=1, lwd=3, code=2)
text(3.1, 8.2e-4,"levee height relative\     \nto local mean sea level",cex=1.1,srt=0)

dev.off()





# ooh. also make a plot of the disintegration onset timing.
library(ncdf4)

filename.brick.gamma    = '/Users/axw322/codes/BRICK/output_model/BRICK-fastdyn_physical_gamma_01Jun2017.nc'

ncdata <- nc_open(filename.brick.gamma)
  time.proj = ncvar_get(ncdata, 'time_proj')
  gsl.rcp26 = ncvar_get(ncdata, 'GlobalSeaLevel_RCP26')
  gsl.rcp45 = ncvar_get(ncdata, 'GlobalSeaLevel_RCP45')
  gsl.rcp85 = ncvar_get(ncdata, 'GlobalSeaLevel_RCP85')
  ais.rcp26 = ncvar_get(ncdata, 'AIS_RCP26')
  ais.rcp45 = ncvar_get(ncdata, 'AIS_RCP45')
  ais.rcp85 = ncvar_get(ncdata, 'AIS_RCP85')
  gsl.nofd.rcp26 = ncvar_get(ncdata, 'GlobalSeaLevel_nofd_RCP26')
  gsl.nofd.rcp45 = ncvar_get(ncdata, 'GlobalSeaLevel_nofd_RCP45')
  gsl.nofd.rcp85 = ncvar_get(ncdata, 'GlobalSeaLevel_nofd_RCP85')
nc_close(ncdata)

n.ensemble <- ncol(gsl.rcp26)

disint.rcp26 <- gsl.rcp26 - gsl.nofd.rcp26
disint.rcp45 <- gsl.rcp45 - gsl.nofd.rcp45
disint.rcp85 <- gsl.rcp85 - gsl.nofd.rcp85

onset.rcp26 <- rep(NA, n.ensemble)
onset.rcp45 <- rep(NA, n.ensemble)
onset.rcp85 <- rep(NA, n.ensemble)

for (i in 1:n.ensemble) {
  i.disint <- which(disint.rcp26[,i] > 0)
  if(length(i.disint) > 0) {onset.rcp26[i] <- time.proj[i.disint[1]]}
  i.disint <- which(disint.rcp45[,i] > 0)
  if(length(i.disint) > 0) {onset.rcp45[i] <- time.proj[i.disint[1]]}
  i.disint <- which(disint.rcp85[,i] > 0)
  if(length(i.disint) > 0) {onset.rcp85[i] <- time.proj[i.disint[1]]}

}

# kernel smoothing (don't do RCP2.6 because not enough of those disintegrate to get a good estimate of pdf)
x.onset <- seq(from=2000, to=2200, by=1)
ks.rcp45 <- density(x=onset.rcp45[which(!is.na(onset.rcp45))], from=x.onset[1], to=max(x.onset))
ks.rcp85 <- density(x=onset.rcp85[which(!is.na(onset.rcp85))], from=x.onset[1], to=max(x.onset), bw=2.8)

# and plotting

pdf('~/Downloads/ais_disint_onset_timing_pdfs.pdf', width=6,height=4,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(.8,.4,.1,.32))

plot(ks.rcp45$x, ks.rcp45$y, type='l', lwd=2, lty=2, axes=FALSE, xlab='', ylab='',
     xlim=c(2000,2200), ylim=c(0,0.04), xaxt='n', yaxt='n', yaxs='i', xaxs='i')
lines(ks.rcp85$x, ks.rcp85$y, lwd=2, lty=1)
lines(c(2065,2065), c(-1000,1000), lty=1, lwd=2, col='red')
text(2077, 0.038, '2065', cex=1.3)
u <- par("usr")
arrows(2000, u[3], 2000, .95*u[4], code = 2, length=.15, xpd = TRUE)
mtext('Probability density', side=2, line=0.4, cex=1.3);
mtext('Disintegration onset year', side=1, line=2.5, cex=1.3);
axis(1,seq(2000,2200,50), cex.axis=1.3)

legend(2130, .04, c('RCP4.5','RCP8.5'), lty=c(2,1), lwd=2, bty='n', cex=1.3)

dev.off()

#
#
#
#
#

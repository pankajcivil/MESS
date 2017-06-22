#===============================================================================
# plotting
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================

plot.dir <- '../figures/'

#===============================================================================
# get some colorblind-friendly colors to plot with
source('colorblindPalette.R')

#===============================================================================



#===============================================================================
# preliminary plots of projections... what's going on?
#=====================================================

# stationary case
xtmp <- seq(from=0, to=10000, by=10)
mle.gev3 <- amcmc_out$gev3[[1]]$samples[which.max(amcmc_out$gev3[[1]]$log.p),]
mle.nav3 <- amcmc_out$nav3[[1]]$samples[which.max(amcmc_out$nav3[[1]]$log.p),]
pdf.nav3 <- naveau_pdf(x=xtmp, kappa=exp(mle.nav3[1]), sigma=mle.nav3[2], xi=mle.nav3[3])
pdf.gev3 <- devd(x=xtmp, loc=mle.gev3[1], scale=mle.gev3[2], shape=mle.gev3[3])

plot(xtmp, pdf.gev3, type='l', col='blue'); lines(xtmp, pdf.nav3, type='l', col='red')

# above should verify that the stationary case looks the same (check this below)

# gev/nav-4 nonstationary, look at 2040, 2060, 2080 and 2100
mle.gev4 <- amcmc_out$gev4[[1]]$samples[which.max(amcmc_out$gev4[[1]]$log.p),]
mle.nav4 <- amcmc_out$nav4[[1]]$samples[which.max(amcmc_out$nav4[[1]]$log.p),]
time_beg <- 2000
time_end <- 2100
time_proj <- time_forc[which(time_forc==time_beg):which(time_forc==time_end)]
temperature_proj <- temperature_forc[which(time_forc==time_beg):which(time_forc==time_end)]
par.gev4 <- project_gev(parameters=mle.gev4, parnames=parnames_all$gev4, auxiliary=temperature_proj)
par.nav4 <- project_naveau(parameters=mle.nav4, parnames=parnames_all$nav4, auxiliary=temperature_proj)
par.tmp <- par.gev4[which(time_proj==2020),]; pdf.gev4.2020 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.gev4[which(time_proj==2040),]; pdf.gev4.2040 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.gev4[which(time_proj==2060),]; pdf.gev4.2060 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.gev4[which(time_proj==2080),]; pdf.gev4.2080 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.gev4[which(time_proj==2100),]; pdf.gev4.2100 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.nav4[which(time_proj==2020),]; pdf.nav4.2020 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])
par.tmp <- par.nav4[which(time_proj==2040),]; pdf.nav4.2040 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])
par.tmp <- par.nav4[which(time_proj==2060),]; pdf.nav4.2060 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])
par.tmp <- par.nav4[which(time_proj==2080),]; pdf.nav4.2080 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])
par.tmp <- par.nav4[which(time_proj==2100),]; pdf.nav4.2100 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])


# gev/nav-5 nonstationary, look at 2040, 2060, 2080 and 2100
mle.gev5 <- amcmc_out$gev5[[1]]$samples[which.max(amcmc_out$gev5[[1]]$log.p),]
mle.nav5 <- amcmc_out$nav5[[1]]$samples[which.max(amcmc_out$nav5[[1]]$log.p),]
time_beg <- 2000
time_end <- 2100
time_proj <- time_forc[which(time_forc==time_beg):which(time_forc==time_end)]
temperature_proj <- temperature_forc[which(time_forc==time_beg):which(time_forc==time_end)]
par.gev5 <- project_gev(parameters=mle.gev5, parnames=parnames_all$gev5, auxiliary=temperature_proj)
par.nav5 <- project_naveau(parameters=mle.nav5, parnames=parnames_all$nav5, auxiliary=temperature_proj)
par.tmp <- par.gev5[which(time_proj==2020),]; pdf.gev5.2020 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.gev5[which(time_proj==2040),]; pdf.gev5.2040 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.gev5[which(time_proj==2060),]; pdf.gev5.2060 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.gev5[which(time_proj==2080),]; pdf.gev5.2080 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.gev5[which(time_proj==2100),]; pdf.gev5.2100 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.nav5[which(time_proj==2020),]; pdf.nav5.2020 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])
par.tmp <- par.nav5[which(time_proj==2040),]; pdf.nav5.2040 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])
par.tmp <- par.nav5[which(time_proj==2060),]; pdf.nav5.2060 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])
par.tmp <- par.nav5[which(time_proj==2080),]; pdf.nav5.2080 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])
par.tmp <- par.nav5[which(time_proj==2100),]; pdf.nav5.2100 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])


# gev/nav-6 nonstationary, look at 2040, 2060, 2080 and 2100
mle.gev6 <- amcmc_out$gev6[[1]]$samples[which.max(amcmc_out$gev6[[1]]$log.p),]
mle.nav6 <- amcmc_out$nav6[[1]]$samples[which.max(amcmc_out$nav6[[1]]$log.p),]
time_beg <- 2000
time_end <- 2100
time_proj <- time_forc[which(time_forc==time_beg):which(time_forc==time_end)]
temperature_proj <- temperature_forc[which(time_forc==time_beg):which(time_forc==time_end)]
par.gev6 <- project_gev(parameters=mle.gev6, parnames=parnames_all$gev6, auxiliary=temperature_proj)
par.nav6 <- project_naveau(parameters=mle.nav6, parnames=parnames_all$nav6, auxiliary=temperature_proj)
par.tmp <- par.gev6[which(time_proj==2020),]; pdf.gev6.2020 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.gev6[which(time_proj==2040),]; pdf.gev6.2040 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.gev6[which(time_proj==2060),]; pdf.gev6.2060 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.gev6[which(time_proj==2080),]; pdf.gev6.2080 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.gev6[which(time_proj==2100),]; pdf.gev6.2100 <- devd(x=xtmp, loc=par.tmp[1], scale=par.tmp[2], shape=par.tmp[3])
par.tmp <- par.nav6[which(time_proj==2020),]; pdf.nav6.2020 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])
par.tmp <- par.nav6[which(time_proj==2040),]; pdf.nav6.2040 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])
par.tmp <- par.nav6[which(time_proj==2060),]; pdf.nav6.2060 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])
par.tmp <- par.nav6[which(time_proj==2080),]; pdf.nav6.2080 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])
par.tmp <- par.nav6[which(time_proj==2100),]; pdf.nav6.2100 <- naveau_pdf(x=xtmp, kappa=par.tmp[1], sigma=par.tmp[2], xi=par.tmp[3])

par(mfrow=c(3,2))
plot(xtmp/1000, pdf.gev4.2020, type='l', xlim=c(0,10), ylim=c(0, 1e-3), axes=FALSE,
     xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', col='black', main='GEV, location nonstationary')
  lines(xtmp/1000, pdf.gev4.2040, col='blue'); lines(xtmp/1000, pdf.gev4.2060, col='purple');
  lines(xtmp/1000, pdf.gev4.2080, col='red');  lines(xtmp/1000, pdf.gev4.2100, col='orange');
  u <- par("usr")
  arrows(0, u[3],0, .95*u[4], code = 2, length=.15, xpd = TRUE)
  mtext('Probability density', side=2, line=0.1, cex=1);
  mtext('Surge level [m]', side=1, line=2.3, cex=1);
  axis(1,seq(0,10), cex.axis=1.3)
plot(xtmp/1000, pdf.nav4.2020, type='l', xlim=c(0,10), ylim=c(0, 1e-3),axes=FALSE,
     xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', col='black', main='Naveau, lower tail nonstationary')
  lines(xtmp/1000, pdf.nav4.2040, col='blue'); lines(xtmp/1000, pdf.nav4.2060, col='purple');
  lines(xtmp/1000, pdf.nav4.2080, col='red');  lines(xtmp/1000, pdf.nav4.2100, col='orange');
  u <- par("usr")
  arrows(0, u[3],0, .95*u[4], code = 2, length=.15, xpd = TRUE)
  mtext('Probability density', side=2, line=0.1, cex=1);
  mtext('Surge level [m]', side=1, line=2.3, cex=1);
  axis(1,seq(0,10), cex.axis=1.3)
plot(xtmp/1000, pdf.gev5.2020, type='l', xlim=c(0,10), ylim=c(0, 2e-3),axes=FALSE,
     xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', col='black', main='GEV, location and scale nonstationary')
  lines(xtmp/1000, pdf.gev5.2040, col='blue'); lines(xtmp/1000, pdf.gev5.2060, col='purple');
  lines(xtmp/1000, pdf.gev5.2080, col='red');  lines(xtmp/1000, pdf.gev5.2100, col='orange');
  u <- par("usr")
  arrows(0, u[3],0, .95*u[4], code = 2, length=.15, xpd = TRUE)
  mtext('Probability density', side=2, line=0.1, cex=1);
  mtext('Surge level [m]', side=1, line=2.3, cex=1);
  axis(1,seq(0,10), cex.axis=1.3)
plot(xtmp/1000, pdf.nav5.2020, type='l', xlim=c(0,10), ylim=c(0, 1e-3),axes=FALSE,
     xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', col='black', main='Naveau, lower tail and scale nonstationary')
  lines(xtmp/1000, pdf.nav5.2040, col='blue'); lines(xtmp/1000, pdf.nav5.2060, col='purple');
  lines(xtmp/1000, pdf.nav5.2080, col='red');  lines(xtmp/1000, pdf.nav5.2100, col='orange');
  u <- par("usr")
  arrows(0, u[3],0, .95*u[4], code = 2, length=.15, xpd = TRUE)
  mtext('Probability density', side=2, line=0.1, cex=1);
  mtext('Surge level [m]', side=1, line=2.3, cex=1);
  axis(1,seq(0,10), cex.axis=1.3)

plot(xtmp/1000, pdf.gev6.2020, type='l', xlim=c(0,10), ylim=c(0, 3e-3),axes=FALSE,
     xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', col='black', main='GEV, all nonstationary')
  lines(xtmp/1000, pdf.gev6.2040, col='blue'); lines(xtmp/1000, pdf.gev6.2060, col='purple');
  lines(xtmp/1000, pdf.gev6.2080, col='red');  lines(xtmp/1000, pdf.gev6.2100, col='orange');
  u <- par("usr")
  arrows(0, u[3],0, .95*u[4], code = 2, length=.15, xpd = TRUE)
  mtext('Probability density', side=2, line=0.1, cex=1);
  mtext('Surge level [m]', side=1, line=2.3, cex=1);
  axis(1,seq(0,10), cex.axis=1.3)
plot(xtmp/1000, pdf.nav6.2020, type='l', xlim=c(0,10), ylim=c(0, 1e-3),axes=FALSE,
     xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i', col='black', main='Naveau, all nonstationary')
  lines(xtmp/1000, pdf.nav6.2040, col='blue'); lines(xtmp/1000, pdf.nav6.2060, col='purple');
  lines(xtmp/1000, pdf.nav6.2080, col='red');  lines(xtmp/1000, pdf.nav6.2100, col='orange');
  u <- par("usr")
  arrows(0, u[3],0, .95*u[4], code = 2, length=.15, xpd = TRUE)
  mtext('Probability density', side=2, line=0.1, cex=1);
  mtext('Surge level [m]', side=1, line=2.3, cex=1);
  axis(1,seq(0,10), cex.axis=1.3)

#===============================================================================



#===============================================================================
# FIGURE
#
# how does the 2000-year return level change throughout 2000-2100?
# want the 1:2000 event, ensemble 5%, 50%, and 95% quantiles

time_beg <- 2000
time_end <- 2100
time_proj <- time_forc[which(time_forc==time_beg):which(time_forc==time_end)]
temperature_proj <- temperature_forc[which(time_forc==time_beg):which(time_forc==time_end)]

protection.target <- 1/2000

# save the quantiles for the ensemble for each year
return_period_target <- vector('list', length(types.of.model)); names(return_period_target) <- types.of.model
for (model in types.of.model) {
  return_period_target[[model]] <- mat.or.vec(length(time_proj), 3)
  colnames(return_period_target[[model]]) <- c('q5','q50','q95')
}


for (model in types.of.gev) {
  print(paste('starting projections for ',model,' now...',sep=''))
  pb <- txtProgressBar(min=0,max=length(time_proj),initial=0,style=3)
  for (t in 1:length(time_proj)) {
    parameters_project <- t(sapply(1:nrow(parameters[[model]]), function(i) {project_gev(parameters=parameters[[model]][i,], parnames=parnames_all[[model]], auxiliary=temperature_proj[t])}))
    colnames(parameters_project) <- c('mu','sigma','xi')
    level_target <- as.numeric(sapply(1:nrow(parameters_project), function(i) {qevd(p=1-protection.target, loc=parameters_project[i,'mu'], scale=parameters_project[i,'sigma'], shape=parameters_project[i,'xi'])}))
    return_period_target[[model]][t,] <- quantile(level_target, c(.05, .5, .95))
    setTxtProgressBar(pb, t)
  }
  close(pb)
}
for (model in types.of.nav) {
  print(paste('starting projections for ',model,' now...',sep=''))
  pb <- txtProgressBar(min=0,max=length(time_proj),initial=0,style=3)
  for (t in 1:length(time_proj)) {
    parameters_project <- t(sapply(1:nrow(parameters[[model]]), function(i) {project_naveau(parameters=parameters[[model]][i,], parnames=parnames_all[[model]], auxiliary=temperature_proj[t])}))
    colnames(parameters_project) <- c('kappa','sigma','xi')
    level_target <- as.numeric(sapply(1:nrow(parameters_project), function(i) {naveau_invcdf(q=1-protection.target, kappa=parameters_project[i,'kappa'], sigma=parameters_project[i,'sigma'], xi=parameters_project[i,'xi'])}))
    return_period_target[[model]][t,] <- quantile(level_target, c(.05, .5, .95), na.rm=TRUE)
    setTxtProgressBar(pb, t)
  }
  close(pb)
}

# convert from mm to m
for (model in types.of.model) {return_period_target[[model]] <- return_period_target[[model]]/1000}

#
# The actual figure
#

pdf(paste(plot.dir,'2000yReturnPeriod_projections.pdf',sep=''),width=8,height=3.5,colormodel='cmyk')
par(mfrow=c(1,2))
par(mai=c(.65,.65,.20,.2))
model <- 'gev3'
plot(time_proj, return_period_target[[model]][,'q50'], type='l', lwd=2, lty=1,
     ylim=c(4.5,14.5), xlab='', ylab='', xaxt='n', yaxt='n', col='black', xaxs='i')
polygon(c(time_proj,rev(time_proj)), c(return_period_target[[model]][,'q5'], rev(return_period_target[[model]][,'q95'])),
        col=rgb(.5,.5,.5,.5), border=NA)
mtext('1/2000 surge level [m]', side=2, line=2.3, cex=1);
mtext('Year', side=1, line=2, cex=1);
axis(1,seq(2000,2100,by=20), cex.axis=1)
axis(2,seq(5,13, by=1), label=c('5','','7','','9','','11','','13'), cex.axis=1)

model <- 'gev4'; ic <- 6
lines(time_proj, return_period_target[[model]][,'q50'], type='l', lwd=2, lty=1,
      col=rgb(mycol[ic,1], mycol[ic,2], mycol[ic,3]) )
polygon(c(time_proj,rev(time_proj)), c(return_period_target[[model]][,'q5'], rev(return_period_target[[model]][,'q95'])),
        col=rgb(mycol[ic,1], mycol[ic,2], mycol[ic,3], 0.5), border=NA)

model <- 'gev5'; ic <- 9
lines(time_proj, return_period_target[[model]][,'q50'], type='l', lwd=2, lty=1,
      col=rgb(mycol[ic,1], mycol[ic,2], mycol[ic,3]) )
polygon(c(time_proj,rev(time_proj)), c(return_period_target[[model]][,'q5'], rev(return_period_target[[model]][,'q95'])),
        col=rgb(mycol[ic,1], mycol[ic,2], mycol[ic,3], 0.5), border=NA)

legend(2000, 15, c('stationary','location nonstationary', 'location, scale nonstationary', '5-95% credible range'),
       lty=c(1,1,1,NA), pch=c(NA,NA,NA,15), lwd=c(2,2,2,10), cex=1, col=c('black', mycol.rgb[6], mycol.rgb[9], rgb(.5,.5,.5)), bty='n' )

# Too crazy ... need to check?
#model <- 'gev6'; ic <- 7
#lines(time_proj, return_period_target[[model]][,'q50'], type='l', lwd=2, lty=1,
#      col=rgb(mycol[ic,1], mycol[ic,2], mycol[ic,3]) )
#polygon(c(time_proj,rev(time_proj)), c(return_period_target[[model]][,'q5'], rev(return_period_target[[model]][,'q95'])),
#        col=rgb(mycol[ic,1], mycol[ic,2], mycol[ic,3], 0.5), border=NA)

par(mai=c(.65,.65,.20,.2))
model <- 'nav3'
plot(time_proj, return_period_target[[model]][,'q50'], type='l', lwd=2, lty=1,
     ylim=c(4.5,14.5), xlab='', ylab='', xaxt='n', yaxt='n', col='black', xaxs='i')

model <- 'nav6'; ic <- 14
lines(time_proj, return_period_target[[model]][,'q50'], type='l', lwd=2, lty=1,
      col=rgb(mycol[ic,1], mycol[ic,2], mycol[ic,3]) )
polygon(c(time_proj,rev(time_proj)), c(return_period_target[[model]][,'q5'], rev(return_period_target[[model]][,'q95'])),
        col=rgb(mycol[ic,1], mycol[ic,2], mycol[ic,3], 0.5), border=NA)

polygon(c(time_proj,rev(time_proj)), c(return_period_target[[model]][,'q5'], rev(return_period_target[[model]][,'q95'])),
        col=rgb(.5,.5,.5,.5), border=NA)
mtext('1/2000 surge level [m]', side=2, line=2.3, cex=1);
mtext('Year', side=1, line=2, cex=1);
axis(1,seq(2000,2100,by=20), cex.axis=1)
axis(2,seq(5,13, by=1), label=c('5','','7','','9','','11','','13'), cex.axis=1)

model <- 'nav4'; ic <- 6
lines(time_proj, return_period_target[[model]][,'q50'], type='l', lwd=2, lty=1,
      col=rgb(mycol[ic,1], mycol[ic,2], mycol[ic,3]) )
polygon(c(time_proj,rev(time_proj)), c(return_period_target[[model]][,'q5'], rev(return_period_target[[model]][,'q95'])),
        col=rgb(mycol[ic,1], mycol[ic,2], mycol[ic,3], 0.5), border=NA)

model <- 'nav5'; ic <- 9
lines(time_proj, return_period_target[[model]][,'q50'], type='l', lwd=2, lty=1,
      col=rgb(mycol[ic,1], mycol[ic,2], mycol[ic,3]) )
polygon(c(time_proj,rev(time_proj)), c(return_period_target[[model]][,'q5'], rev(return_period_target[[model]][,'q95'])),
        col=rgb(mycol[ic,1], mycol[ic,2], mycol[ic,3], 0.5), border=NA)

legend(2000, 15, c('stationary','lower tail nonstationary', 'lower tail, scale nonstationary', 'both tails, scale nonstationary', '5-95% credible range'),
       lty=c(1,1,1,1,NA), pch=c(NA,NA,NA,NA,15), lwd=c(2,2,2,2,10), cex=1, col=c('black', mycol.rgb[6], mycol.rgb[9], mycol.rgb[14], rgb(.5,.5,.5)), bty='n' )

dev.off()
#===============================================================================



#===============================================================================
# End
#===============================================================================

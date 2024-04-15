nBin = 20

bgmu2 = bgmE2 = bgmZ2 = kcsu2 = kcsE2 = kcsZ2 =
  matrix(NA, ncol= length(setList), nrow = nBin)
bgmu2_ref = bgmE2_ref = bgmZ2_ref =
  kcsu2_ref = kcsE2_ref = kcsZ2_ref = c()
for(i in seq_along(setList)) {
  D2 = dataList[[paste0(setList[i],'_cal')]]
  uE = D2$uE
  io = order(uE, decreasing = FALSE)
  uE = uE[io]
  E  = D2$E[io]
  Z  = E / uE
  M  = length(uE)

  bgmu2_ref[i] = skewgm(uE^2)
  bgmE2_ref[i] = skewgm(E^2)
  bgmZ2_ref[i] = skewgm(Z^2)
  kcsu2_ref[i] = kurtcs(uE^2)
  kcsE2_ref[i] = kurtcs(E^2)
  kcsZ2_ref[i] = kurtcs(Z^2)

  intrv = ErrViewLib::genIntervals(1:M, nBin)
  for(j in 1:intrv$nbr) {
    sel      = intrv$lwindx[j]:intrv$upindx[j]
    bgmu2[j,i] = skewgm(uE[sel]^2)
    bgmE2[j,i] = skewgm(E[sel]^2)
    bgmZ2[j,i] = skewgm(Z[sel]^2)
    kcsu2[j,i] = kurtcs(uE[sel]^2)
    kcsE2[j,i] = kurtcs(E[sel]^2)
    kcsZ2[j,i] = kurtcs(Z[sel]^2)
  }
}

nS = length(setList)
png(
  file = file.path(figDir, paste0('fig_08.png')),
  width  = 2*gPars$reso,
  height = 3*gPars$reso
)
par(
  mfrow = c(3,2),
  mar = c(3,2,2,0), #gPars$mar,
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)
boxplot(bgmu2, outcex=0.3,
        xlab = 'Set #',
        ylab = expression(beta[GM]({u[E]}^2)),
        main = '' )
grid()
abline(h=0.6, lty = 2, col = gPars$cols[2])
boxplot(bgmu2, outcex=0.3, col = gPars$cols_tr2[5], add = TRUE)
points(1:nS, bgmu2_ref, pch = 16, col = gPars$cols[2])
box()

boxplot(kcsu2, outcex=0.3,
        xlab = 'Set #',
        ylim = c(-2,6),
        ylab = expression(kappa[cs]({u[E]}^2)))
grid()
abline(h=3, lty = 2, col = gPars$cols[2])
boxplot(kcsu2, outcex=0.3, col = gPars$cols_tr2[5], add = TRUE)
points(1:nS, kcsu2_ref, pch = 16, col = gPars$cols[2])
box()

boxplot(bgmE2, outcex=0.3,
        xlab = 'Set #',
        ylab = expression(beta[GM](E^2)))
grid()
abline(h=0.8, lty = 2, col = gPars$cols[2])
boxplot(bgmE2, outcex=0.3, col = gPars$cols_tr2[5], add = TRUE)
points(1:nS, bgmE2_ref, pch = 16, col = gPars$cols[2])
box()

boxplot(kcsE2, outcex=0.3,
        xlab = 'Set #',
        ylim = c(0,24),
        ylab = expression(kappa[cs](E^2)))
grid()
abline(h=5, lty = 2, col = gPars$cols[2])
boxplot(kcsE2, outcex=0.3, col = gPars$cols_tr2[5], add = TRUE)
points(1:nS, kcsE2_ref, pch = 16, col = gPars$cols[2])
box()

boxplot(bgmZ2, outcex=0.3,
        xlab = 'Set #',
        ylab = expression(beta[GM](Z^2)))
grid()
abline(h=0.8, lty = 2, col = gPars$cols[2])
boxplot(bgmZ2, outcex=0.3, col = gPars$cols_tr2[5], add = TRUE)
points(1:nS, bgmZ2_ref, pch = 16, col = gPars$cols[2])
box()

boxplot(kcsZ2, outcex=0.3,
        xlab = 'Set #',
        ylim = c(0,24),
        ylab = expression(kappa[cs](Z^2)))
grid()
abline(h=5, lty = 2, col = gPars$cols[2])
boxplot(kcsZ2, outcex=0.3, col = gPars$cols_tr2[5], add = TRUE)
points(1:nS, kcsZ2_ref, pch = 16, col = gPars$cols[2])
box()

dev.off()

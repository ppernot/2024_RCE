
# Skewness

M = 1e6
nuSeq = 10 ^ (seq(log10(1), log10(100), length.out = 50))
res = res2 = res3 = res4 = res5 = rep(NA_real_,length(nuSeq))
for (i in seq_along(nuSeq)) {
  cat(i, '/')
  nu = nuSeq[i]
  X = rf(M, 1, nu)
  res[i] = ErrViewLib::skewgm(X)
  X = MCMCpack::rinvgamma(M, nu/2, nu/2)
  res2[i] = ErrViewLib::skewgm(X)
  # res3[i] = cv(X)
  # X = rchisq(M, nu)
  # res4[i] = ErrViewLib::skewgm(X)
  # res5[i] = cv(X)
}
cat('\n')

png(
  file = file.path(figDir, paste0('fig_05a.png')),
  width  = 1*gPars$reso,
  height = 1*gPars$reso
)
par(
  mfrow = c(1, 1),
  mar = c(3,2,2,0), #gPars$mar,
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)
plot(
  res2, res, type = 'b', pch = 16, col = gPars$cols[1],
  # xlim = c(0,1), xaxs = 'i',
  xlim = c(0.1,1), xaxs = 'i',
  xlab = expression(beta[GM](u[E]^2)),
  # ylim = c(0,1), yaxs = 'i',
  ylim = c(0.6,1), yaxs = 'i',
  ylab = expression(list(beta[GM](E^2),beta[GM](Z^2))),
  panel.first = grid()
)
polygon(c(0.6,0.6,1,1), c(0,2,2,0), border = NA, col= gPars$cols_tr[2])
polygon(c(0,0,1,1), c(1,0.8,0.8,1), border = NA, col= gPars$cols_tr[2])

for(i in seq_along(setList)) {
  D2 = dataList[[paste0(setList[i],'_cal')]]
  uE = D2$uE
  E = D2$E
  x = ErrViewLib::skewgm(uE^2)
  # y1 = cv(uE)
  y = ErrViewLib::skewgm(E^2)
  z = ErrViewLib::skewgm((E/uE)^2)
  points(x, y, pch = 17, col = gPars$cols[2])
  text(x, y, labels = i, col = gPars$cols[2], pos = 2)
  points(x, z, pch = 15, col = gPars$cols[5])
  text(x, z, labels = i, col = gPars$cols[5], pos = 2)
}
legend(
  'topleft', bty = 'n', cex = 0.85,
  legend = c(
    'NIG model',
    expression(beta[GM](E^2)),
    expression(beta[GM](Z^2))
  ),
  col = gPars$cols[c(1,2,5)],
  pch = c(16, 17, 15),
  lty = c(1,0,0)
)
box()
dev.off()

## Kurtosis ####

M = 1e5
nuSeq = 10 ^ (seq(log10(1), log10(100), length.out = 50))
res = res2 = res3 = res4 = res5 = rep(NA_real_,length(nuSeq))
for (i in seq_along(nuSeq)) {
  cat(i, '/')
  nu = nuSeq[i]
  X = rf(M, 1, nu)
  res[i] = ErrViewLib::kurtcs(X)
  X = MCMCpack::rinvgamma(M, nu/2, nu/2)
  res2[i] = ErrViewLib::kurtcs(X)
  # res3[i] = cv(X)
  # X = rchisq(M, nu)
  # res4[i] = ErrViewLib::kurtcs(X)
  # res5[i] = cv(X)
}
cat('\n')

png(
  file = file.path(figDir, paste0('fig_05b.png')),
  width  = 1*gPars$reso,
  height = 1*gPars$reso
)
par(
  mfrow = c(1, 1),
  mar = c(3,2,2,0), #gPars$mar,
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)
plot(
  res2, res, type = 'b', log = 'y',
  pch = 16, col = gPars$cols[1],
  # xlim = c(0,1), xaxs = 'i',
  xlim = c(-2, 6), xaxs = 'i',
  xlab = expression(kappa[CS](u[E]^2)),
  # ylim = c(0,1), yaxs = 'i',
  ylim = c(1, 40), yaxs = 'i',
  ylab = expression(list(kappa[CS](E^2),kappa[CS](Z^2))),
  panel.first = grid()
)
polygon(c(3,3,10,10), c(0.1,50,50,0.1), border = NA, col= gPars$cols_tr[2])
polygon(c(-2,-2,10,10), c(50,5,5,50), border = NA, col= gPars$cols_tr[2])

for(i in seq_along(setList)) {
  D2 = dataList[[paste0(setList[i],'_cal')]]
  uE = D2$uE
  E = D2$E
  x = ErrViewLib::kurtcs(uE^2)
  # y1 = cv(uE)
  y = ErrViewLib::kurtcs(E^2)
  z = ErrViewLib::kurtcs((E/uE)^2)
  points(x, y, pch = 17, col = gPars$cols[2])
  text(x, y, labels = i, col = gPars$cols[2], pos = 2)
  points(x, z, pch = 15, col = gPars$cols[5])
  text(x, z, labels = i, col = gPars$cols[5], pos = 2)
}
legend(
  'topleft', bty = 'n', cex = 0.85, ncol=1,
  legend = c(
    'NIG model',
    expression(kappa[CS](E^2)),
    expression(kappa[CS](Z^2))
  ),
  col = gPars$cols[c(1,2,5)],
  pch = c(16, 17, 15),
  lty = c(1,0,0)
)
box()
dev.off()

M = 500000

## Skewness bGM ####

nuSeq = 10 ^ (seq(log10(2), log10(50), length.out = 30))
res = res2 = rep(NA_real_,length(nuSeq))
for (i in seq_along(nuSeq)) {
  nu = nuSeq[i]
  X = rf(M, 1, nu)
  res[i] = ErrViewLib::skewgm(X)
  X = MCMCpack::rinvgamma(M, nu/2, nu/2)
  res2[i] = ErrViewLib::skewgm(X)
}
X = rnorm(M)^2
resN = ErrViewLib::skewgm(X)

png(
  file = file.path(figDir, paste0('fig_03a.png')),
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

plot(nuSeq, res, type = 'b', log = 'x',
     pch = 19, col = gPars$cols[5],
     xlab = expression(nu),
     ylab = expression(beta[GM](X^2)), ylim = c(0,1.01), yaxs='i',
     panel.first = grid(equilogs = FALSE),
     panel.last = box())
points(nuSeq, res2, type = 'b',
     pch = 17, col = gPars$cols[2])
abline(h = resN, lty = 2, col = gPars$cols[5])
abline(v = 4, lty = 3, col = gPars$cols[1])
mtext(4, 1, at = 4, col = gPars$cols[1], cex = 0.85*gPars$cex)
legend(
  'topright', bty = 'n', cex = 0.85,
  legend = c(expression(X^2%~%F(1,nu)),
             expression(X^2%~%{Gamma^-1}(nu/2,nu/2))),
  pch = c(19,17),
  col = gPars$cols[c(5,2)]
)
dev.off()

## Kurtosis kCS ####

nuSeq = 10 ^ (seq(log10(2), log10(50), length.out = 30))
res = res2 = rep(NA_real_,length(nuSeq))
for (i in seq_along(nuSeq)) {
  nu = nuSeq[i]
  X = rf(M, 1, nu)
  res[i] = ErrViewLib::kurtcs(X)
  X = MCMCpack::rinvgamma(M, nu/2, nu/2)
  res2[i] = ErrViewLib::kurtcs(X)
}
X = rnorm(M)^2
resN = ErrViewLib::kurtcs(X)

png(
  file = file.path(figDir, paste0('fig_03b.png')),
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

plot(nuSeq, res, type = 'b', log = 'x',
     pch = 19, col = gPars$cols[5],
     xlab = expression(nu),
     ylab = expression(kappa[CS](X^2)),  yaxs='i',
     ylim = c(0,13),
     panel.first = grid(equilogs = FALSE),
     panel.last = box())
points(nuSeq, res2, type = 'b',
       pch = 17, col = gPars$cols[2])
abline(h = resN, lty = 2, col = gPars$cols[5])
abline(v = 4, lty = 3, col = gPars$cols[1])
mtext(4, 1, at = 4, col = gPars$cols[1], cex = 0.85*gPars$cex)
legend(
  'topright', bty = 'n', cex = 0.85,
  legend = c(expression(X^2%~%F(1,nu)),
             expression(X^2%~%{Gamma^-1}(nu/2,nu/2))),
  pch = c(19,17),
  col = gPars$cols[c(5,2)]
)
dev.off()

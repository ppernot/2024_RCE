getCase <- function(case, filterZerE = FALSE) {
  uE = unlist(
    read.csv(
      file.path(
        '..','Data','JAC2024',
        case,
        'model_errors_leaveout_calibrated.csv')
    ),
    use.names = FALSE)
  E  = unlist(
    read.csv(
      file.path(
        '..','Data','JAC2024',
        case,
        'residuals_leaveout.csv')
    ),
    use.names = FALSE)

  # Clean up
  zeruE = uE <= 0
  zerE  = E == 0
  zer = zeruE
  if(filterZerE)
    zer = zer | zerE
  uE = uE[!zer]
  E  = E[!zer]

  return(
    list(
      E = E,
      uE = uE,
      zeruE = zeruE,
      zerE = zerE
    )
  )
}

# stats ####
cases = sort(
  c('debyeT_aflow','dielectric','RPV_TTS','perovskite_conductivity',
          'bandgap_expt','perovskite_Opband','piezoelectric','heusler',
          'concrete','phonon_freq','semiconductor_lvls','double_perovskite_gap',
          'perovskite_workfunction','metallicglass_Rc','metallicglass_Dmax',
          'perovskite_stability','thermal_conductivity','perovskite_Habs',
          'exfoliation_E','steel_yield','perovskite_tec',
          'diffusion','Li_conductivity','metallicglass_Rc_LLM',
          'hea_hardness','oxide_vacancy','perovskite_ASR','Mg_alloy',
          'superconductivity','thermalcond_aflow','thermalexp_aflow',
          'perovskite_formationE','elastic_tensor')
)


if(doCalc) {
  # scores ####
  nBoot = 5000
  cl <- makeCluster(detectCores())
  stats   = names(calScoresBS(1:length(E),cbind(E,uE)))

  scores = bias = cilo = ciup = ciString = zmat =
    matrix(NA, nrow = length(cases), ncol = length(stats))

  targets = c(1,0,1)

  for(i in seq_along(cases)) {
    D2 = getCase(cases[i]); print(cases[i])
    uE = D2$uE
    E  = D2$E
    M  = length(uE)

    # BS scores and CIs
    bs = fPredBS(cbind(E,uE), calScoresBS,
                 nBoot = nBoot, cl = cl)
    scores[i,]   = bs$t0
    bias[i,]     = bs$bias
    cilo[i,]     = bs$bca[1,]
    ciup[i,]     = bs$bca[2,]
    ciString[i,] = paste0('[', signif(bs$bca[1,],3),', ',
                          signif(bs$bca[2,],3), ']')
    zmat[i,]     = fZetaBS(bs, targets)
  }
  stopCluster(cl)
  colnames(cilo) = colnames(ciup) = colnames(scores) = stats
  save(
    stats, scores, bias,
    cilo, ciup, ciString, zmat,
    file = 'Jac2024Scores.Rda'
  )
} else {
  load(file = 'JAC2024Scores.Rda')
}

df = data.frame(
  Name = NA,
  M = NA, nzeruE = NA, nzerE = NA,
  bgmuE2 = NA, bgmE2 = NA, bgmZ2 = NA,
  zms = NA, rce = NA)
for(i in seq_along(cases)) {
  case = cases[i]
  D = getCase(case); print(case)

  uE = D$uE
  E  = D$E
  Z  = E / uE
  M = length(Z)

  df = rbind(
    df,
    c(
      Name   = case,
      M      = M,
      nzeruE = sum(D$zeruE),
      nzerE  = sum(D$zerE),
      bgmuE2 = round(ErrViewLib::skewgm(uE^2),3),
      bgmE2  = round(ErrViewLib::skewgm(E^2),3),
      bgmZ2  = round(ErrViewLib::skewgm(Z^2),3),
      zms    = round(scores[i,1],2),
      rce    = round(scores[i,2],2)
    )
  )
}
df = df[-1,]
rownames(df) = 1:length(cases)

sink(file = file.path(tabDir,'JAC2024Scores.tex'))
print(knitr::kable(df, 'latex'))
sink()

# Fig. 10 ####
png(
  file = file.path(figDir, 'fig_10a.png'),
  width  = 1 * gPars$reso,
  height = 1 * gPars$reso
)
par(
  mfrow = c(2,1),
  mar = c(3,3,0.5,0.2),
  mgp = gPars$mgp,
  pty = 'm',
  tcl = gPars$tcl,
  cex = 0.8 * gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)
x = 1:length(cases)
ylim = c(0.3,1.5)

stat = 'ZMS'
lo = cilo[,stat]
up = ciup[,stat]
val = lo <= 1 & 1 <= up
cols = ifelse(val, 5, 2)
plot(x, scores[,stat],
     pch = 16, col = gPars$cols[cols],
     xlab = 'Set #',
     ylim = ylim, ylab = "ZMS",
     panel.first = {grid(); abline(h = 1, lty = 2)},
     panel.last = box()
)
segments(x, lo, x, up, col = gPars$cols[cols],
         lwd = 1.5*gPars$lwd)
legend(
  'bottomright', bty = 'n',
  legend = c('valid','invalid'),
  pch = 16, lty = 1, lwd = 2*gPars$lwd,
  col = gPars$cols[c(5,2)]
)

ylim = c(-0.5, 0.2)
stat = 'RCE'
lo = cilo[,stat]
up = ciup[,stat]
val = lo <= 0 & 0 <= up
cols = ifelse(val, 5, 2)
plot(x, scores[,stat],
     pch = 16, col = gPars$cols[cols],
     xlab = 'Set #',
     ylim = ylim, ylab = "RCE",
     panel.first = {grid(); abline(h = 0, lty = 2)},
     panel.last = box()
)
segments(x, lo, x, up, col = gPars$cols[cols],
         lwd = 1.5*gPars$lwd)
dev.off()

png(
  file = file.path(figDir, paste0('fig_10b.png')),
  width  = 1*gPars$reso,
  height = 1*gPars$reso
)
par(
  mfrow = c(1, 1),
  mar = gPars$mar,
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)
xlim = c(-0.05, 6.5)
plot(
  abs(zmat[,1]), abs(zmat[,2]),
  pch = 16, col = gPars$cols[1], cex = 0.5,
  xlab = expression(group("|",zeta[ZMS],"|")),
  xlim = xlim, xaxs = 'i',
  ylab = expression(group("|",zeta[RCE],"|")),
  ylim = xlim, yaxs = 'i'
)
polygon(c(0,0,1,1), c(1,10,10,1), border = NA,
        col= gPars$cols_tr2[2])
polygon(c(1,10,10,1), c(0,0,1,1), border = NA,
        col= gPars$cols_tr2[2])
abline(a=0,b=1, lty = 2)
text(abs(zmat[,1]), abs(zmat[,2]), labels = 1:NROW(zmat),
     cex = 0.5, col = gPars$cols[1], pos = 4, offset = 0.25)
dev.off()

# Fig. 9 ####
M = 1e5
nuSeq = 10 ^ (seq(log10(1), log10(100), length.out = 50))
res = res2 = rep(NA_real_,length(nuSeq))
for (i in seq_along(nuSeq)) {
  cat(i, '/')
  nu = nuSeq[i]
  X = rf(M, 1, nu)
  res[i] = ErrViewLib::skewgm(X)
  X = MCMCpack::rinvgamma(M, nu/2, nu/2)
  res2[i] = ErrViewLib::skewgm(X)
}

png(
  file = file.path(figDir, paste0('fig_09.png')),
  width  = 1*gPars$reso,
  height = 1*gPars$reso
)
par(
  mfrow = c(1, 1),
  mar = gPars$mar,
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)
plot(
  res2, res, type = 'b', pch = 16, cex = 0.5, col = gPars$cols[1],
  xlim = c(0.1,1.025), xaxs = 'i',
  xlab = expression(beta[GM](u[E]^2)),
  ylim = c(0.6,1.01), yaxs = 'i',
  ylab = expression(list(beta[GM](E^2),beta[GM](Z^2))),
  panel.first = grid()
)
polygon(c(0.6,0.6,1.1,1.1), c(0,2,2,0), border = NA, col= gPars$cols_tr[2])
polygon(c(0,0,1.1,1.1), c(1.1,0.8,0.8,1.1), border = NA, col= gPars$cols_tr[2])

for(i in seq_along(cases)) {
  x = as.numeric(df$bgmuE2[i])
  y = as.numeric(df$bgmE2[i])
  z = as.numeric(df$bgmZ2[i])
  points(x, y, pch = 17, col = gPars$cols[2], cex = 0.5)
  text(x, y, labels = i, col = gPars$cols[2], cex = 0.65, pos = 2)
  points(x, z, pch = 15, col = gPars$cols[5], cex = 0.5)
  text(x, z, labels = i, col = gPars$cols[5], cex = 0.65, pos = 1)
}
legend(
  'topleft', bty = 'n', cex = 0.65,
  legend = c(
    'NIG model',
    expression(beta[GM](E^2)),
    expression(beta[GM](Z^2))
  ),
  col = gPars$cols[c(1,2,5)],
  pch = c(16, 17, 15),
  lty = c(1,0,0)
)
dev.off()


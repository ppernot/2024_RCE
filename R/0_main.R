figDir = '../Figs'
tabDir = '../Tabs'

library(nptest)
library(parallel)
library(ErrViewLib)
gPars = ErrViewLib::setgPars(type = 'publish')

# Load functions ####
source('functions.R')

# Load datasets ####
source('getData.R')

# Results ####
nBoot =  10000

## Data summary and properties ####
cl <- makeCluster(detectCores())
stats   = names(calScoresBS1(1:length(E),cbind(E,uE)))

scores = bias = cilo = ciup = ciString = zmat =
  matrix(NA, nrow = length(setList), ncol = length(stats))
size =
  betaGM_uE2 = betaGM_E2 = betaGM_Z2 =
  kappaCS_uE2 = kappaCS_E2 = kappaCS_Z2 =
  cv_uE = cv_uE2 = c()

targets = c(1,0)

for(i in seq_along(setList)) {
  D2 = dataList[[paste0(setList[i],'_cal')]]; print(setList[i])
  uE = D2$uE
  E  = D2$E
  M  = length(uE)
  size[i] = M

  cv_uE[i]       = cv(uE)
  cv_uE2[i]      = cv(uE^2)
  betaGM_uE2[i]  = ErrViewLib::skewgm(uE^2)
  betaGM_E2[i]   = ErrViewLib::skewgm(E^2)
  betaGM_Z2[i]   = ErrViewLib::skewgm((E/uE)^2)
  kappaCS_uE2[i] = ErrViewLib::kurtcs(uE^2)
  kappaCS_E2[i]  = ErrViewLib::kurtcs(E^2)
  kappaCS_Z2[i]  = ErrViewLib::kurtcs((E/uE)^2)

  # BS scores and CIs
  bs = fPredBS(cbind(E,uE), calScoresBS1, nBoot = nBoot, cl = cl)
  scores[i,]   = bs$t0
  bias[i,]     = bs$bias
  cilo[i,]     = bs$bca[1,]
  ciup[i,]     = bs$bca[2,]
  ciString[i,] = paste0('[', signif(bs$bca[1,],3),', ',
                        signif(bs$bca[2,],3), ']')
  zmat[i,]     = fZetaBS(bs, targets)
}
stopCluster(cl)

X = round(
  cbind(
    cv_uE, cv_uE2,
    betaGM_uE2, betaGM_E2, betaGM_Z2,
    kappaCS_uE2, kappaCS_E2, kappaCS_Z2
  ),
  2
)
sink(file = file.path(tabDir,'betagm.tex'))
print(knitr::kable(X, 'latex'))
sink()

save(
  stats, setList, scores, bias,
  cilo, ciup, ciString, zmat,
  cv_uE, cv_uE2,
  betaGM_uE2, betaGM_E2, betaGM_Z2,
  kappaCS_uE2, kappaCS_E2, kappaCS_Z2,
  file = 'scores.Rda'
)
# load('scores.Rda')

## Table 2 ####
for(i in seq_along(stats)) {
  df = data.frame(
    set      = 1:length(setList),
    stat     = signif(scores[,i],3),
    bias     = signif(bias[,i],2),
    CI       = ciString[,i],
    zeta     = round(zmat[,i],2)
  )
  print(knitr::kable(df))
  sink(file = file.path(tabDir,paste0('tab',stats[i],'.tex')))
  print(knitr::kable(df, 'latex'))
  sink()
}

## Figs. 1 & 2 / Table 1 ####

source('fitDist.R')

## Fig. 3 ####

source('beta_kappa.R')

## Fig. 4 ####

### Fig. 4a ####
doCalc = FALSE
if(doCalc) {
  # Takes about 1/2 day to complete on 8 cores...
  source("testValidRCE.R")
} else {
  load('testValidRCE.Rda')
}

pro = cilo = ciup = matrix(NA,nrow=length(dfList),ncol = 2)
for(j in seq_along(dfList)) {
  nu = dfList[j]
  # zm = zList[[paste0(nu)]]
  tm = tList[[paste0(nu)]]
  pro[j,]  = colMeans(tm)
  success  = colMeans(tm) * nTry
  trials   = c(nTry, nTry)
  ci       = DescTools::BinomCI(success, trials, method = "wilsoncc")
  cilo[j,] = ci[,2]
  ciup[j,] = ci[,3]
}

png(
  file = file.path(figDir, paste0('fig_04a.png')),
  width  = gPars$reso,
  height = gPars$reso
)
par(
  mfrow = c(1, 1),
  mar = c(3, 3, 4, 0.5), #gPars$mar,
  mgp = gPars$mgp,
  pty = 'm',
  tcl = gPars$tcl,
  cex = gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)
matplot(dfList, pro, type = 'b',
        pch=16:17, lty =1, col=gPars$cols[1:2],
        log = 'x',
        xlab = expression(nu[IG]), xlim = c(2,10),
        ylab = expression(p[val]), ylim = c(0.6,1.0))
grid(equilogs = FALSE)
abline(h=0.95,lty=2)
for(i in 1:2)
  segments(dfList,cilo[,i],dfList,ciup[,i],col=i)
box()
mtext(paste0('- ',signif(kappaCS,2)),side = 3, at = dfList,
      col = 'blue', cex = 0.75 * par()$cex, las = 2)
mtext(expression(kappa[CS](u[E]^2)), side = 3, at = 1.6,
      padj = 0, col = 'blue', cex = par()$cex)
mtext(paste0(signif(betaGM,2)),side = 3, at = dfList,
      line = 2, col = 'red', cex = 0.75 * par()$cex, las = 2)
mtext(expression(beta[GM](u[E]^2)), side = 3, at = 1.6,
      padj = -1.6, col = 'red', cex = par()$cex)
legend(
  'bottomright', bty = 'n',
  legend = stats,
  lty = 1,
  pch = 16:17,
  col = gPars$cols[1:2]
)
dev.off()

### Fig. 4b ####
doCalc = FALSE
if(doCalc) {
  # Takes about 1/2 day to complete on 8 cores...
  source("testValidRCE2.R")
} else {
  load('testValidRCE2.Rda')
}

pro = cilo = ciup = matrix(NA,nrow=length(dfList),ncol = 2)
for(j in seq_along(dfList)) {
  nu = dfList[j]
  # zm = zList[[paste0(nu)]]
  tm = tList[[paste0(nu)]]
  pro[j,]  = colMeans(tm)
  success  = colMeans(tm) * nTry
  trials   = c(nTry, nTry)
  ci       = DescTools::BinomCI(success, trials, method = "wilsoncc")
  cilo[j,] = ci[,2]
  ciup[j,] = ci[,3]
}

png(
  file = file.path(figDir, paste0('fig_04b.png')),
  width  = gPars$reso,
  height = gPars$reso
)
par(
  mfrow = c(1, 1),
  mar = c(3, 3, 4, 0.5), #gPars$mar,
  mgp = gPars$mgp,
  pty = 'm',
  tcl = gPars$tcl,
  cex = gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)
matplot(dfList, pro, type = 'b',
        pch=16:17, lty =1, col=gPars$cols[1:2],
        log = 'x',
        xlab = expression(nu[D]), xlim = c(2,10),
        ylab = expression(p[val]), ylim = c(0.6,1.0))
grid(equilogs = FALSE)
abline(h=0.95,lty=2)
for(i in 1:2)
  segments(dfList,cilo[,i],dfList,ciup[,i],col=i)
box()
mtext(paste0('- ',signif(kappaCS,2)),side = 3, at = dfList,
      col = 'blue', cex = 0.75 * par()$cex, las = 2)
mtext(expression(kappa[CS](u[E]^2)), side = 3, at = 1.6,
      padj = 0, col = 'blue', cex = par()$cex)
mtext(paste0(signif(betaGM,2)),side = 3, at = dfList,
      line = 2, col = 'red', cex = 0.75 * par()$cex, las = 2)
mtext(expression(beta[GM](u[E]^2)), side = 3, at = 1.6,
      padj = -1.6, col = 'red', cex = par()$cex)
legend(
  'bottomright', bty = 'n',
  legend = stats,
  lty = 1,
  pch = 16:17,
  col = gPars$cols[1:2]
)
dev.off()

## Fig. 5 ####
# Might take some time...
source('beta_kappa_Sets.R')


## Fig. 6a ####
png(
  file = file.path(figDir, paste0('fig_06a.png')),
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

xlim = c(-0.1,4)# range(abs(zmat))
plot(abs(zmat[,1]), abs(zmat[,2]),
     type = 'p', pch = 16, cex=0.75, col = gPars$cols[1],
     xlab = expression(group("|",zeta[ZMS],"|")),
     xlim = xlim, xaxs = 'i',
     ylab = expression(group("|",zeta[RCE],"|")),
     ylim = xlim, yaxs = 'i')
grid()
polygon(c(-1,-1,1,1), c(1,6,6,1), border = NA, col= gPars$cols_tr2[2])
polygon(c(1,6,6,1), c(-1,-1,1,1), border = NA, col= gPars$cols_tr2[2])
# abline(h=1, lty=2)
# abline(v=1, lty = 2)
abline(a=0, b=1, lty = 2)
text(abs(zmat[,1]), abs(zmat[,2]),1:nrow(zmat),
     pos = 3, adj = 0, offset = 0.25)
box()

dev.off()

## Fig. 6b ####

cl <- makeCluster(detectCores())
stats   = names(calScoresBS1(1:length(E),cbind(E,uE)))

scores = bias = cilo = ciup = ciString = zmat =
  matrix(NA, nrow = length(setList), ncol = length(stats))
size = betagm = c()

targets = c(1,0)

for(i in seq_along(setList)) {
  D2 = dataList[[paste0(setList[i],'_cal')]]; print(setList[i])
  uE = D2$uE
  E  = D2$E
  M  = length(uE)
  size[i] = M

  # Remove top 5% uE
  io = order(uE, decreasing = TRUE)
  uE = uE[io]
  E  = E[io]
  M  = length(E)
  t5 = round(M*0.05)
  uE = uE[(t5+1):M]
  E  = E[(t5+1):M]
  M  = length(uE)


  betagm[i] = ErrViewLib::skewgm(uE)

  # BS scores and CIs
  bs = fPredBS(cbind(E,uE), calScoresBS1, nBoot = nBoot, cl = cl)
  scores[i,]   = bs$t0
  bias[i,]     = bs$bias
  cilo[i,]     = bs$bca[1,]
  ciup[i,]     = bs$bca[2,]
  ciString[i,] = paste0('[', signif(bs$bca[1,],3),', ',
                        signif(bs$bca[2,],3), ']')
  zmat[i,]     = fZetaBS(bs, targets)
}
stopCluster(cl)

png(
  file = file.path(figDir, paste0('fig_06b.png')),
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

xlim = c(-0.1,4)# range(abs(zmat))
plot(abs(zmat[,1]), abs(zmat[,2]),
     type = 'p', pch = 16, cex=0.75, col = gPars$cols[1],
     xlab = expression(group("|",zeta[ZMS],"|")),
     xlim = xlim, xaxs = 'i',
     ylab = expression(group("|",zeta[RCE],"|")),
     ylim = xlim, yaxs = 'i')
grid()
polygon(c(-1,-1,1,1), c(1,6,6,1), border = NA, col= gPars$cols_tr2[2])
polygon(c(1,6,6,1), c(-1,-1,1,1), border = NA, col= gPars$cols_tr2[2])
# abline(h=1, lty=2)
# abline(v=1, lty = 2)
abline(a=0, b=1, lty = 2)
text(abs(zmat[,1]), abs(zmat[,2]),1:nrow(zmat),
     pos = 3, adj = 0, offset = 0.25)
box()

dev.off()

## Fig. 7 ####

png(
  file = file.path(figDir, paste0('fig_07.png')),
  width  = 3*gPars$reso,
  height = 3*gPars$reso
)
par(
  mfrow = c(3,3),
  mar = c(3,2,2,0), #gPars$mar,
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = 1.5*gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)
cl <- makeCluster(detectCores())
for(i in seq_along(setList)) {
  D2 = dataList[[paste0(setList[i],'_cal')]]
  uE = D2$uE
  E  = D2$E
  M = length(E)
  Z  = E / uE

  pcVec = seq(0,10,by=0.1)
  io = order(uE, E, decreasing = TRUE)
  X = cbind(E[io],uE[io])
  ccRCE = Sconf0(X , RCEfun, pcVec)
  ccRCE = ccRCE-ccRCE[1]
  ccZMS = Sconf0(X , ZMSfun, pcVec)
  ccZMS = ccZMS-ccZMS[1]

  ciZMS = c(cilo[i,1]-scores[i,1],ciup[i,1]-scores[i,1])
  ciRCE = c(cilo[i,2]-scores[i,2],ciup[i,2]-scores[i,2])

  ylim = range(c(ccRCE,ccZMS,ciRCE,ciZMS))
  plot(
    pcVec, ccZMS, log = '', type = 'b',
    col = gPars$cols[5], lwd = 2*gPars$lwd,
    pch = 16, cex = 0.5,
    xlab = 'k% discarded', xlim = c(-1,10),
    ylab = expression(group("",list(Delta[ZMS],Delta[RCE]),"")),
    ylim = ylim,
    main = paste0('Set ', i)
  )
  grid(equilogs = FALSE)
  points(pcVec, ccRCE, col = gPars$cols[2], type = 'b',
         pch = 16, lwd = 2*gPars$lwd, cex = 0.5)
  abline(h=c(0,1), lty =2, col = gPars$cols[1])

  segments(-1, ciRCE[1], -1, ciRCE[2],
           lwd=40, lend=2, col=gPars$cols_tr2[2])
  segments(-0.5, ciZMS[1], -0.5, ciZMS[2],
           lwd=40, lend=2, col=gPars$cols_tr2[5])

  if(i==1)
    legend(
      'topright', bty = 'n',
      legend = c(expression(Delta[ZMS]),expression(Delta[RCE])),
      pch = 16, lty = 0,
      col = gPars$cols[c(5,2)]
    )
}
stopCluster(cl)
dev.off()

## Fig. 8 ####

source('bin_beta_kappa.R')

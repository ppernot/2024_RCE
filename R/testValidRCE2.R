library(nptest)
library(parallel)
library(ErrViewLib)

source('functions.R')

load('testValidRCE2.Rda')

nBoot = 5000
M     = 5000
nTry  = 1000

uE = sqrt(MCMCpack::rinvgamma(M, 3, 3 ))
E  = uE * rnorm(M)
stats = names(calScoresBS1(1:M,cbind(E,uE)))

dfList = c(2.5, 3:10, 20)

if(FALSE) {
  tList = zList = list()
  cl <- makeCluster(detectCores())
  for(j in seq_along(dfList)) {
    # for(j in c(1,10)) {
    nu = dfList[j]
    cat('\n',nu,': ')
    scores = tm = bias = muBS = zmatBS =
      matrix(NA, nrow = nTry, ncol = length(stats))
    for(i in 1:nTry) {
      cat(i,'/ ')
      E = uE * rT4(M, nu)
      X  = cbind(E,uE)
      scores[i,] = calScoresBS1(1:M,X)
      muBS[i,] = c(1, 0)

      bs = nptest::np.boot(
        x=1:M, data = X, statistic = calScoresBS1,
        R = nBoot, level = 0.95, method = "bca",
        parallel = TRUE, cl = cl)

      ci       = bs$bca
      bias[i,] = bs$bias
      tm[i,]   = muBS[i,] >= ci[1,] & muBS[i,] <= ci[2,]

      delta = scores[i,] - muBS[i,]
      lim = (delta >= 0) * ci[1,] + (delta < 0) *ci[2,]
      width = abs(scores[i,] - lim)
      zmatBS[i,] = abs(delta) / width
    }
    zList[[paste0(nu)]] = zmatBS <= 1
    tList[[paste0(nu)]] = tm
    betaGM[j] = mean(tmp)
  }
  stopCluster(cl)
}

# Estimate average skewness for simulated uncertainties
betaGM = kappaCS = tmp = tmp2 = c()
for(j in seq_along(dfList)) {
  nu = dfList[j]
  cat('\n',nu,': ')
  scores = tm = bias = muBS = zmatBS =
    matrix(NA, nrow = nTry, ncol = length(stats))
  for(i in 1:500) {
    cat(i,'/ ')
    E       = uE * rT4(M, nu)
    tmp[i]  = ErrViewLib::skewgm(E^2)
    tmp2[i] = ErrViewLib::kurtcs(E^2)
  }
  cat('\n')
  betaGM[j]  = mean(tmp)
  kappaCS[j] = mean(tmp2)
}

save(nTry, stats, dfList, zList, tList, betaGM, kappaCS,
     file = 'testValidRCE2.Rda')

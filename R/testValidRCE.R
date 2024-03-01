cl <- makeCluster(detectCores())

nBoot = 5000
M     = 5000
nTry  = 1000

dfList = c(2:10,15,20)

uE = sqrt(MCMCpack::rinvgamma(M, 2, 2 ))
E  = uE * rnorm(M)
stats = names(calScoresBS1(1:M,cbind(E,uE)))

tList = zList = list()

for(j in seq_along(dfList)) {
  nu = dfList[j]
  cat('\n',nu,': ')
  scores = tm = bias = muBS = zmatBS =
    matrix(NA, nrow = nTry, ncol = length(stats))
  for(i in 1:nTry) {
    cat(i,'/ ')
    uE = sqrt(MCMCpack::rinvgamma(M, nu/2, nu/2 ))
    E = uE * rnorm(M)
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
}
stopCluster(cl)

# Estimate average skewness for simulated uncertainties
betaGM = tmp = c()
for(j in seq_along(dfList)) {
  nu = dfList[j]
  cat('\n',nu,': ')
  scores = tm = bias = muBS = zmatBS =
    matrix(NA, nrow = nTry, ncol = length(stats))
  for(i in 1:nTry) {
    uE = sqrt(MCMCpack::rinvgamma(M, nu/2, nu/2 ))
    tmp[i] =  ErrViewLib::skewgm(uE)
  }
  betaGM[j] = mean(tmp)
}

save(nTry,stats,dfList,zList,tList,betaGM,file='testValidRCE.Rda')

# Sensitivity of RCE, RCE2 and ZMS to nuIG and nuD

stats = c('RCE','RCE2','1-ZMS')

fStats = function(x, data) {

  # Binning-dependent statistics with unknown targets
  uE = data[x,2]
  E  = data[x,1]

  MV   = mean(uE^2)
  MSE  = mean(E^2)
  RMV  = sqrt(MV)
  RMSE = sqrt(MSE)
  ZMS = mean((E/uE)^2)

  c(
    1 - RMSE/ RMV,
    1 - MSE/ MV,
    1 - ZMS
  )
}

# set.seed(123)

nMC  = 1000
M    = 5000
nuIG = 2*c(1, 1.5, 2, 3, 4, 6, 10)
nuTS = c(3, 4, 5, 6, 10, 100)
nuCases = expand.grid(nuIG, nuTS)

smc = matrix(NA, nrow = nMC, ncol = 3)
scores = uscores = matrix(NA, nrow = nrow(nuCases), ncol = 3)

nuEff = beta = betaZ = c()
for(j in 1:nrow(nuCases)) {
  df1 = nuCases[j,1]
  df2 = nuCases[j,2]
  for(k in 1:nMC) {
    uE = sort(sqrt(MCMCpack::rinvgamma(M,df1/2,df1/2)))
    Ep = uE * rT4(M, df = df2)
    smc[k,] = fStats(1:M, cbind(Ep,uE))
  }
  scores[j,]  = apply(smc, 2, mean, na.rm = TRUE)
  uscores[j,] = apply(smc, 2, sd, na.rm = TRUE)/sqrt(nMC)

  # Estimate nu_eff
  uE = sort(sqrt(MCMCpack::rinvgamma(10*M,df1/2,df1/2)))
  Ep = uE * rT4(10*M, df = df2)
  fit.t = fitdistrplus::fitdist(
    Ep, "t_ls",
    start = list(df = df1,mu = mean(Ep),sigma = sd(Ep)),
    keepdata = FALSE)
  nuEff[j] = signif(summary(fit.t)$estimate[1],2)
  beta[j] = round(ErrViewLib::skewgm(Ep^2),2)
  betaZ[j] = round(ErrViewLib::skewgm((Ep/uE)^2),2)
}
save(stats, scores, uscores, beta, betaZ,
     nuEff, nuIG, nuTS, nuCases,
     file = 'altRCE.Rda')





cl <- makeCluster(detectCores())
stats   = names(calScoresBS1(1:length(E),cbind(E,uE)))

scores = bias = cilo = ciup = ciString = zmat =
  matrix(NA, nrow = length(setList), ncol = length(stats))
size = c()

targets = c(1,0)

for(i in seq_along(setList)) {
  D2 = dataList[[paste0(setList[i],'_cal')]]; print(setList[i])
  uE = D2$uE
  E  = D2$E
  M  = length(uE)
  size[i] = M

  # BS scores and CIs
  bs = fPredBS(cbind(E,uE), calScoresBS1, nBoot = 10000, cl = cl)
  scores[i,]   = bs$t0
  bias[i,]     = bs$bias
  cilo[i,]     = bs$bca[1,]
  ciup[i,]     = bs$bca[2,]
  ciString[i,] = paste0('[', signif(bs$bca[1,],3),', ',
                             signif(bs$bca[2,],3), ']')
  # valid[i,]    = targets >= cilo[i,] & targets <= ciup[i,]
  zmat[i,]     = fZetaBS(bs, targets)
}
stopCluster(cl)

save(
  stats, setList, scores, bias, cilo, ciup, ciString, zmat,
  file = 'scores.Rda'
)


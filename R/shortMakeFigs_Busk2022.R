figDir = '../Figs'
library(ErrViewLib)
library(CHNOSZ)
gPars = ErrViewLib::setgPars(type = 'publish'); gPars$cex = 4.5
scalePoints = 0.2

set.seed(123) # Reproducibility

fractHetero = function(S) {
  fHetero = c()
  for(i in seq_along(S)) {
    compo = CHNOSZ::count.elements(S[i])
    tot   = sum(compo, na.rm = TRUE)
    CandH = sum(compo['C'], compo['H'], na.rm = TRUE)
    fHetero[i] = 1 - CandH / tot
  }
  fHetero
}

testfVal = function(res, p = 0.95) {
  fv  = res$fVal
  lfv = res$lofVal
  ufv = res$upfVal
  msg = paste0('Percent of valid intervals :',signif(fv,2) ,'\n',
      'Lower limit                :',signif(lfv,2),'\n',
      'Upper limit                :',signif(ufv,2),'\n',
      ifelse((lfv-p)*(ufv-p) < 0,
             paste0('Valid  : interval contains ',p),
             paste0('Invalid: interval does not contain ',p)),
      '\n'
  )
  return(msg)
}

# Busk2022 ####
D = read.table(
  '../Data/BUS2022/qm9_U0_test_Orig.csv',
  sep = ',',
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

S  = D[, "formula"]
R  = D[, "U0"]
C  = D[, "prediction"]
uC = D[, "uncertainty_total"] ^ 0.5
E  = R - C
uE = uC

M = length(uE)

unit = '[eV]'

# subsample
if (M > 25000) {
  sel = sample(1:M, 25000)
  C  = C[sel]
  E  = E[sel]
  uE = uE[sel]
}
Z = E / uE

resVarZ = ErrViewLib::varZCI(Z, method = 'cho')
VarZ    = resVarZ$mean
uVarZ   = resVarZ$sd
print(ErrViewLib::prettyUnc(VarZ, uVarZ, numDig = 1))

# Get mass from formula; use as Input feature
masses = CHNOSZ::mass(S)
fHetero = fractHetero(S)
print(signif(cor(cbind(masses,fHetero,abs(E),uE),method = "spearman"),2))

# Figs ####

## Fig 01a ####
nBin = 100
png(
  file = file.path(figDir, 'FigShort_01a.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotEvsPU(
  uE,
  Z,
  nBin = nBin,
  type = 'horiz',
  logX = TRUE,
  runQuant = FALSE,
  runMean = TRUE,
  runMS = TRUE,
  scalePoints = scalePoints,
  ylim = c(-5, 5),
  xlab = "Uncertainty, uE [eV]",
  # ylab = "Error [eV]",
  label = 1,
  legLoc = 'topright',
  gPars = gPars
)
dev.off()

## Fig 01b ####
png(
  file = file.path(figDir, 'FigShort_01b.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotEvsPU(
  masses,
  Z,
  nBin = nBin,
  type = 'horiz',
  logX = FALSE,
  runQuant = FALSE,
  runMean = TRUE,
  runMS = TRUE,
  scalePoints = scalePoints,
  xlim = c(50, 150),
  ylim = c(-5, 5),
  xlab = "Molecular mass, X1 [Da]",
  label = 2,
  showLegend = FALSE,
  gPars = gPars
)
dev.off()

## Fig 01c ####
png(
  file = file.path(figDir, 'FigShort_01c.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotEvsPU(
  fHetero,
  Z,
  nBin = nBin,
  type = 'horiz',
  logX = FALSE,
  runQuant = FALSE,
  runMean = TRUE,
  runMS = TRUE,
  scalePoints = scalePoints,
  xlim = c(-0.05,0.8),
  ylim = c(-5, 5),
  xlab = "Fraction of heteroatoms, X2",
  label = 3,
  showLegend = FALSE,
  gPars = gPars
)
dev.off()

## Fig 01d ####
png(
  file = file.path(figDir, 'FigShort_01d.png'),
  width  = gPars$reso,
  height = gPars$reso
)
sam = sample(M,M)
ErrViewLib::plotEvsPU(
  fHetero[sam],
  Z[sam],
  type = 'horiz',
  logX = FALSE,
  runQuant = TRUE,
  runMS = TRUE,
  scalePoints = scalePoints,
  xlim = c(-0.05,0.8),
  ylim = c(-5, 5),
  xlab = "Fraction of heteroatoms",
  label = 3,
  showLegend = FALSE,
  gPars = gPars
)
dev.off()

# ## Fig 02a ####

slide = FALSE
nBin  = 100

## Fig 02a ####
ylim = c(-0.5,0.5)
png(
  file = file.path(figDir, 'FigShort_02a.png'),
  width  = gPars$reso,
  height = gPars$reso
)
resZMu = ErrViewLib::plotLZM(
  uE,
  Z,
  ylim = ylim,
  xlab = 'Uncertainty, uE [eV]',
  logX = TRUE,
  nBin = nBin,
  slide = slide,
  equiPop = TRUE,
  label = 1,
  gPars = gPars
)
# title(main = paste0('fv = ',round(res$fVal,2),
#                     ', CI = [',round(res$lofVal,2),', ',
#                     round(res$upfVal,2),']'),
#       line = 1, font.main = 1, cex.main = 1.1)
dev.off()


## Fig 02b ####
png(
  file = file.path(figDir, 'FigShort_02b.png'),
  width  = gPars$reso,
  height = gPars$reso
)
resZMmass = ErrViewLib::plotLZM(
  masses,
  Z,
  xlim = c(50, 150),
  ylim = ylim,
  xlab = 'Molecular Mass [Da]',
  logX = FALSE,
  nBin = nBin,
  equiPop = TRUE,
  slide = slide,
  label = 2,
  gPars = gPars
)
# title(main = paste0('fv = ',round(res$fVal,2),
#                     ', CI = [',round(res$lofVal,2),', ',
#                     round(res$upfVal,2),']'),
#       line = 1, font.main = 1, cex.main = 1.1)
dev.off()

## Fig 02c ####
png(
  file = file.path(figDir, 'FigShort_02c.png'),
  width  = gPars$reso,
  height = gPars$reso
)
resZMfHet = ErrViewLib::plotLZM(
  fHetero,
  Z,
  xlim = c(-0.05, 0.8),
  ylim = ylim,
  xlab = 'Fraction of heteroatoms',
  logX = FALSE,
  nBin = nBin,
  equiPop = TRUE,
  slide = slide,
  label = 3,
  gPars = gPars
)
# title(main = paste0('fv = ',round(res$fVal,2),
#                     ', CI = [',round(res$lofVal,2),', ',
#                     round(res$upfVal,2),']'),
#       line = 1, font.main = 1, cex.main = 1.1)
dev.off()

## Fig 02d ####
ylim = c(0.0,2.5)
png(
  file = file.path(figDir, 'FigShort_02d.png'),
  width  = gPars$reso,
  height = gPars$reso
)
resZMSu = ErrViewLib::plotLZMS(
  uE,
  Z,
  ylim = ylim,
  xlab = 'Uncertainty, uE [eV]',
  logX = TRUE,
  nBin = nBin,
  method = 'bootstrap',
  slide = slide,
  equiPop = TRUE,
  score = TRUE,
  label = 4,
  gPars = gPars
)
# title(main = paste0('fv = ',round(res1$fVal,2),
#                     ', CI = [',round(res1$lofVal,2),', ',
#                     round(res1$upfVal,2),']'),
#       line = 1, font.main = 1, cex.main = 1.1)
dev.off()

## Fig 02e ####
png(
  file = file.path(figDir, 'FigShort_02e.png'),
  width  = gPars$reso,
  height = gPars$reso
)
resZMSmass = ErrViewLib::plotLZMS(
  masses,
  Z,
  xlim = c(50, 150),
  ylim = ylim,
  xlab = 'Molecular Mass [Da]',
  logX = FALSE,
  nBin = nBin,
  method = 'bootstrap',
  equiPop = TRUE,
  slide = slide,
  logBin = TRUE,
  score = TRUE,
  label = 5,
  gPars = gPars
)
# title(main = paste0('fv = ',round(res2$fVal,2),
#                     ', CI = [',round(res2$lofVal,2),', ',
#                     round(res2$upfVal,2),']'),
#       line = 1, font.main = 1, cex.main = 1.1)
dev.off()

## Fig 02f ####
png(
  file = file.path(figDir, 'FigShort_02f.png'),
  width  = gPars$reso,
  height = gPars$reso
)
resZMSfHet = ErrViewLib::plotLZMS(
  fHetero,
  Z,
  xlim = c(-0.05, 0.8),
  ylim = ylim,
  xlab = 'Fraction of heteroatoms',
  logX = FALSE,
  nBin = nBin,
  method = 'bootstrap',
  equiPop = TRUE,
  slide = slide,
  logBin = TRUE,
  score = TRUE,
  label = 6,
  gPars = gPars
)
# title(main = paste0('fv = ',round(res3$fVal,2),
#                     ', CI = [',round(res3$lofVal,2),', ',
#                     round(res3$upfVal,2),']'),
#       line = 1, font.main = 1, cex.main = 1.1)
dev.off()

# fvRes = rbind(
#   data.frame(
#     cv = 'uE',
#     fvZM = round(resZMu$fVal,2),
#     ciZM = paste0('[',round(resZMu$lofVal,2),', ',
#                       round(resZMu$upfVal,2),']'),
#     fvZMS = round(resZMSu$fVal,2),
#     ciZMS = paste0('[',round(resZMSu$lofVal,2),', ',
#                   round(resZMSu$upfVal,2),']')
#   ),
#   data.frame(
#     cv = 'X1',
#     fvZM = round(resZMmass$fVal,2),
#     ciZM = paste0('[',round(resZMmass$lofVal,2),', ',
#                   round(resZMmass$upfVal,2),']'),
#     fvZMS = round(resZMSmass$fVal,2),
#     ciZMS = paste0('[',round(resZMSmass$lofVal,2),', ',
#                    round(resZMSmass$upfVal,2),']')
#   ),
#   data.frame(
#     cv = 'X2',
#     fvZM = round(resZMfHet$fVal,2),
#     ciZM = paste0('[',round(resZMfHet$lofVal,2),', ',
#                   round(resZMfHet$upfVal,2),']'),
#     fvZMS = round(resZMSfHet$fVal,2),
#     ciZMS = paste0('[',round(resZMSfHet$lofVal,2),', ',
#                    round(resZMSfHet$upfVal,2),']')
#   )
# )

## Fig 02g ####
png(
  file = file.path(figDir, 'FigShort_02g.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotACF(res1$pc, label = 7, gPars = gPars)
dev.off()

## Fig 02h ####
png(
  file = file.path(figDir, 'FigShort_02h.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotACF(res2$pc, label = 8, gPars = gPars)
dev.off()

## Fig 02i ####
png(
  file = file.path(figDir, 'FigShort_02i.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotACF(res3$pc, label = 9, gPars = gPars)
dev.off()


## Fig 04a ####
ylim = c(-1, 1)
popMin = 100
png(
  file = file.path(figDir, 'FigShort_04a.png'),
  width  = gPars$reso,
  height = gPars$reso
)
resZMustr = ErrViewLib::plotStratZM(
  uE, Z,
  popMin = popMin,
  xlab = 'Uncertainty [eV]',
  ylim = ylim,
  label = 1,
  gPars = gPars)
# title(main = paste0('fv = ',round(res$fVal,2),
#                     ', CI = [',round(res$lofVal,2),', ',
#                     round(res$upfVal,2),']'),
#       line = 1, font.main = 1, cex.main = 1.1)
# resZMustr = res
dev.off()

## Fig 04b ####
png(
  file = file.path(figDir, 'FigShort_04b.png'),
  width  = gPars$reso,
  height = gPars$reso
)
resZMmassstr = ErrViewLib::plotStratZM(
  masses, Z,
  popMin = popMin,
  xlab = 'Molecular Mass [Da]',
  ylim = ylim,
  label = 2,
  gPars = gPars)
# title(main = paste0('fv = ',round(res$fVal,2),
#                     ', CI = [',round(res$lofVal,2),', ',
#                     round(res$upfVal,2),']'),
#       line = 1, font.main = 1, cex.main = 1.1)
# resmassZM = res
dev.off()

## Fig 04c ####
png(
  file = file.path(figDir, 'FigShort_04c.png'),
  width  = gPars$reso,
  height = gPars$reso
)
resZMfHetstr = ErrViewLib::plotStratZM(
  fHetero, Z,
  popMin = popMin,
  xlab = 'Fraction of heteroatoms',
  ylim = ylim,
  label = 3,
  gPars = gPars)
# title(main = paste0('fv = ',round(res$fVal,2),
#                     ', CI = [',round(res$lofVal,2),', ',
#                     round(res$upfVal,2),']'),
#       line = 1, font.main = 1, cex.main = 1.1)
# resfHetZM = res
dev.off()

## Fig 04d ####
ylim = c(0, 2)
png(
  file = file.path(figDir, 'FigShort_04d.png'),
  width  = gPars$reso,
  height = gPars$reso
)
resZMSustr = ErrViewLib::plotStratZMS(
  uE, Z,
  popMin = popMin,
  method = 'bootstrap',
  xlab = 'Uncertainty [eV]',
  ylim = ylim,
  label = 4,
  gPars = gPars)
# title(main = paste0('fv = ',round(res$fVal,2),
#                     ', CI = [',round(res$lofVal,2),', ',
#                     round(res$upfVal,2),']'),
#       line = 1, font.main = 1, cex.main = 1.1)
# resZMS = res
dev.off()

## Fig 04e ####
png(
  file = file.path(figDir, 'FigShort_04e.png'),
  width  = gPars$reso,
  height = gPars$reso
)
resZMSmassstr = ErrViewLib::plotStratZMS(
  masses, Z,
  popMin = popMin,
  method = 'bootstrap',
  xlab = 'Molecular Mass [Da]',
  ylim = ylim,
  label = 5,
  gPars = gPars)
# title(main = paste0('fv = ',round(res$fVal,2),
#                     ', CI = [',round(res$lofVal,2),', ',
#                     round(res$upfVal,2),']'),
#       line = 1, font.main = 1, cex.main = 1.1)
# resmassZMS = res
dev.off()

## Fig 04f ####
png(
  file = file.path(figDir, 'FigShort_04f.png'),
  width  = gPars$reso,
  height = gPars$reso
)
resZMSfHetstr = ErrViewLib::plotStratZMS(
  fHetero, Z,
  popMin = popMin,
  method = 'bootstrap',
  xlab = 'Fraction of heteroatoms',
  ylim = ylim,
  label = 6,
  gPars = gPars)
# title(main = paste0('fv = ',round(res$fVal,2),
#                     ', CI = [',round(res$lofVal,2),', ',
#                     round(res$upfVal,2),']'),
#       line = 1, font.main = 1, cex.main = 1.1)
# resfHetZMS = res
dev.off()

# fvResstr = rbind(
#   data.frame(
#     cv = 'uE',
#     fvZM = round(resZMustr$fVal,2),
#     ciZM = paste0('[',round(resZMustr$lofVal,2),', ',
#                   round(resZMustr$upfVal,2),']'),
#     fvZMS = round(resZMSustr$fVal,2),
#     ciZMS = paste0('[',round(resZMSustr$lofVal,2),', ',
#                    round(resZMSustr$upfVal,2),']')
#   ),
#   data.frame(
#     cv = 'X1',
#     fvZM = round(resZMmassstr$fVal,2),
#     ciZM = paste0('[',round(resZMmassstr$lofVal,2),', ',
#                   round(resZMmassstr$upfVal,2),']'),
#     fvZMS = round(resZMSmassstr$fVal,2),
#     ciZMS = paste0('[',round(resZMSmassstr$lofVal,2),', ',
#                    round(resZMSmassstr$upfVal,2),']')
#   ),
#   data.frame(
#     cv = 'X2',
#     fvZM = round(resZMfHetstr$fVal,2),
#     ciZM = paste0('[',round(resZMfHetstr$lofVal,2),', ',
#                   round(resZMfHetstr$upfVal,2),']'),
#     fvZMS = round(resZMSfHetstr$fVal,2),
#     ciZMS = paste0('[',round(resZMSfHetstr$lofVal,2),', ',
#                    round(resZMSfHetstr$upfVal,2),']')
#   )
# )

## Fig03 ####

nMC = 1000
nBin = 100

# resZM0 = ErrViewLib::plotLZM(
#   uE, Z,
#   nBin = nBin, plot = FALSE)
# resZMS0 = ErrViewLib::plotLZMS(
#   uE, Z,
#   nBin = nBin, plot = FALSE, score = TRUE)
# resmassZM0 = ErrViewLib::plotLZM(
#   masses, Z,
#   nBin = nBin, plot = FALSE)
# resmassZMS0 = ErrViewLib::plotLZMS(
#   masses, Z,
#   nBin = nBin, plot = FALSE, score = TRUE)
# resfHetZM0 = ErrViewLib::plotLZM(
#   fHetero, Z,
#   nBin = nBin, plot = FALSE)
# resfHetZMS0 = ErrViewLib::plotLZMS(
#   fHetero, Z,
#   nBin = nBin, plot = FALSE, score = TRUE)

fRes = paste0(figDir,'/res_fig03.Rda')
if(!file.exists(fRes)) {

  fVal = resZMu$isVal
  tabZM = matrix(NA,nrow = nMC, ncol = length(fVal))
  for(i in 1:nMC) {
    cat(i,'/ ')
    res = ErrViewLib::plotLZM(
      uE, Z, aux = sample(M,M),
      nBin = nBin, plot = FALSE)
    tabZM[i,] = res$isVal
  }
  tabZMS = matrix(NA,nrow = nMC, ncol = nBin)
  for(i in 1:nMC) {
    cat(i,'/ ')
    res = ErrViewLib::plotLZMS(
      uE, Z, aux = sample(M,M),
      nBin = nBin, plot = FALSE, score = TRUE)
    tabZMS[i,] = res$isVal
  }
  fVal = resZMmass$isVal
  tabmassZM = matrix(NA,nrow = nMC, ncol = length(fVal))
  for(i in 1:nMC) {
    cat(i,'/ ')
    res = ErrViewLib::plotLZM(
      masses, Z, aux = sample(M,M),
      nBin = nBin, plot = FALSE)
    tabmassZM[i,] = res$isVal
  }
  tabmassZMS = matrix(NA,nrow = nMC, ncol = nBin)
  for(i in 1:nMC) {
    cat(i,'/ ')
    res = ErrViewLib::plotLZMS(
      masses, Z, aux = sample(M,M),
      nBin = nBin, plot = FALSE, score = TRUE)
    tabmassZMS[i,] = res$isVal
  }
  fVal = resZMfHet$isVal
  tabfHetZM = matrix(NA,nrow = nMC, ncol = length(fVal))
  for(i in 1:nMC) {
    cat(i,'/ ')
    res = ErrViewLib::plotLZM(
      fHetero, Z, aux = sample(M,M),
      nBin = nBin, plot = FALSE)
    tabfHetZM[i,] = res$isVal
  }
  tabfHetZMS = matrix(NA,nrow = nMC, ncol = nBin)
  for(i in 1:nMC) {
    cat(i,'/ ')
    res = ErrViewLib::plotLZMS(
      fHetero, Z, aux = sample(M,M),
      nBin = nBin, plot = FALSE, score = TRUE)
    tabfHetZMS[i,] = res$isVal
  }

  save(tabZM, tabZMS, tabmassZM, tabmassZMS, tabfHetZM, tabfHetZMS,
       file = fRes)
} else {
  load(file = fRes)

}

resCI = list()

mZM     = rowMeans(tabZM)
resMCZM        = list()
resMCZM$fVal   = mean(mZM)
resMCZM$lofVal = quantile(mZM,0.025)
resMCZM$upfVal = quantile(mZM,0.975)
sZM     = rbinom(mZM,nBin,mZM) / nBin
resMCBIZM        = list()
resMCBIZM$fVal   = mean(sZM)
resMCBIZM$lofVal = quantile(sZM,0.025)
resMCBIZM$upfVal = quantile(sZM,0.975)
mZM     = rowMeans(tabZMS)
resMCZMS        = list()
resMCZMS$fVal   = mean(mZM)
resMCZMS$lofVal = quantile(mZM,0.025)
resMCZMS$upfVal = quantile(mZM,0.975)
sZM     = rbinom(mZM,nBin,mZM) / nBin
resMCBIZMS        = list()
resMCBIZMS$fVal   = mean(sZM)
resMCBIZMS$lofVal = quantile(sZM,0.025)
resMCBIZMS$upfVal = quantile(sZM,0.975)
resCI[['uE']]$title    = "Uncertainty"
resCI[['uE']]$mZM      = resMCZM
resCI[['uE']]$sZM      = resMCBIZM
resCI[['uE']]$mZMS     = resMCZMS
resCI[['uE']]$sZMS     = resMCBIZMS
resCI[['uE']]$ZM0      = resZMu
resCI[['uE']]$ZMS0     = resZMSu
resCI[['uE']]$stratZM  = resZMustr
resCI[['uE']]$stratZMS = resZMSustr

mZM     = rowMeans(tabmassZM)
resMCZM        = list()
resMCZM$fVal   = mean(mZM)
resMCZM$lofVal = quantile(mZM,0.025)
resMCZM$upfVal = quantile(mZM,0.975)
sZM     = rbinom(mZM,nBin,mZM) / nBin
resMCBIZM        = list()
resMCBIZM$fVal   = mean(sZM)
resMCBIZM$lofVal = quantile(sZM,0.025)
resMCBIZM$upfVal = quantile(sZM,0.975)
mZM     = rowMeans(tabmassZMS)
resMCZMS        = list()
resMCZMS$fVal   = mean(mZM)
resMCZMS$lofVal = quantile(mZM,0.025)
resMCZMS$upfVal = quantile(mZM,0.975)
sZM     = rbinom(mZM,nBin,mZM) / nBin
resMCBIZMS        = list()
resMCBIZMS$fVal   = mean(sZM)
resMCBIZMS$lofVal = quantile(sZM,0.025)
resMCBIZMS$upfVal = quantile(sZM,0.975)
resCI[['mass']]$title = "Molecular mass"
resCI[['mass']]$mZM      = resMCZM
resCI[['mass']]$sZM      = resMCBIZM
resCI[['mass']]$mZMS     = resMCZMS
resCI[['mass']]$sZMS     = resMCBIZMS
resCI[['mass']]$ZM0      = resZMmass
resCI[['mass']]$ZMS0     = resZMSmass
resCI[['mass']]$stratZM  = resZMmass
resCI[['mass']]$stratZMS = resZMSmassstr

mZM     = rowMeans(tabfHetZM)
resMCZM        = list()
resMCZM$fVal   = mean(mZM)
resMCZM$lofVal = quantile(mZM,0.025)
resMCZM$upfVal = quantile(mZM,0.975)
sZM     = rbinom(mZM,nBin,mZM) / nBin
resMCBIZM        = list()
resMCBIZM$fVal   = mean(sZM)
resMCBIZM$lofVal = quantile(sZM,0.025)
resMCBIZM$upfVal = quantile(sZM,0.975)
mZM     = rowMeans(tabfHetZMS)
resMCZMS        = list()
resMCZMS$fVal   = mean(mZM)
resMCZMS$lofVal = quantile(mZM,0.025)
resMCZMS$upfVal = quantile(mZM,0.975)
sZM     = rbinom(mZM,nBin,mZM) / nBin
resMCBIZMS        = list()
resMCBIZMS$fVal   = mean(sZM)
resMCBIZMS$lofVal = quantile(sZM,0.025)
resMCBIZMS$upfVal = quantile(sZM,0.975)
resCI[['fHet']]$title = "Fraction of heteroatoms"
resCI[['fHet']]$mZM      = resMCZM
resCI[['fHet']]$sZM      = resMCBIZM
resCI[['fHet']]$mZMS     = resMCZMS
resCI[['fHet']]$sZMS     = resMCBIZMS
resCI[['fHet']]$ZM0      = resZMfHet
resCI[['fHet']]$ZMS0     = resZMSfHet
resCI[['fHet']]$stratZM  = resZMfHetstr
resCI[['fHet']]$stratZMS = resZMSfHetstr

png(
  file = file.path(figDir, 'FigShort_03.png'),
  width  = 2.4*gPars$reso,
  height = gPars$reso
)
par(
  mfrow=c(1,3),
  mar = gPars$mar,
  mgp = gPars$mgp,
  pty = 'm',
  tcl = gPars$tcl,
  cex = gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)

ylim = c(0.2,1)
cols = gPars$cols

del = 0.15
xm2 = 1:2 - 1.5*del
xm1 = 1:2 - 0.5*del
xp1 = 1:2 + 0.5*del
xp2 = 1:2 + 1.5*del
pch0 = 20
for(i in seq_along(names(resCI))) {
  res = resCI[[names(resCI)[i]]]
  # Nominal stats
  plot(
    xm2,c(res$ZM0$fVal, res$ZMS0$fVal),
    type = 'p', pch = pch0 + 1, col = cols[1],
    xlim = c(0.5,2.5),
    xaxt = 'n',
    xlab ='',
    yaxs = 'i',
    ylim = ylim,
    ylab = 'Fraction of validated bins',
    main = res$title)
  grid(nx=0, ny = 4)
  abline(h=0.95, lty =2)
  mtext('0.95', side = 2, at = 0.95,
        cex = 0.75*gPars$cex, las = 2, line=0.25)
  axis(1,at = 1:2, labels = c('< Z >','< Z^2 >'))
  segments(
    xm2,c(res$ZM0$lofVal, res$ZMS0$lofVal),
    xm2,c(res$ZM0$upfVal, res$ZMS0$upfVal), col = cols[1], lwd = 8)

  # MC stats
  points(
    xm1,c(res$mZM$fVal, res$mZMS$fVal),pch = pch0 + 2, col = cols[2])
  segments(
    xm1,c(res$mZM$lofVal, res$mZMS$lofVal),
    xm1,c(res$mZM$upfVal, res$mZMS$upfVal), col = cols[2], lwd = 8)

  # MC + BI stats
  points(
    xp1, c(res$sZM$fVal,res$sZMS$fVal), col = cols[3], pch = pch0 + 3)
  segments(
    xp1,c(res$sZM$lofVal, res$sZMS$lofVal),
    xp1,c(res$sZM$upfVal, res$sZMS$upfVal), col = cols[3], lwd = 8)

  # Strat stats
  points(
    xp2,c(res$stratZM$fVal, res$stratZMS$fVal),pch = pch0 + 4, col = cols[4])
  segments(
    xp2,c(res$stratZM$lofVal, res$stratZMS$lofVal),
    xp2,c(res$stratZM$upfVal, res$stratZMS$upfVal), col = cols[4], lwd = 8)

  if(i==1)
    legend(
      'bottom', bty = 'n',
      legend = c(
        'Nominal',
        'Random order',
        'Random + Binomial',
        'Stratified'),
      lty = 1, lwd = 8,
      pch = pch0 + 1:4,
      col = cols[1:4]
    )
  box()

  label = i
  mtext(
    text = paste0('(', letters[label], ')'),
    side = 3,
    adj = 1,
    cex = gPars$cex,
    line = 0.25
  )
}

dev.off()



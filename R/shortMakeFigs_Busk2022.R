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

## Fig 02 ####
nBin = 100
png(
  file = file.path(figDir, 'Fig_02a.png'),
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

png(
  file = file.path(figDir, 'Fig_02b.png'),
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

png(
  file = file.path(figDir, 'Fig_02c.png'),
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


# ## Fig 03 ####

slide = FALSE
nBin  = 100

ylim = c(-0.5,0.5)
png(
  file = file.path(figDir, 'Fig_03a.png'),
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
dev.off()

png(
  file = file.path(figDir, 'Fig_03b.png'),
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
dev.off()

png(
  file = file.path(figDir, 'Fig_03c.png'),
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
dev.off()

ylim = c(0.0,2.5)
png(
  file = file.path(figDir, 'Fig_03d.png'),
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
dev.off()

png(
  file = file.path(figDir, 'Fig_03e.png'),
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
dev.off()

png(
  file = file.path(figDir, 'Fig_03f.png'),
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
dev.off()

png(
  file = file.path(figDir, 'Fig_03g.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotACF(resZMSu$pc, label = 7, gPars = gPars)
dev.off()

png(
  file = file.path(figDir, 'Fig_03h.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotACF(resZMSmass$pc, label = 8, gPars = gPars)
dev.off()

png(
  file = file.path(figDir, 'Fig_03i.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotACF(resZMSfHet$pc, label = 9, gPars = gPars)
dev.off()


## Fig 05 ####
ylim = c(-1, 1)
popMin = 100
png(
  file = file.path(figDir, 'Fig_05a.png'),
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
dev.off()

png(
  file = file.path(figDir, 'Fig_05b.png'),
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
dev.off()

png(
  file = file.path(figDir, 'Fig_05c.png'),
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
dev.off()

ylim = c(0, 2)
png(
  file = file.path(figDir, 'Fig_05d.png'),
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
dev.off()

png(
  file = file.path(figDir, 'Fig_05e.png'),
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
dev.off()

png(
  file = file.path(figDir, 'Fig_05f.png'),
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
dev.off()


## Fig04 ####

nMC = 1000
nBin = 100

fRes = paste0(figDir,'/res_fig04.Rda')
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
  file = file.path(figDir, 'Fig_04.png'),
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



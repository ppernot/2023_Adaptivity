
png(
  file = file.path(figDir, 'FigShort_05.png'),
  width  = gPars$reso,
  height = gPars$reso
)
# Z1 = rt(Z,df=5)/sqrt(5/3)
res = plotAGVZMS(
  Z,
  nMC = 1000,
  verbose = TRUE,
  nBoot = 500,
  method = 'bootstrap',
  col = 6,
  control = TRUE
)
dev.off()












# Reference AGC ####
fAGC <- function(
  mSeq,
  D,
  nMC = 100,
  nGran = 100,
  replace = FALSE,
  fscore = 'zms',
  zms_target = 1
) {

  Mn = NROW(D)

  agvm = agvs = c()
  for(i in seq_along(mSeq)) {
    M1 = mSeq[i]
    wv = c()
    for(j in 1:nMC) {
      worst = 0
      for(k in 1:nGran) {
        s = sample(Mn, M1, replace = replace)
        if(fscore == 'rce'){
          score = abs(ErrViewLib::rce(D,s))
        } else {
          z = D[s,2]/D[s,1]
          score = sqrt(abs(mean(z^2) - zms_target))
        }
        worst = max(worst, score)
      }
      wv[j] = worst
    }
    agvm[i] = mean(wv)
    agvs[i] = sd(wv)
  }
  return(list(
    agvm = agvm,
    agvs = agvs
  ))
}

mSeq = c(50, 100, 250, 500, 1000, 2000, 5000, M)
Mn   = M
x    = signif(mSeq/Mn,2)

uEn = rchisq(Mn, df = 5)

# Normal
En   = uEn * rnorm(Mn)
D    = cbind(uEn, En)
set.seed(123)
resn = fAGC(mSeq, D)
set.seed(123)
resn1 = fAGC(mSeq, D, fscore = 'rce')

# # Student
# nu    = 9
# En   = uEn * rt(Mn, df = nu) / sqrt(nu / (nu - 2))
# D    = cbind(uEn, En)
# rest = fAGC(mSeq, D)
# rest1 = fAGC(mSeq, D, fscore = 'rce')

# QM9
D = cbind(uE, E)
set.seed(123)
resq  = fAGC(mSeq, D)
set.seed(123)
resq1 = fAGC(mSeq, D, fscore = 'rce')

# Plot
par(mfrow=c(1,2))

## ZMS
agvm = resn$agvm
agvs = resn$agvs
x    = signif(mSeq / Mn, 2)
plot(
  x,
  agvm,
  type = 'l',
  col = gPars$cols[2],
  log = 'x',
  ylim = c(0, 2),
  yaxs = 'i',
  main = 'ZMS'
)
grid()
polygon(c(x, rev(x)),
        c(agvm - 2 * agvs, rev(agvm + 2 * agvs)),
        border = NA,
        col = gPars$cols_tr[2])

# agvm = rest$agvm
# agvs = rest$agvs
# lines(x,agvm, col = 2)
# polygon(c(x,rev(x)),c(agvm-2*agvs,rev(agvm+2*agvs)),border = NA, col = gPars$cols_tr[2])

agvm = resq$agvm
agvs = resq$agvs
lines(x, agvm, col = gPars$cols[6])
polygon(c(x, rev(x)),
        c(agvm - 2 * agvs, rev(agvm + 2 * agvs)),
        border = NA,
        col = gPars$cols_tr[6])

## RCE
agvm = resn1$agvm
agvs = resn1$agvs
x    = signif(mSeq / Mn, 2)
plot(
  x,
  agvm,
  type = 'l',
  col = gPars$cols[2],
  log = 'x',
  ylim = c(0, 2),
  yaxs = 'i',
  main = "RCE"
)
grid()
polygon(c(x, rev(x)),
        c(agvm - 2 * agvs, rev(agvm + 2 * agvs)),
        border = NA,
        col = gPars$cols_tr[2])

# agvm = rest1$agvm
# agvs = rest1$agvs
# lines(x,agvm, col = 2)
# polygon(c(x,rev(x)),c(agvm-2*agvs,rev(agvm+2*agvs)),border = NA, col = gPars$cols_tr[2])

agvm = resq1$agvm
agvs = resq1$agvs
lines(x, agvm, col = gPars$cols[6])
polygon(c(x, rev(x)),
        c(agvm - 2 * agvs, rev(agvm + 2 * agvs)),
        border = NA,
        col = gPars$cols_tr[6])















# Dispersion of stats in AGV ####
uE1 = rchisq(Z,df = 5)
E1  = uE1 * rt(Z, df = 5) / sqrt(5/3)
Y1 = cbind(uE1,E1)
Z1 = E1 / uE1

res = ErrViewLib::plotAGVZMS(
  Z1,
  method = 'stud',
  control = TRUE
)

Y = cbind(uE,E)

popMin = 50; popMax = 1000
nMC = 10000
M = length(Z)
sizes  = M / 2^{1:20}
sizes  = round(sizes[sizes >= popMin & sizes <= popMax])

mz = matrix(NA,nrow = nMC, ncol = length(sizes))
mz1 = matrix(NA,nrow = nMC, ncol = length(sizes))
for (i in seq_along(sizes)) {
  for (j in 1:nMC) {
    sam = sample(M, sizes[i])
    mz[j,i]  = abs(ErrViewLib::rce(Y ,sam)) # mean(Z[sam]^2)
    mz1[j,i] = abs(ErrViewLib::rce(Y1,sam)) # mean(Z1[sam]^2)
  }
}
library(vioplot)
x = rev(1:length(sizes))
vioplot::vioplot(mz1,side = 'left', col = cols[6], plotCentre = 'line',
                 rectCol = cols[6], xlim = c(0,6), ylim = c(0,2.5), at = x)
vioplot::vioplot(mz,side = 'right', col = cols[2], plotCentre = 'line',
                 rectCol = cols[2], add = TRUE, at = x)
# abline(h=1, lty = 2)
grid()
vioplot::vioplot(mz1,side = 'left', col = cols[6], plotCentre = 'line',
                 rectCol = cols[6], ylim = c(0,2.5), add = TRUE, at = x)
vioplot::vioplot(mz,side = 'right', col = cols[2], plotCentre = 'line',
                 rectCol = cols[2], add = TRUE, at = x)
box()

# Coverage of statistics CIs ####

nMC   = 1000
nBoot = 1500
mSeq  = c(50,100,250,500)
nu = 5

# Check coverage of <Z> and <Z^2>
vec = c()
for(i in seq_along(mSeq)) {
  M1 = mSeq[i]
  cat(M1,'/')
  OK = 0
  for(j in 1:nMC) {
    # Z1 = rnorm(M1,0,1)
    # u1  = rchisq(M1, df = 4)
    E1  = rnorm(M1, 0, u1)
    E1  = u1 * rt(M1, df = nu) / sqrt(nu / (nu - 2))
    Z1 = E1/u1
    m   = mean(Z1)
    loM = m + sd(Z1)/sqrt(M1) * qt(0.025, df = M1-1)
    upM = m + sd(Z1)/sqrt(M1) * qt(0.975, df = M1-1)
    OK  = OK + as.numeric(loM * upM <= 0)
  }
  vec[i] = OK / nMC; print(vec[i])
}
vecMZ = vec



MSZw = function(X, w) {
  weighted.mean(X^2,w)
}
MSZ = function(X, indx = 1:length(X)) {
  mean(X[indx]^2)
}
nMC   = 1000
nBoot = 1500
v1 = c()
b = 0
nu = 5
for(i in seq_along(mSeq)) {
  M1 = mSeq[i]
  cat(M1,'/')
  OK = 0
  for(j in 1:nMC) {
    u1  = rchisq(M1, df = 4)
    # E1  = rnorm(M1, 0, u1)
    E1  = b + u1 * rt(M1, df = nu) / sqrt(nu / (nu - 2))
    u2  = sqrt(u1^2 + b^2)
    Z1  = E1 / u2

    bs  = boot::boot(Z1, MSZ, R = nBoot)
    bci = boot::boot.ci(bs, conf = 0.95, type = 'bca')
    lo  = bci$bca[1, 4]
    up  = bci$bca[1, 5]
    # bci = boot::abc.ci(Z1, MSZw, conf = 0.95)
    # lo  = bci[2]
    # up  = bci[3]

    OK  = OK + as.numeric((lo-1) * (up-1) <= 0)
  }
  v1[i] = OK / nMC; print(v1[i])
}
vecMZS1 = v1

v2 = c()
b = 1
for(i in seq_along(mSeq)) {
  M1 = mSeq[i]
  cat(M1,'/')
  OK = 0
  for(j in 1:nMC) {
    u1  = rchisq(M1, df = 4)
    E1  = rnorm(M1, b, u1)
    u2  = sqrt(u1^2 + b^2)
    Z1  = E1 / u2

    m = mean(Z1^2)
    s = sd(Z1^2) / sqrt(M1)
    lo = m + s * qt(0.025, df = M-1)
    up = m + s * qt(0.975, df = M-1)
    OK  = OK + as.numeric((lo-1) * (up-1) <= 0)
  }
  v2[i] = OK / nMC; print(v2[i])
}
vecMZS2 = v2

v3 = c()
b = 1
for(i in seq_along(mSeq)) {
  M1 = mSeq[i]
  cat(M1,'/')
  OK = 0
  nu = 5
  mZ = c()
  for(j in 1:nMC) {
    u1  = rchisq(M1, df = 4)
    E1  = b + u1 * rt(M1, df=nu) / sqrt(nu/(nu-2))
    u2  = sqrt(u1^2 + b^2)
    Z1  = E1 / u2

    m = mean(Z1^2); mZ[j] = m
    s = sd(Z1^2) / sqrt(M1)
    lo = m + s * qt(0.025, df = M-1)
    up = m + s * qt(0.975, df = M-1)
    OK  = OK + as.numeric((lo-1) * (up-1) <= 0)
  }
  v3[i] = OK / nMC; print(v3[i])
}
vecMZS3 = v3

matplot(
  mSeq, vecMZ, log = 'x',
  type = 'b', pch = 1, lty = 1, col = 1,
  xlab = "Sample size",
  ylab = "Coverage", ylim = c(0.8,1.0))
for(j in 1:length(vecMZ)) {
  ci      = DescTools::BinomCI(vecMZ[j]*nMC, nMC, method = "wilsoncc")
  lofVal  = ci[,2]
  upfVal  = ci[,3]
  segments(mSeq[j], lofVal, mSeq[j], upfVal, col = 1)
}

grid(equilogs = FALSE)
abline(h=0.95, lty = 2)
legend(
  'bottomright', bty = 'n',
  legend = c('< Z >', '< Z^2 >'),
  lty = 1, pch = 1, col = c(1,2,4)
)
box()

points(
  mSeq, vecMZS2,
  type = 'b', pch = 1, lty = 1, col = 2)
for(j in 1:length(vecMZS2)) {
  ci      = DescTools::BinomCI(vecMZS2[j]*nMC, nMC, method = "wilsoncc")
  lofVal  = ci[,2]
  upfVal  = ci[,3]
  segments(mSeq[j], lofVal, mSeq[j], upfVal, col = 2)
}

points(
  mSeq, vecMZS1,
  type = 'b', pch = 1, lty = 1, col = 4)
for(j in 1:length(vecMZS1)) {
  ci      = DescTools::BinomCI(vecMZS1[j]*nMC, nMC, method = "wilsoncc")
  lofVal  = ci[,2]
  upfVal  = ci[,3]
  segments(mSeq[j], lofVal, mSeq[j], upfVal, col = 4)
}

points(
  mSeq, vecMZS3,
  type = 'b', pch = 1, lty = 1, col = 3)
for(j in 1:length(vecMZS3)) {
  ci      = DescTools::BinomCI(vecMZS3[j]*nMC, nMC, method = "wilsoncc")
  lofVal  = ci[,2]
  upfVal  = ci[,3]
  segments(mSeq[j], lofVal, mSeq[j], upfVal, col = 3)
}

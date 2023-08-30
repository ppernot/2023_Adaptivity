nMC = 10000

Mn = 10000
mZ2 = 1

nSeq = c(50, 100, 250, 500, 1000, 2000)

# Normal sample
Zn = rnorm(Mn)
nBin = 100
res0 = ErrViewLib::plotLZMS(
  1:Mn, Zn,
  nBin = nBin,
  method = 'bootstrap'
)

eZMS0 = c()
for(j in seq_along(nSeq)) {
  N = nSeq[j]; print(N)
  vm = 0
  for(i in 1:nMC) {
    z = Zn[sample(Mn,N)]
    s = abs(mean(z^2) - mZ2)
    vm = max(vm,s)
  }
  eZMS0[j] = vm
}

# Non-normal sample
nu = 5
Zn = rt(Mn, df = nu) / sqrt(nu / (nu - 2))
nBin = 100
res1 = ErrViewLib::plotLZMS(
  1:Mn, Zn,
  nBin = nBin,
  method = 'bootstrap'
)

eZMS1 = c()
for(j in seq_along(nSeq)) {
  N = nSeq[j]; print(N)
  vm = 0
  for(i in 1:nMC) {
    z = Zn[sample(Mn,N)]
    s = abs(mean(z^2) - mZ2)
    vm = max(vm,s)
  }
  eZMS1[j] = vm
}

# QM9 sample
Mn = M
Zn = Z

eZMS2 = c()
for(j in seq_along(nSeq)) {
  N = nSeq[j]; print(N)
  vm = 0
  for(i in 1:nMC) {
    z = Zn[sample(Mn,N)]
    s = abs(mean(z^2) - mZ2)
    vm = max(vm,s)
  }
  eZMS2[j] = vm
}

# Plot
par(mfrow =c (1,3))
del  = 1.02
xlab = 'Group size'
ylab = 'Max <Z^2> abs. error'

par(mfrow=c(1,1))
matplot(nSeq, cbind(eZMS0,eZMS1,eZMS2), ylim = c(0,3),
        type = 'b', pch = 16, lty = 1, log = 'x',
        xlab = xlab, ylab = ylab)
grid()
box()



# ACF
nBin = 100

uEn = rchisq(Mn,df=4)
En  = uEn * rnorm(Mn)
Zn = En / uEn
res0 = ErrViewLib::plotLZMS(1:Mn, Zn, nBin = nBin,
                            method = 'stud', plot = FALSE)
rce0 = ErrViewLib::plotLRCE(1:Mn, uEn, En, nBin = nBin, plot = TRUE)

# nu = 5
En = E #uEn * rt(Mn, df = nu) / sqrt(nu / (nu - 2))
uEn = uE
Zn = En / uEn
res1 = ErrViewLib::plotLZMS(1:Mn, Zn, nBin = nBin,
                            method = 'stud', plot = FALSE)
rce1 = ErrViewLib::plotLRCE(1:Mn, uEn, En, nBin = nBin, plot = FALSE)


par(mfrow=c(2,4))
hist(res0$pc-1,xlim = c(-1,1), main = 'ZMS-1')
acf(res0$pc, main = 'ZMS')
hist(rce0$pc,xlim = c(-1,1), main = signif(max(abs(rce0$pc)),2))
acf(rce0$pc, main = 'RCE')
hist(res1$pc-1,xlim = c(-1,1), main = 'ZMS-1')
acf(res1$pc, main = 'ZMS')
hist(rce1$pc,xlim = c(-1,1), main = signif(max(abs(rce1$pc)),2))
acf(rce1$pc,main = 'RCE')


rce1 = ErrViewLib::plotLRCE(uEn, uEn, En, nBin = nBin, logX = TRUE, plot = TRUE)
rce1 = ErrViewLib::plotLRCE(masses, uEn, En, nBin = nBin, plot = TRUE)
rce1 = ErrViewLib::plotLRCE(fHetero, uEn, En, nBin = nBin, plot = TRUE)

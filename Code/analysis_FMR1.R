##########################################################################################
##########################################################################################
###  Analysis of the MCMC output under the hidden Markov model (HMM) 
###  for the FMR1 data at sites 1-22  
###
###  Audrey Qiuyan Fu
###  2010.10.24
###
###  Fu et al. In vivo properties of human DNA methyltransferases inferred from 
###  methylation patterns. 
###
###  See User Manual.
##########################################################################################
##########################################################################################

###############################################################
# read in FMR1 data and other info
###############################################################
dataDir = "../Data/"
file.in = paste (dataDir, "FMR1through22SEPrev.txt", sep="")

fmr1 = read.table (file.in, header = FALSE, sep=",")
dim (fmr1)
fmr1.data = as.matrix (fmr1[,-1])
n.site = ncol (fmr1.data)
n.proc = 3
n.start = n.proc*2

# nucleotide positions of CpG sites
loc = c(0, 13, 15, 25, 36, 52, 68, 70, 74, 78, 83, 85, 87, 98, 103, 111, 115, 128, 130, 136, 138, 140)
dist = diff (loc)

# generate prior density of hemi-preference ratio
x1 = runif (1e6)
x2 = runif (1e6)
x = log (x1/x2)
x.den = density (x)
plot (x.den$x, x.den$y, type="l")

###############################################################
# analyze MCMC output
###############################################################
mcmcDir = "../Data/HMM/"
file.pars.mcmc = paste (mcmcDir, "FMR1Human22_hmm_v6_1_pars_mcmc.out", sep="")
par.mcmc1 = as.matrix(read.table (file.pars.mcmc, header = FALSE))
file.pars.mcmc = paste (mcmcDir, "FMR1Human22_hmm_v6_1b_pars_mcmc.out", sep="")
par.mcmc2 = as.matrix(read.table (file.pars.mcmc, header = FALSE))
file.pars.mcmc = paste (mcmcDir, "FMR1Human22_hmm_v6_1c_pars_mcmc.out", sep="")
par.mcmc3 = as.matrix(read.table (file.pars.mcmc, header = FALSE))

par.mcmc = as.matrix (rbind (par.mcmc1, par.mcmc2, par.mcmc3))
mcmc.length = nrow (par.mcmc)

tau.mcmc = par.mcmc[,1:3]
rho.mcmc = par.mcmc[,4:6]
m.mcmc = par.mcmc[,(n.start+1):(n.start+n.site)]
rm.mcmc = par.mcmc[,n.start+n.site+1]
gm.mcmc = par.mcmc[,n.start+n.site+2]
c.mcmc = par.mcmc[,n.start+n.site+3]
md.mcmc = par.mcmc[,n.start+n.site+4]
mm.mcmc = par.mcmc[,n.start+n.site+5]

pi = tau.mcmc / (tau.mcmc + rho.mcmc * (1-tau.mcmc))
proc1 = 1/rho.mcmc
proc0 = 1/tau.mcmc
ratio = mm.mcmc / md.mcmc
lratio = log (ratio)

apply (pi, 2, quantile, probs=c(0.1,0.5,0.9))
#10% 0.979017 0.003530928 0.01172001
#50% 0.987352 0.020477935 0.04661654
#90% 0.993037 0.054837175 0.08312758
apply (tau.mcmc, 2, quantile, probs=c(0.1,0.5,0.9))
#10% 0.06864755 0.001501866 0.002493081
#50% 0.12071875 0.010276640 0.010909495
#90% 0.20450711 0.030724513 0.031640273
apply (rho.mcmc, 2, quantile, probs=c(0.1,0.5,0.9))
#10% 0.0006576621 0.1601562 0.07171082
#50% 0.0016748520 0.6466910 0.27632355
#90% 0.0045614420 0.9343106 0.79098728
apply (par.mcmc[,n.start+n.site+1:5], 2, quantile, probs=c(0.1,0.5,0.9))
#			rm		  gm         c           md         mm
#10% 0.8320747 0.07237764 0.01230639 0.00325899 0.9744363
#50% 0.8712380 0.10855620 0.01929070 0.01689118 0.9869354
#90% 0.8976006 0.17049204 0.02510150 0.05100420 0.9973514
apply (proc1, 2, quantile, probs=c(0.1,0.5,0.9))
#10%  219.2317 1.070308  1.264246
#50%  597.0677 1.546334  3.619019
#90% 1520.5413 6.244106 13.945165
apply (proc0, 2, quantile, probs=c(0.1,0.5,0.9))
#10%  4.889806  32.54733  31.60570
#50%  8.283719  97.30893  91.66498
#90% 14.567318 665.83832 401.13400
quantile (ratio, probs=c(0.1,0.5,0.9))
#19.46348  58.48406 303.28738 

# maintenance and de novo methylation event probabilities
m.maint = pi[,1]*mm.mcmc + pi[,3] - pi[,1]*mm.mcmc*pi[,3]
d.denovo = pi[,1]*md.mcmc + pi[,3] - pi[,1]*md.mcmc*pi[,3]

quantile (m.maint, c(0.1, 0.5, 0.9))
#0.9622987 0.9739020 0.9855840 
quantile (d.denovo, c(0.1, 0.5, 0.9))
#0.03497444 0.06944821 0.10388080 

# Generate histograms of key parameters
plotDir = "../Plots/"
# Fig 3 in paper
plotfilename = paste (plotDir, "Plot_FMR1Human22_hmm_v6_1_hist_rates.ps", sep="")
postscript (plotfilename)
cexlab.tmp = 2.5
cexaxis.tmp = 2
par (mfcol=c(3,2), mar=c(4.5,3.5,4.5,2.5), oma=c(0.5, 0.5, 0.5, 5))
par (xpd=TRUE)
hist (rho.mcmc[,1], breaks=20, xlim=c(0,0.02), freq = FALSE, col="red", main="", xlab=expression(rho[M]), ylab="", cex.axis = cexaxis.tmp, cex.lab = cexlab.tmp)
segments (0, 1, 0.02, 1, lwd=2)
text (-0.002, 510, "a", cex=2)
hist (rho.mcmc[,2], breaks=20, xlim=c(0,1), ylim=c(0,2), freq = FALSE, col="red", main="", xlab=expression(rho[RP]), ylab="", cex.axis = cexaxis.tmp, cex.lab = cexlab.tmp)
segments (0, 1, 1, 1, lwd=2)
text (-0.1, 2.8, "b", cex=2)
hist (rho.mcmc[,3], breaks=20, xlim=c(0,1), ylim=c(0,2), freq = FALSE, col="red", main="", xlab=expression(rho[RD]), ylab="", cex.axis = cexaxis.tmp, cex.lab = cexlab.tmp)
segments (0, 1, 1, 1, lwd=2)
text (-0.1, 2.8, "c", cex=2)
par (mar=c(4.5,5.5,4.5,0.5))
hist (mm.mcmc, breaks=20, xlim=c(0.9,1), freq=FALSE, col="blue", main="", xlab=expression(mu[M]), ylab="", cex.axis=cexaxis.tmp, cex.lab=cexlab.tmp)
text (0.89, 72, "d", cex=2)
hist (md.mcmc, breaks=20, xlim=c(0,0.15), freq=FALSE, col="blue", main="", xlab=expression(delta[M]), ylab="", cex.axis=cexaxis.tmp, cex.lab=cexlab.tmp)
text (-0.015, 44, "e", cex=2)
hist (lratio[lratio<10], breaks=20, freq = FALSE, col="blue", xlim=c(-10, 10), ylim=c(0,0.4), main="", xlab=expression(ln(mu[M]/delta[M])), ylab="", cex.axis = cexaxis.tmp, cex.lab = cexlab.tmp)
lines (x.den$x[x.den$x<10 & x.den$x>-10], x.den$y[x.den$x<10 & x.den$x>-10], lwd=2)
text (-12, 0.57, "f", cex=2)
dev.off ()

# Supplementary Fig 2
plotfilename = paste (plotDir, "Plot_FMR1Human22_hmm_v6_discrete_1_hist_c.ps", sep="")
postscript (plotfilename)
par (mfrow=c(1,1))
hist (c.mcmc, breaks=seq(0,0.06,0.01/3), xlim=c(0,0.06), freq=FALSE, main="", xlab=expression(c), ylab="", cex.axis=1.6, cex.lab=1.6)
dev.off ()



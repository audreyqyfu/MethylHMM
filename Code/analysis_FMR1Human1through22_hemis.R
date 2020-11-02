####################################################################
# R code for inference of hemimethylated CpG dyads in the FMR1 data 
# under the hidden Markov model (HMM) and under the indendpent events 
# model (IEM).
#################################################################### 

############################################
# compile R code
############################################
codeDir = "Tools/"
source (paste(codeDir,"computeLatentEmiss.R", sep=""))
source (paste(codeDir,"computeLatentMarginal.R", sep=""))
source (paste(codeDir,"computePostProbEventHemi.R", sep=""))
source (paste(codeDir,"findHemiCoord.R", sep=""))

############################################
# read in data and MCMC output
############################################
# read in data
dataDir = "../Data/"
file.in <- paste (dataDir, "FMR1through22SEPrev.txt", sep="")

fmr1 = read.table (file.in, header = FALSE, sep=",")
dim (fmr1)
fmr1.data = as.matrix (fmr1[,2:ncol(fmr1)])
n.site = ncol (fmr1.data)
n.proc = 3
n.start = n.proc*2

loc = c(0, 13, 15, 25, 36, 52, 68, 70, 74, 78, 83, 85, 87, 98, 103, 111, 115, 128, 130, 136, 138, 140)
dist = diff (loc)

############################################
# under the HMM
############################################
# read in mcmc output
mcmcDir = "../Data/HMM/"
file.pars.mcmc = paste (mcmcDir, "FMR1Human22_hmm_v6_1_pars_mcmc.out", sep="")
par.mcmc1 = as.matrix(read.table (file.pars.mcmc, header = FALSE))
file.pars.mcmc = paste (mcmcDir, "FMR1Human22_hmm_v6_1b_pars_mcmc.out", sep="")
par.mcmc2 = as.matrix(read.table (file.pars.mcmc, header = FALSE))
file.pars.mcmc = paste (mcmcDir, "FMR1Human22_hmm_v6_1c_pars_mcmc.out", sep="")
par.mcmc3 = as.matrix(read.table (file.pars.mcmc, header = FALSE))
par.mcmc = as.matrix (rbind (par.mcmc1, par.mcmc2, par.mcmc3))
tau.mcmc = par.mcmc[,1:3]
rho.mcmc = par.mcmc[,4:6]
m.mcmc = par.mcmc[,(n.start+1):(n.start+n.site)]
rm.mcmc = par.mcmc[,n.start+n.site+1]
gm.mcmc = par.mcmc[,n.start+n.site+2]
c.mcmc = par.mcmc[,n.start+n.site+3]
md.mcmc = par.mcmc[,n.start+n.site+4]
mm.mcmc = par.mcmc[,n.start+n.site+5]

# read in strand assignment probabilities
file.pr = paste (mcmcDir, "FMR1Human22_hmm_v6_1_strandtype_mcmc.out", sep="")
assign.pr1 = as.matrix (read.table (file.pr, header = FALSE))
file.pr = paste (mcmcDir, "FMR1Human22_hmm_v6_1b_strandtype_mcmc.out", sep="")
assign.pr2 = as.matrix (read.table (file.pr, header = FALSE))
file.pr = paste (mcmcDir, "FMR1Human22_hmm_v6_1c_strandtype_mcmc.out", sep="")
assign.pr3 = as.matrix (read.table (file.pr, header = FALSE))
assign.pr = as.matrix (rbind (assign.pr1, assign.pr2, assign.pr3))

######################################################
# compute Pr(event | hemi)
######################################################
# process 600 MCMC samples
# can take up to a minute
postprobeventHMM = computePostProbEventHemi (ds.data=fmr1.data, assign.pr=assign.pr, m=m.mcmc, c=c.mcmc, tau=tau.mcmc, rho=rho.mcmc, md=md.mcmc, mm=mm.mcmc, model="HMM")

postprobeventHMM.mean = apply (postprobeventHMM, c(2,3), mean)


############################################
# under the independent events model
############################################
# MCMC samples
mcmcDir = "../Data/IEM/"
file.pars.mcmc = paste (mcmcDir, "1through22_ext_est_site_estc_1.out", sep="")
par.mcmc1 = as.matrix(read.table (file.pars.mcmc, header = FALSE))
file.pars.mcmc = paste (mcmcDir, "1through22_ext_est_site_estc_2.out", sep="")
par.mcmc2 = as.matrix(read.table (file.pars.mcmc, header = FALSE))
file.pars.mcmc = paste (mcmcDir, "1through22_ext_est_site_estc_3.out", sep="")
par.mcmc3 = as.matrix(read.table (file.pars.mcmc, header = FALSE))
par.mcmc = as.matrix (rbind (par.mcmc1, par.mcmc2, par.mcmc3))

m.mcmc <- par.mcmc[, 1:n.site]
mu.mcmc <- par.mcmc[, (n.site+1):(2*n.site)]
dp.mcmc <- par.mcmc[, (2*n.site+1):(3*n.site)]
dd.mcmc <- par.mcmc[, (3*n.site+1):(4*n.site)]
c.mcmc <- par.mcmc[, (4*n.site+1):(5*n.site)]

# read in strand assignment probabilities
file.pr = paste (mcmcDir, "1through22_ext_est_site_estc_strandtype_1.out", sep="")
assign.pr1 = as.matrix (read.table (file.pr, header = FALSE))
file.pr = paste (mcmcDir, "1through22_ext_est_site_estc_strandtype_2.out", sep="")
assign.pr2 = as.matrix (read.table (file.pr, header = FALSE))
file.pr = paste (mcmcDir, "1through22_ext_est_site_estc_strandtype_3.out", sep="")
assign.pr3 = as.matrix (read.table (file.pr, header = FALSE))
assign.pr = as.matrix (rbind (assign.pr1, assign.pr2, assign.pr3))

# process 1728 MCMC samples
# can take a minute or two
postprobeventIEM = computePostProbEventHemi (ds.data=fmr1.data, assign.pr=assign.pr, m=m.mcmc, c=c.mcmc, mu=mu.mcmc, dp=dp.mcmc, dd=dd.mcmc, model="IEM")

postprobeventIEM.mean = apply (postprobeventIEM, c(2,3), mean)

###########################################################
# scatterplots of posterior probabilities under two models
###########################################################
# scatter plot for posterior probabilities under two models
# coloring hemis on 4 most informative patterns
# find indices of these hemis
h.coord = findHemiCoord (fmr1.data)
patts4 = c(60, 82, 159, 165)
hemis.patts4.index = which (h.coord[,1]==60) 
for (i in 2:4)
	hemis.patts4.index = c(hemis.patts4.index, which (h.coord[,1]==patts4[i]))

# also coloring hemis at outlying CpG sites identified by the independent events model
# find indices of these hemis
sites4 = c(10, 14, 15, 16)
hemis.sites4.index = which (h.coord[,2]==sites4[1]) 
for (i in 2:4)
	hemis.sites4.index = c(hemis.sites4.index, which (h.coord[,2]==sites4[i]))

# graphical parameters
title.vec = c ("Failure of maintenance", "De novo on parent", "De novo on daughter", "Measurement error")
lower.vec = c (0, 0, 0, 0)
upper.vec = c (1, 1, 1, 1)
lwd.tmp = 0.4
cex.tmp = 0.6

# uncomment following commented lines to save the plot
plotDir = "../Plots/"
postscript (paste (plotDir,"Plot_FMR1Human22_HMM_IEM_postprobhemis.ps", sep=""))
par (mfrow=c(2,2), mar=c(4.5,4.5,4.5,0.5), oma=c(0.5, 0.5, 0.5, 5))
for (j in 1:4)
{
	plot (postprobeventIEM.mean[,j], postprobeventHMM.mean[,j], xlim=c(lower.vec[j], upper.vec[j]), ylim=c(lower.vec[j], upper.vec[j]), pch=16, lwd=lwd.tmp, cex=cex.tmp, xlab="Indendpent events model", ylab="Hidden Markov model", main=title.vec[j])
	abline (a=0,b=1,lty=2, xlim=c(lower.vec[j], upper.vec[j]), ylim=c(lower.vec[j], upper.vec[j]))
	points (postprobeventIEM.mean[hemis.sites4.index,j], postprobeventHMM.mean[hemis.sites4.index,j], pch=16, cex=cex.tmp, lwd=lwd.tmp, col="blue")
	points (postprobeventIEM.mean[hemis.patts4.index,j], postprobeventHMM.mean[hemis.patts4.index,j], pch=16, cex=cex.tmp, lwd=lwd.tmp, col="red")
}
dev.off ()

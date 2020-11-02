##########################################################################################
# main_test_invitro.R
#
# Audrey Qiuyan Fu 
# 2010.10.24
#
# Fu AQ, Genereux DP, Stoeger R, Laird CD and Stephens M.
# In vivo properties of DNA methyltransferases inferred from methylation patterns. 
#
# THE MCMC SETUP IN THIS FILE RUNS FOR ABOUT 30 SECONDS.
#
# Command lines to run the MCMC program to estimate several properties of DNA methyltransferases
# using double-stranded DNA methylation patterns.  

# The in vitro data from Goayl et al. (2006) Fig. 2B are used here.
#
# See User Manual for instructions.
###########################################################################################

# to compile C code
# open terminal window, go to where this file (main_test.R) is located
# enter each of following two lines in terminal window
#
# R CMD SHLIB ./Source/StrandAssignmentProbs.c -o ./Source/StrandAssignmentProbs.so -lm
# R CMD SHLIB ./Source/loglikHMM.c -o ./Source/loglikHMM.so -lm

#################################################
### Basic setups (MODIFY FOR YOUR OWN RUN)
#################################################
## input file
codeDir = "Source/"
dataDir = "../Data/"
file.in <- paste (dataDir, "Goyal2B.dat", sep="")
 
data.tmp = as.matrix (read.delim (file.in, header=FALSE, sep=","))
N = nrow (data.tmp)
S = ncol (data.tmp)
data = matrix (1, nrow=2*N, ncol=S)
data[2*(1:N),] = data.tmp
dim (data)

## output files
file.pars.mcmc = paste (dataDir, "test_Goayl2B_pars_mcmc.out", sep="")
file.strandtype = paste (dataDir, "test_Goyal2B_strandtype_mcmc.out", sep="")

## MCMC specs
n.iter = 100
burn.in = 0.2
step.size = 5
seed.value = 22

## Nucleotide positions of CpG sites
# substrate 1 in Goyal et al 2006
loc = c(0, 4, 25, 27, 34, 39, 48, 54, 60, 62, 98, 103, 114, 121, 127, 130, 133, 160, 162, 171, 188, 194, 214, 219, 222, 234, 240, 250, 263, 267, 278, 281, 289, 297, 319, 327, 339, 355, 361, 368, 370, 402, 421, 434, 451, 465, 491, 501, 516, 534, 541, 550, 554, 581)
# substrate for Vilkaitis et al 2005 Fig. 7A
# loc = c(0, 16, 38, 48, 58, 68, 88, 98, 108, 118, 128, 138, 148, 161, 198, 201, 211, 214, 234, 250, 273, 283, 291)
dist = diff (loc)

## Measurement error due to failure of bisulfite conversion
## (fixed throughout)
error.b = 0.00

## Print MCMC iteration numbers 
PRINT.ITER=TRUE


#################################################
### Tuning the MCMC run
#################################################
## Initial values of parameters
# site-specific methylation probability m
m.init = apply (data, 2, mean)
# associating probability tau for M, RP, RD
tau.curr = c(0.8, 0, 0)
# dissociating probability rho for M, RP, RD
rho.curr = c(0.1, 1, 1)
# mean parameter in the beta distribution 
# for site-specific methylation probability m
rm.curr = 0.8
# scaled variance parameter in the beta distribution 
# for site-specific methylation probability m
gm.curr = 0.01
# measurement error due to inappropriate bisulfite conversion 
c.curr = 0.01
# de novo activity rate of DNMT1 (on daughter strand)
m.denovo.curr = 0.05
# maintenace activity rate of DNMT1 (on daughter strand)
m.maint.curr = 0.9
# de novo activity rate of DNMT3s on daughter strand
d.denovo.curr = 0.8
# maintenance activity rate of DNMT3s on daughter strand
d.maint.curr = 0.8

## Standard deviations used to generate proposals of parameters
# sd for site-specific methylation probability m
sd.m = 0.01
# sd for associating probability tau; (M, RP, RD)
sd.tau = c(0.01, 0, 0)
# sd for dissociating probability rho; (M, RP, RD)
sd.rho = c(0.05, 0, 0)
# sd for mean parameter in beta distribution of m
sd.rm = 0.04
# sd for scaled variance parameter in beta distribution of m
sd.gm = 0.3
# sd for inappropriate bisulfite conversion error
sd.c = 0.005
# sd for de novo activity rate of DNMT1
sd.m.denovo = 0.08
# sd for maintenance activity rate of DNMT1
sd.m.maint = 0.04
# sd for de novo activity rate of DNMT3 on daughter strand
sd.d.denovo = 2
# sd for maintenance activity rate of DNMT3 on daughter strand
sd.d.maint = 2


#################################################
### Advanced options
#################################################
# whether (TRUE) or not (FALSE) to estimate 
# measurement error due to inappropriate bisulfite conversion
ESTIMATE.C = TRUE
# Lower and upper bound in the uniform prior 
# measurement error due to inappropriate bisulfite conversion
c.ub = 0.06
c.lb = 0
# Whether or not to estimate de novo and maintenance activity rates
# for a class of enzymes
DNMT1.EST = 0	# 1 for estimating these rates for DNMT1; 0 otherwise
DNMT3.EST = 0	# 1 for estimating these rates for the DNMT3s; 0 otherwise
## Priors 
# Choose type of prior for rho
# (options: uniform; log10 uniform, jeffreys)
RHO.PRIOR = "unif"
rho.lower=c(0,0,0)
rho.upper=c(1,1,1)
# Choose type of prior for tau
# (options: unif, logunif, jeffreys)
# and specify parameters in the prior if needed
TAU.PRIOR = "unif"
tau.lower = c(0, 0, 0)
tau.upper = c(1, 1, 1)



#################################################
### Compile R files and run program
#################################################
## compile the R files
source (paste(codeDir,"computeTransProbs.R", sep=""))
source (paste(codeDir,"computeTransProbsInd.R", sep=""))
source (paste(codeDir,"HMMD_MCMC.R", sep=""))
source (paste(codeDir,"loglikHMM.R", sep=""))
dyn.load (paste(codeDir,"loglikHMM.so", sep=""))
source (paste(codeDir,"mtx_exp.R", sep=""))
source (paste(codeDir,"reflect.R", sep=""))
source (paste(codeDir,"Sampler_c.R", sep=""))
source (paste(codeDir,"Sampler_ddenovo.R", sep=""))
source (paste(codeDir,"Sampler_dmaint.R", sep=""))
source (paste(codeDir,"Sampler_g.R", sep=""))
source (paste(codeDir,"Sampler_lambda.R", sep=""))
source (paste(codeDir,"Sampler_m.R", sep=""))
source (paste(codeDir,"Sampler_mdenovo.R", sep=""))
source (paste(codeDir,"Sampler_mmaint.R", sep=""))
source (paste(codeDir,"Sampler_r.R", sep=""))
source (paste(codeDir,"Sampler_rho_ind_jeffreys.R", sep=""))
source (paste(codeDir,"Sampler_rho_ind_log.R", sep=""))
source (paste(codeDir,"Sampler_rho_ind_unif.R", sep=""))
source (paste(codeDir,"Sampler_rho.R", sep=""))
source (paste(codeDir,"Sampler_tau_ind_jeffreys.R", sep=""))
source (paste(codeDir,"Sampler_tau_ind_log.R", sep=""))
source (paste(codeDir,"Sampler_tau_ind.R", sep=""))
source (paste(codeDir,"Sampler_tau.R", sep=""))
source (paste(codeDir,"StrandAssignmentProbs.R", sep=""))
dyn.load (paste(codeDir,"StrandAssignmentProbs.so", sep=""))

## Start the MCMC run
HMMD_MCMC (file.pars.mcmc, file.strandtype, data, dist, tau.curr, rho.curr, m.init, rm.curr, gm.curr, c.curr, m.denovo.curr, m.maint.curr, d.denovo.curr, d.maint.curr, error.b, n.iter, burn.in, step.size, sd.tau, sd.rho, sd.m, sd.rm, sd.gm, sd.c, sd.m.denovo, sd.m.maint, sd.d.denovo, sd.d.maint, tau.lower, tau.upper, rho.lower, rho.upper, c.ub, c.lb, seed.value, TAU.PRIOR, RHO.PRIOR, DNMT1.EST, DNMT3.EST, ERROR.c=ESTIMATE.C, DE.NOVO=1, PRINT.FLAG=0, PRINT.ITER)
##########################################################################################
# Command lines to run the MCMC program to estimate several properties of DNA methyltransferases
# using double-stranded DNA methylation patterns.  
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
file.in <- paste (dataDir, "Leptin1through22Fat.dat", sep="")
 
tmp = read.table (file.in, header = FALSE, sep=",")
dim (tmp)
meth.data = as.matrix (tmp[,-1])

## output files
file.pars.mcmc = paste (dataDir, "Leptin1through22Fat_pars_mcmc_1.out", sep="")
file.strandtype = paste (dataDir, "Leptin1through22Fat_strandtype_mcmc_1.out", sep="")

## MCMC specs
n.iter = 100
burn.in = 0.2
step.size = 5
seed.value = 18

## Nucleotide positions of CpG sites
loc = c(0, 13, 17, 22, 25, 40, 45, 55, 66, 69, 78, 89, 102, 107, 121, 139, 148, 153, 159, 178, 181, 183)
dist = diff (loc)

## Measurement error due to failure of bisulfite conversion
## (fixed throughout)
error.b = 0.003

## Print MCMC iteration numbers 
PRINT.ITER=TRUE


#################################################
### Tuning the MCMC run
#################################################
## Initial values of parameters
# site-specific methylation probability m
m.init = apply (meth.data, 2, mean)
# associating probability tau for M, RP, RD
tau.curr = c(0.5, 0.5, 0.5)
# dissociating probability rho for M, RP, RD
rho.curr = c(0.1, 0.5, 0.5)
# mean parameter in the beta distribution 
# for site-specific methylation probability m
rm.curr = 0.5
# scaled variance parameter in the beta distribution 
# for site-specific methylation probability m
gm.curr = 0.01
# measurement error due to inappropriate bisulfite conversion 
c.curr = 0.03
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
sd.tau = c(0.15, 0.05, 0.03)
# sd for dissociating probability rho; (M, RP, RD)
sd.rho = c(0.05, 0.9, 0.6)
# sd for mean parameter in beta distribution of m
sd.rm = 0.12
# sd for scaled variance parameter in beta distribution of m
sd.gm = 0.6
# sd for inappropriate bisulfite conversion error
sd.c = 0.01
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
c.ub = 0.1
c.lb = 0
# Whether or not to estimate de novo and maintenance activity rates
# for a class of enzymes
DNMT1.EST = 1	# 1 for estimating these rates for DNMT1; 0 otherwise
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
HMMD_MCMC (file.pars.mcmc, file.strandtype, meth.data, dist, tau.curr, rho.curr, m.init, rm.curr, gm.curr, c.curr, m.denovo.curr, m.maint.curr, d.denovo.curr, d.maint.curr, error.b, n.iter, burn.in, step.size, sd.tau, sd.rho, sd.m, sd.rm, sd.gm, sd.c, sd.m.denovo, sd.m.maint, sd.d.denovo, sd.d.maint, tau.lower, tau.upper, rho.lower, rho.upper, c.ub, c.lb, seed.value, TAU.PRIOR, RHO.PRIOR, DNMT1.EST, DNMT3.EST, ERROR.c=ESTIMATE.C, DE.NOVO=1, PRINT.FLAG=0, PRINT.ITER)
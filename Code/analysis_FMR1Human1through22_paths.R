####################################################################
# R code for inference of top two most likely explanations for each 
# pattern in the FMR1 data under the hidden Markov model (HMM).
#################################################################### 

############################################
# compile R code
############################################
codeDir = "Tools/"
source (paste(codeDir,"computeLatentEmiss.R", sep=""))
source (paste(codeDir,"computeLatentEmissError.R", sep=""))
source (paste(codeDir,"computeLatentMarginal.R", sep=""))
source (paste(codeDir,"computeLatentTrans.R", sep=""))
source (paste(codeDir,"computeTransProbs.R", sep=""))
source (paste(codeDir,"computeTransProbsind.R", sep=""))
source (paste(codeDir,"forward_backward.R", sep=""))
source (paste(codeDir,"integer_base_b.R", sep=""))
source (paste(codeDir,"loglikHMMlatent.R", sep=""))
source (paste(codeDir,"loglikHMMmethyl.R", sep=""))
source (paste(codeDir,"mtx_exp.R", sep=""))
source (paste(codeDir,"summary_onepath.R", sep=""))
source (paste(codeDir,"summary_paths.R", sep=""))
source (paste(codeDir,"Viterbi.R", sep=""))
source (paste(codeDir,"Viterbi2nd.R", sep=""))
source (paste(codeDir,"Viterbi2ndIndHMMmethyl.R", sep=""))
source (paste(codeDir,"Viterbi2ndPathProb.R", sep=""))
source (paste(codeDir,"ViterbiHMMmethyl.R", sep=""))



##################################################
# read in FMR1 data at sites 1-22 and MCMC output
##################################################
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

# read in mcmc output
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

tau.med = apply (tau.mcmc, 2, median)
rho.med = apply (rho.mcmc, 2, median)
m.med = apply (m.mcmc, 2, median)
mm.med = median (mm.mcmc)
md.med = median (md.mcmc)
c.med = median (c.mcmc)

file.pr = paste (mcmcDir, "FMR1Human22_hmm_v6_1_strandtype_mcmc.out", sep="")
assign.pr1 = as.matrix (read.table (file.pr, header = FALSE))
file.pr = paste (mcmcDir, "FMR1Human22_hmm_v6_1b_strandtype_mcmc.out", sep="")
assign.pr2 = as.matrix (read.table (file.pr, header = FALSE))
file.pr = paste (mcmcDir, "FMR1Human22_hmm_v6_1c_strandtype_mcmc.out", sep="")
assign.pr3 = as.matrix (read.table (file.pr, header = FALSE))
assign.pr = as.matrix (rbind (assign.pr1, assign.pr2, assign.pr3))
tpr = apply (assign.pr, 2, mean)


########################################################
# compute log likelihood and infer top two most likely
# explanations for each pattern
########################################################

##### compute log likelihood of each pattern
# take a few seconds
# h64: 64 hidden states at each site; each state is a binary vector (Q, D, P, M, RP, RD)
# Q: parent-strand CpG after replication being methylated (1) or unmethylated (0)
# D: daughter-strand CpG being methylated (1) or unmethylated (0)
# P: parent-strand CpG before replication being methylated (1) or unmethylated (0)
# M: DNMT1 associated (1) or unassociated (0) with DNA at this site
# RP: DNMT3s associated (1) or unassociated (0) with parent strand at this site
# RD: DNMT3s associated (1) or unassociated (0) with daughter strand at this site
logLh64 = loglikHMMmethyl (fmr1.data, dist, tau.med, rho.med, m.med, mm.med, md.med, c.med, b=0.003, model="h64")

# double-check results by computing the same quantity 
# under two other models with different numbers of hidden states
# h16: 16 hidden states at each site; each state is a binary vector (P, M, RP, RD)
# P: parent-strand CpG before replication being methylated (1) or unmethylated (0)
# M: DNMT1 associated (1) or unassociated (0) with DNA at this site
# RP: DNMT3s associated (1) or unassociated (0) with parent strand at this site
# RD: DNMT3s associated (1) or unassociated (0) with daughter strand at this site
#logLh16 = loglikHMMmethyl (fmr1.data, dist, tau.med, rho.med, m.med, mm.med, md.med, c.med, b=0.003, model="h16")

# h8: 8 hidden states at each site; each state is a binary vector (M, RP, RD)
# M: DNMT1 associated (1) or unassociated (0) with DNA at this site
# RP: DNMT3s associated (1) or unassociated (0) with parent strand at this site
# RD: DNMT3s associated (1) or unassociated (0) with daughter strand at this site
#logLh8 = loglikHMMmethyl (fmr1.data, dist, tau.med, rho.med, m.med, mm.med, md.med, c.med, b=0.003, model="h8")

##### compute Viterbi path (most likely) of each pattern
# 169 patterns
# take a few seconds
vpathh64 = ViterbiHMMmethyl (fmr1.data, dist, tau.med, rho.med, m.med, mm.med, md.med, c.med, b=0.003, loglik=logLh64, model="h64")

##### find 2nd best
# for pattern 159
k=159
v2nd.path = Viterbi2ndIndHMMmethyl (fmr1.data[(2*k-1):(2*k),], vpathh64$path[k,], dist, tau.med, rho.med, m.med, mm.med, md.med, c.med, b=0.003, logLh64[k], tpr[k], model="h64")
v2nd.path
top2.summary = summaryPaths (fmr1.data[(2*k-1):(2*k),], vpathh64$path[k,], vpathh64$condprob[k], v2nd.path)
top2.summary

##### results
# see annotation of output below
#######################################################################
k = 159
$parent		
# which strand, top or bottom, is inferred as the parent strand
[1] "b"

$stprob		
# the probability of strand assignment
V159 
0.0004179539 

$condprob	
# probability of an explanation given the pattern
# top two most likely probabilities are shown
[1] 0.10931474 0.03899977

$paths		
# top two most likely explanations
	[,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22]
1   61    5    5   61   61   61   41   41   41    41    18    18    58    41    61    61    61    61     5    61    61    61
2   61    5    5   61   61   61   41   41   41    41    61    61    61    61    61    61    61    61     5    61    61    61

$converted			
# explanations converted to interpretable format
# each interpretation is converted to (Q, D, P, M, RP, RD)
# Q,D,P are methylation states on post-replication parent strand, daughter strand and pre-replication parent strand
# M, RP, RD are association states of DNMT1 on daughter strand, DNMT3s on parent strand and DNMT3s on daughter strand
   1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
x  1 0 0 1 1 1 0 0 0  0  1  1  1  0  1  1  1  1  0  1  1  1
y  1 0 0 1 1 1 1 1 1  1  0  0  1  1  1  1  1  1  0  1  1  1
Q  1 0 0 1 1 1 1 1 1  1  0  0  1  1  1  1  1  1  0  1  1  1
D  1 0 0 1 1 1 0 0 0  0  1  1  1  0  1  1  1  1  0  1  1  1
P  1 0 0 1 1 1 1 1 1  1  0  0  1  1  1  1  1  1  0  1  1  1
M  1 1 1 1 1 1 0 0 0  0  0  0  0  0  1  1  1  1  1  1  1  1
RP 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
RD 0 0 0 0 0 0 0 0 0  0  1  1  1  0  0  0  0  0  0  0  0  0
Q  1 0 0 1 1 1 1 1 1  1  1  1  1  1  1  1  1  1  0  1  1  1
D  1 0 0 1 1 1 0 0 0  0  1  1  1  1  1  1  1  1  0  1  1  1
P  1 0 0 1 1 1 1 1 1  1  1  1  1  1  1  1  1  1  0  1  1  1
M  1 1 1 1 1 1 0 0 0  0  1  1  1  1  1  1  1  1  1  1  1  1
RP 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
RD 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
#######################################################################

#######################################################################
k = 165
$parent
[1] "b"

$stprob
V165 
0.002359515 

$condprob
[1] 0.34228170 0.04796458

$paths
  [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22]
1   61   61   61   61   61   41   41   41   41    61     5    61    61    61    61    61     5    61    61    61    61     5
2   61   61   61   61   61   61   41   41   41    61     5    61    61    61    61    61     5    61    61    61    61     5

$converted
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
x  1 1 1 1 1 0 0 0 0  1  0  1  1  1  1  1  0  1  1  1  1  0
y  1 1 1 1 1 1 1 1 1  1  0  1  1  1  1  1  0  1  1  1  1  0
Q  1 1 1 1 1 1 1 1 1  1  0  1  1  1  1  1  0  1  1  1  1  0
D  1 1 1 1 1 0 0 0 0  1  0  1  1  1  1  1  0  1  1  1  1  0
P  1 1 1 1 1 1 1 1 1  1  0  1  1  1  1  1  0  1  1  1  1  0
M  1 1 1 1 1 0 0 0 0  1  1  1  1  1  1  1  1  1  1  1  1  1
RP 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
RD 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
Q  1 1 1 1 1 1 1 1 1  1  0  1  1  1  1  1  0  1  1  1  1  0
D  1 1 1 1 1 1 0 0 0  1  0  1  1  1  1  1  0  1  1  1  1  0
P  1 1 1 1 1 1 1 1 1  1  0  1  1  1  1  1  0  1  1  1  1  0
M  1 1 1 1 1 1 0 0 0  1  1  1  1  1  1  1  1  1  1  1  1  1
RP 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
RD 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
#######################################################################

#######################################################################
k=60
$parent
[1] "b"

$stprob
V60 
0.0003057218 

$condprob
[1] 0.05522696 0.03778380

$paths
  [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22]
1   41   41   61    5   61   61   61   61   61    61     5    61    61    41    41    41    41    61    61    61    61    61
2   41   41   61    5   61   61   61   61   45    61     5    61    61    41    41    41    41    61    61    61    61    61

$converted
   1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
x  0 0 1 0 1 1 1 1 0  1  0  1  1  0  0  0  0  1  1  1  1  1
y  1 1 1 0 1 1 1 1 1  1  0  1  1  1  1  1  1  1  1  1  1  1
Q  1 1 1 0 1 1 1 1 1  1  0  1  1  1  1  1  1  1  1  1  1  1
D  0 0 1 0 1 1 1 1 1  1  0  1  1  0  0  0  0  1  1  1  1  1
P  1 1 1 0 1 1 1 1 1  1  0  1  1  1  1  1  1  1  1  1  1  1
M  0 0 1 1 1 1 1 1 1  1  1  1  1  0  0  0  0  1  1  1  1  1
RP 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
RD 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
Q  1 1 1 0 1 1 1 1 1  1  0  1  1  1  1  1  1  1  1  1  1  1
D  0 0 1 0 1 1 1 1 0  1  0  1  1  0  0  0  0  1  1  1  1  1
P  1 1 1 0 1 1 1 1 1  1  0  1  1  1  1  1  1  1  1  1  1  1
M  0 0 1 1 1 1 1 1 1  1  1  1  1  0  0  0  0  1  1  1  1  1
RP 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
RD 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
#######################################################################

#######################################################################
k = 82
$parent
[1] "b"

$stprob
V82 
0.003562172 

$condprob
[1] 0.42700386 0.02181886

$paths
  [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22]
1   61    5   61    5    5   61   61   61   61    61    61     5    61    61    61    61    61    41     1    41    41    41
2   61    5   61    5    5   61   61   61   61    61    61     5    61    61    61    61    61    41    41    41    41    41

$converted
   1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
x  1 0 1 0 0 1 1 1 1  1  1  0  1  1  1  1  1  0  0  0  0  0
y  1 0 1 0 0 1 1 1 1  1  1  0  1  1  1  1  1  1  0  1  1  1
Q  1 0 1 0 0 1 1 1 1  1  1  0  1  1  1  1  1  1  0  1  1  1
D  1 0 1 0 0 1 1 1 1  1  1  0  1  1  1  1  1  0  0  0  0  0
P  1 0 1 0 0 1 1 1 1  1  1  0  1  1  1  1  1  1  0  1  1  1
M  1 1 1 1 1 1 1 1 1  1  1  1  1  1  1  1  1  0  0  0  0  0
RP 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
RD 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
Q  1 0 1 0 0 1 1 1 1  1  1  0  1  1  1  1  1  1  1  1  1  1
D  1 0 1 0 0 1 1 1 1  1  1  0  1  1  1  1  1  0  0  0  0  0
P  1 0 1 0 0 1 1 1 1  1  1  0  1  1  1  1  1  1  1  1  1  1
M  1 1 1 1 1 1 1 1 1  1  1  1  1  1  1  1  1  0  0  0  0  0
RP 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
RD 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
#######################################################################

#######################################################################
k=83
$parent
[1] "t"

$stprob
V83 
0.8597175 

$condprob
[1] 0.02487117 0.01701573

$paths
  [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22]
1   61    5   61   61   61    5   61    5   61    61    61     5    61    61    61    61    61    61     5    61    61    61
2   45    5   61   61   61    5   61    5   61    61    61     5    61    61    61    61    61    61     5    61    61    61
3   61    5   61   61   61    5   61    5   61    45    61     5    61    61    61    61    61    61     5    61    61    61
4   61    5   61   61   61    5   61    5   61    61    61     5    45    61    61    61    61    61     5    61    61    61
5   61    5   61   61   61    5   61    5   61    61    61     5    61    61    61    61    45    61     5    61    61    61

$converted
   1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
x  1 0 1 1 1 0 1 0 1  1  1  0  1  1  1  1  1  1  0  1  1  1
y  0 0 1 1 1 0 1 0 1  0  1  0  0  1  1  1  0  1  0  1  1  1

Q  1 0 1 1 1 0 1 0 1  1  1  0  1  1  1  1  1  1  0  1  1  1
D  1 0 1 1 1 0 1 0 1  1  1  0  1  1  1  1  1  1  0  1  1  1
P  1 0 1 1 1 0 1 0 1  1  1  0  1  1  1  1  1  1  0  1  1  1
M  1 1 1 1 1 1 1 1 1  1  1  1  1  1  1  1  1  1  1  1  1  1
RP 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
RD 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0

Q  1 0 1 1 1 0 1 0 1  1  1  0  1  1  1  1  1  1  0  1  1  1
D  0 0 1 1 1 0 1 0 1  1  1  0  1  1  1  1  1  1  0  1  1  1
P  1 0 1 1 1 0 1 0 1  1  1  0  1  1  1  1  1  1  0  1  1  1
M  1 1 1 1 1 1 1 1 1  1  1  1  1  1  1  1  1  1  1  1  1  1
RP 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
RD 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0

Q  1 0 1 1 1 0 1 0 1  1  1  0  1  1  1  1  1  1  0  1  1  1
D  1 0 1 1 1 0 1 0 1  0  1  0  1  1  1  1  1  1  0  1  1  1
P  1 0 1 1 1 0 1 0 1  1  1  0  1  1  1  1  1  1  0  1  1  1
M  1 1 1 1 1 1 1 1 1  1  1  1  1  1  1  1  1  1  1  1  1  1
RP 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
RD 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0

Q  1 0 1 1 1 0 1 0 1  1  1  0  1  1  1  1  1  1  0  1  1  1
D  1 0 1 1 1 0 1 0 1  1  1  0  0  1  1  1  1  1  0  1  1  1
P  1 0 1 1 1 0 1 0 1  1  1  0  1  1  1  1  1  1  0  1  1  1
M  1 1 1 1 1 1 1 1 1  1  1  1  1  1  1  1  1  1  1  1  1  1
RP 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
RD 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0

Q  1 0 1 1 1 0 1 0 1  1  1  0  1  1  1  1  1  1  0  1  1  1
D  1 0 1 1 1 0 1 0 1  1  1  0  1  1  1  1  0  1  0  1  1  1
P  1 0 1 1 1 0 1 0 1  1  1  0  1  1  1  1  1  1  0  1  1  1
M  1 1 1 1 1 1 1 1 1  1  1  1  1  1  1  1  1  1  1  1  1  1
RP 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
RD 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
#######################################################################


#######################################################################
# number of hemi sites: 8
k=5
$parent
[1] "t"

$stprob
V5 
0.5505707 

$condprob
[1] 0.005893030 0.004031746

$paths
  [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22]
1    5    5   61   61   61   61   61   61   61    61    41    41    61    61    61    61    61    22    22    61    61    61
2    5    5   61   61   61   61   61   61   61    61    41    41    61    61    61    61    61    22    22    61    45    61

$converted
   1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
x  0 0 1 0 1 1 1 0 1  1  1  1  1  1  1  1  1  0  0  1  1  1
y  0 0 0 1 1 1 1 1 1  1  0  0  1  1  1  1  1  1  1  1  0  1

Q  0 0 1 1 1 1 1 1 1  1  1  1  1  1  1  1  1  0  0  1  1  1
D  0 0 1 1 1 1 1 1 1  1  0  0  1  1  1  1  1  1  1  1  1  1
P  0 0 1 1 1 1 1 1 1  1  1  1  1  1  1  1  1  0  0  1  1  1
M  1 1 1 1 1 1 1 1 1  1  0  0  1  1  1  1  1  1  1  1  1  1
RP 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
RD 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  1  1  0  0  0

Q  0 0 1 1 1 1 1 1 1  1  1  1  1  1  1  1  1  0  0  1  1  1
D  0 0 1 1 1 1 1 1 1  1  0  0  1  1  1  1  1  1  1  1  0  1
P  0 0 1 1 1 1 1 1 1  1  1  1  1  1  1  1  1  0  0  1  1  1
M  1 1 1 1 1 1 1 1 1  1  0  0  1  1  1  1  1  1  1  1  1  1
RP 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
RD 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  1  1  0  0  0
#######################################################################

#######################################################################
# number of hemi sites: 1
k=168
$parent
[1] "t"

$stprob
V168 
0.5070984 

$condprob
[1] 0.08869964 0.06068430

$paths
  [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22]
1   61   61   61   61   61   61   61   61   61    61    61    61    61    61    61    61    61    61    61    61    61    61
2   61   45   61   61   61   61   61   61   61    61    61    61    61    61    61    61    61    61    61    61    61    61

$converted
   1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
x  1 1 1 1 1 1 1 1 1  1  1  1  1  1  1  1  1  1  1  1  1  1
y  1 0 1 1 1 1 1 1 1  1  1  1  1  1  1  1  1  1  1  1  1  1
Q  1 1 1 1 1 1 1 1 1  1  1  1  1  1  1  1  1  1  1  1  1  1
D  1 1 1 1 1 1 1 1 1  1  1  1  1  1  1  1  1  1  1  1  1  1
P  1 1 1 1 1 1 1 1 1  1  1  1  1  1  1  1  1  1  1  1  1  1
M  1 1 1 1 1 1 1 1 1  1  1  1  1  1  1  1  1  1  1  1  1  1
RP 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
RD 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
Q  1 1 1 1 1 1 1 1 1  1  1  1  1  1  1  1  1  1  1  1  1  1
D  1 0 1 1 1 1 1 1 1  1  1  1  1  1  1  1  1  1  1  1  1  1
P  1 1 1 1 1 1 1 1 1  1  1  1  1  1  1  1  1  1  1  1  1  1
M  1 1 1 1 1 1 1 1 1  1  1  1  1  1  1  1  1  1  1  1  1  1
RP 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
RD 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0
#######################################################################



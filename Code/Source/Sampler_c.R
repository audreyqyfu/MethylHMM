#################################################################################
# Sampler_c
#
# 2010.09.12: Bug: add argument dist and pass on to loglikHMM.
# 2010.09.10: lambda = c(tau, rho, m) is now pass to loglikHMM.
# 2010.04.16: allow rates of maintanence and de novo activities 
#			  for DNMT3 to be different
# 2009.09.08: modified for a single data file;
#			  reflect proposals outside the support.
# 2009.08.11: argument lambda changes to lambda.long
# 2009.08.09: parametrization as in Note_10_EmissProbs
# 2008.09.08: change to relaxed HMM on maintenance enzyme
#
# Update c, constant across sites.
# 
# input: c.curr: current value of c
#        data: unprocessed double-stranded data; matrix of N by S
#		 dist: distances between adjacent CpG sites; vector of length S-1.
#        lambda: 3 tau (M, RP, RD), 3 rho (M, RP, RD) and S m
#		 m.denovo: de novo rate of Dnmt1
#		 m.maint: rate of maintenance activity of Dnmt1
#		 d.denovo.curr: de novo rate of Dnmt3 on daughter strand
#		 d.maint.curr: rate of maintenance activity of Dnmt3 on daughter strand
#		 b: failure of conversion rate
#        sd.c: sd in the transition kernel to generate 
#              proposal for c
#        lower: lower bound on c
#        upper: upper bound on c
#        DE.NOVO: 0 for assumption 1, 1 for assumption 2
#        PRINT.FLAG: 0 for no output, 1 for output
# output: result: updated value of c
####################################################################################

Sampler_c = function (c.curr, data, dist, lambda, m.denovo, m.maint, d.denovo.curr, d.maint.curr, b, sd.c, lower, upper, DE.NOVO=1, PRINT.FLAG=0)
{
	n.proc = 3
	
	c.tmp = rnorm (1, c.curr, sd.c)
	c.new = ifelse (c.tmp>upper | c.tmp<lower, reflect (c.tmp,lower, upper), c.tmp)

	###################################################
    # compute likelihood 
    ###################################################
    # under current value of lambda
    loglik.curr = loglikHMM (data, dist, lambda, m.denovo, m.maint, d.denovo.curr, d.maint.curr, b, c.curr, DE.NOVO, PRINT.FLAG)
    # under proposed value of lambda
	loglik.new = loglikHMM (data, dist, lambda, m.denovo, m.maint, d.denovo.curr, d.maint.curr, b, c.new, DE.NOVO, PRINT.FLAG)

    #################################################
    # compute Hastings ratio and decide to accept or 
    # reject the proposal
    #################################################
	log.ratio = loglik.new - loglik.curr
	result = ifelse ((log.ratio!='NaN') & (log(runif(1))<log.ratio), c.new, c.curr)

	return (result)
}

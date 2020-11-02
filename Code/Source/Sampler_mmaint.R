#############################################################################
# Sampler_mmaint
#
# 2010.09.12: Bug: add argument dist and pass on to loglikHMM.
# 2010.09.10: lambda = c(tau, rho, m) is now pass to loglikHMM.
# 2010.04.16: allow rates of maintanence and de novo activities 
#			  for DNMT3 to be different
# 2009.09.09: modified to analyze only one input data file;
#			  reflect proposals into the support.
# 2009.08.11: argument lambda changes to lambda.long
# 2009.08.09: parametrization as in Note_10_EmissProbs
# 2008.09.08: MCMC for the two data sets together
# 
# update m.maint, rate of maintenance activity of Dnmt1
# 
# input: m.maint.curr: current value of m.maint
#        data: unprocessed double-stranded data
#		 dist: distances between adjacent CpG sites; vector of length S-1.
#        lambda: 3 tau (M, RP, RD), 3 rho (M, RP, RD) and S m
#		 m.denovo: de novo rate of Dnmt1
#		 d.denovo: de novo rate of Dnmt3 on daughter strand
#		 d.maint: rate of maintenance activity of Dnmt3 on daughter strand
#		 b: failure of conversion rate 
#		 c: inappropriate conversion rate 
#        sd.m.maint: sd in the proposal distribution of m.maint
#		 lower: lower bound on m.maint
#		 upper: upper bound on m.maint
#        DE.NOVO: 0 for assumption 1, 1 for assumption 2
#        PRINT.FLAG: 0 for no output, 1 for output
# output: result: updated value of m.maint
###############################################################################

Sampler_mmaint = function (m.maint.curr, data, dist, lambda, m.denovo, d.denovo, d.maint, b, c, sd.m.maint, lower, upper, DE.NOVO=1, PRINT.FLAG=0)
{
	m.maint.tmp = rnorm (1, m.maint.curr, sd.m.maint)
	m.maint.new = ifelse (m.maint.tmp>upper | m.maint.tmp<lower, reflect (m.maint.tmp, lower, upper), m.maint.tmp)
	
	###################################################
    # compute likelihood 
    ###################################################
    # under current value of lambda
    loglik.curr = loglikHMM (data, dist, lambda, m.denovo, m.maint.curr, d.denovo, d.maint, b, c, DE.NOVO, PRINT.FLAG)
    # under proposed value of lambda
    loglik.new = loglikHMM (data, dist, lambda, m.denovo, m.maint.new, d.denovo, d.maint, b, c, DE.NOVO, PRINT.FLAG)

    #################################################
    # compute Hastings ratio and decide to accept or 
    # reject the proposal
    #################################################
	log.ratio = loglik.new - loglik.curr
	result = ifelse ((log.ratio!='NaN') & (log(runif(1))<log.ratio), m.maint.new, m.maint.curr)

	return (result)
}

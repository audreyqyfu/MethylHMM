#############################################################################
# Sampler_dmaint
#
# 2010.09.12: Bug: add argument dist and pass on to loglikHMM.
# 2010.09.10: lambda = c(tau, rho, m) is now pass to loglikHMM.
# 2010.04.16: allow rates of maintanence and de novo activities 
#			  for DNMT3 to be different
# 
# update d.maint, rate of maintenance activity of Dnmt3 on daughter strand
# 
# input: d.maint.curr: current value of d.maint
#        data: unprocessed double-stranded data
#		 dist: distances between adjacent CpG sites; vector of length S-1.
#        lambda: 3 tau (M, RP, RD), 3 rho (M, RP, RD) and S m
#		 d.denovo: de novo rate of Dnmt3 on daughter strand
#		 m.denovo: de novo rate of Dnmt1
#		 m.maint: rate of maintenance activity of Dnmt1
#		 b: failure of conversion rate 
#		 c: inappropriate conversion rate 
#        sd.m.maint: sd in the proposal distribution of m.maint
#		 lower: lower bound on d.maint
#		 upper: upper bound on d.maint
#        DE.NOVO: 0 for assumption 1, 1 for assumption 2
#        PRINT.FLAG: 0 for no output, 1 for output
# output: result: updated value of d.maint
###############################################################################

Sampler_dmaint = function (d.maint.curr, data, dist, lambda, d.denovo, m.denovo, m.maint, b, c, sd.d.maint, lower, upper, DE.NOVO=1, PRINT.FLAG=0)
{
	d.maint.tmp = rnorm (1, d.maint.curr, sd.d.maint)
	d.maint.new = ifelse (d.maint.tmp>upper | d.maint.tmp<lower, reflect (d.maint.tmp, lower, upper), d.maint.tmp)
	
	###################################################
    # compute likelihood 
    ###################################################
    # under current value of lambda
    loglik.curr = loglikHMM (data, dist, lambda, m.denovo, m.maint, d.denovo, d.maint.curr, b, c, DE.NOVO, PRINT.FLAG)
    # under proposed value of lambda
    loglik.new = loglikHMM (data, dist, lambda, m.denovo, m.maint, d.denovo, d.maint.new, b, c, DE.NOVO, PRINT.FLAG)

    #################################################
    # compute Hastings ratio and decide to accept or 
    # reject the proposal
    #################################################
	log.ratio = loglik.new - loglik.curr
	result = ifelse ((log.ratio!='NaN') & (log(runif(1))<log.ratio), d.maint.new, d.maint.curr)

	return (result)
}

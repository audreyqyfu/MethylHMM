#######################################################
# Sampler_mdenovo
#
# 2010.09.12: Bug: add argument dist and pass on to loglikHMM.
# 2010.09.10: lambda = c(tau, rho, m) is now pass to loglikHMM.
# 2010.04.16: allow rates of maintanence and de novo activities 
#			  for DNMT3 to be different
# 2009.09.09: modified for a single input data file.
# 2009.08.11: argument lambda changes to lambda.long
# 2009.08.09: parametrization as in Note_10_EmissProbs
# 2008.09.08: MCMC for two data sets together
#
# update m.denovo, de novo rate of Dnmt1 at previously 
# unmethylated sites
# 
# input: m.denovo.curr: current value of m.denovo
#        data: unprocessed double-stranded data
#		 dist: distances between adjacent CpG sites; vector of length S-1.
#        lambda: 3 tau (M, RP, RD), 3 rho (M, RP, RD) and S m
#		 m.maint: rate of maintenance activity of Dnmt1
#		 d.denovo: de novo rate of Dnmt3 on daughter strand
#		 d.maint: rate of maintenance activity of Dnmt3 on daughter strand
#		 b: failure of conversion rate 
#		 c: inappropriate conversion rate 
#        sd.m.denovo: sd in the transition kernel to generate 
#              proposal for m.denovo
#		 lower: lower bound on m.denovo
#		 upper: upper bound on m.denovo
#        DE.NOVO: 0 for assumption 1, 1 for assumption 2
#        PRINT.FLAG: 0 for no output, 1 for output
# output: result: updated value of m.denovo
#######################################################

Sampler_mdenovo = function (m.denovo.curr, data, dist, lambda, m.maint, d.denovo, d.maint, b, c, sd.m.denovo, lower, upper, DE.NOVO=1, PRINT.FLAG=0)
{
	n.proc = 3
	
	m.denovo.tmp = rnorm (1, m.denovo.curr, sd.m.denovo)
	m.denovo.new = ifelse (m.denovo.tmp>upper | m.denovo.tmp<lower, reflect (m.denovo.tmp, lower, upper), m.denovo.tmp)

	###################################################
    # compute likelihood 
    ###################################################
    # under current value of lambda
    loglik.curr = loglikHMM (data, dist, lambda, m.denovo.curr, m.maint, d.denovo, d.maint, b, c, DE.NOVO, PRINT.FLAG)
	# under proposed value of lambda
    loglik.new = loglikHMM (data, dist, lambda, m.denovo.new, m.maint, d.denovo, d.maint, b, c, DE.NOVO, PRINT.FLAG)

    #################################################
    # compute Hastings ratio and decide to accept or 
    # reject the proposal
    #################################################
	log.ratio = loglik.new - loglik.curr
	result = ifelse ((log.ratio!='NaN') & (log(runif(1))<log.ratio), m.denovo.new, m.denovo.curr)

	return (result)
}

#######################################################
# Sampler_ddenovo
#
# 2010.09.12: Bug: add argument dist and pass on to loglikHMM.
# 2010.09.10: lambda = c(tau, rho, m) is now pass to loglikHMM.
# 2010.04.16: allow rates of maintanence and de novo activities 
#			  for DNMT3 to be different
#
# update d.denovo, de novo rate of Dnmt3 at previously 
# unmethylated sites on daughter strand
# 
# input: d.denovo.curr: current value of d.denovo
#        data: unprocessed double-stranded data
#		 dist: distances between adjacent CpG sites; vector of length S-1.
#        lambda: 3 tau (M, RP, RD), 3 rho (M, RP, RD) and S m
#		 d.maint: rate of maintenance activity of Dnmt3 on daughter strand
#		 m.denovo: de novo rate of Dnmt1
#		 m.maint: rate of maintenance activity of Dnmt1
#		 b: failure of conversion rate 
#		 c: inappropriate conversion rate 
#        sd.m.denovo: sd in the transition kernel to generate 
#              proposal for m.denovo
#		 lower: lower bound on d.denovo
#		 upper: upper bound on d.denovo
#        DE.NOVO: 0 for assumption 1, 1 for assumption 2
#        PRINT.FLAG: 0 for no output, 1 for output
# output: result: updated value of d.denovo
#######################################################

Sampler_ddenovo = function (d.denovo.curr, data, dist, lambda, d.maint, m.denovo, m.maint, b, c, sd.d.denovo, lower, upper, DE.NOVO=1, PRINT.FLAG=0)
{	
	d.denovo.tmp = rnorm (1, d.denovo.curr, sd.d.denovo)
	d.denovo.new = ifelse (d.denovo.tmp>upper | d.denovo.tmp<lower, reflect (d.denovo.tmp, lower, upper), d.denovo.tmp)

	###################################################
    # compute likelihood 
    ###################################################
    # under current value of lambda
    loglik.curr = loglikHMM (data, dist, lambda, m.denovo, m.maint, d.denovo.curr, d.maint, b, c, DE.NOVO, PRINT.FLAG)
	# under proposed value of lambda
    loglik.new = loglikHMM (data, dist, lambda, m.denovo, m.maint, d.denovo.new, d.maint, b, c, DE.NOVO, PRINT.FLAG)

    #################################################
    # compute Hastings ratio and decide to accept or 
    # reject the proposal
    #################################################
	log.ratio = loglik.new - loglik.curr
	result = ifelse ((log.ratio!='NaN') & (log(runif(1))<log.ratio), d.denovo.new, d.denovo.curr)

	return (result)
}

############################################################################
# Sampler_m
#
# 2010.09.12: Pass argument dist to loglikHMM.
# 2010.09.10: Discrete model representation.
# 2010.04.16: allow rates of maintanence and de novo activities 
#			  for DNMT3 to be different
# 2009.09.09: reflect proposals outside (0,1); lower and upper bounds in the 
#			  uniform prior are now arguments.
# 2009.08.11: return m only, rather than with pi and rho
# 2009.08.09: parametrization as in Note_10_EmissProbs
# 2008.09.08: incorporate de novo and spontaneous failure rates into HMM
#
# Function to update all m simutaneously for one data set.
#
# input: m.curr: current values of m; vector of length S
#		 data: of one region; matrix of N by S
#		 tau: vector of stationary probabilities; length 3
#		 rho: vector of jump rates; length 3
#		 dist: physical distance between CpG sites in nt (nucleotide)
#		 m.denovo: de novo rate of Dnmt1
#		 m.maint: rate of maintenance activity of Dnmt1
#		 d.denovo: de novo rate of Dnmt3 on daughter strand
#		 d.maint: rate of maintenance activity of Dnmt3 on daughter strand
#		 b: failure of conversion rate
#		 c: inappropriate conversion rate
#		 r.m: mean of the beta distn for m
#		 g.m: scaled variance in the beta distn for m
#		 sd.m: sd in the proposal distn for m
#		 DE.NOVO: assumption on de novo events
#		 PRINT.FLAG: 1 for output, and 0 o.w.
#		 accept.m: acceptance rate of m
# output: result: list containing updated m and acceptance rate of m
#############################################################################
Sampler_m = function (m.curr, data, tau, rho, dist, m.denovo, m.maint, d.denovo, d.maint, b, c, r.m, g.m, sd.m, lower, upper, accept.m, DE.NOVO=1, PRINT.FLAG=0)
{
	n.m = ncol(data)
	lambda.curr = c(tau, rho, m.curr)
	
	alpha.m = r.m*(1-g.m)/g.m
	beta.m = (1-r.m)*(1-g.m)/g.m

	# sample new m
	m.tmp = rnorm (n.m, m.curr, sd.m)
	m.new = rep (0, n.m)
	for (i in 1:n.m)
		m.new[i] = ifelse (m.tmp[i]>upper | m.tmp[i]<lower, reflect (m.tmp[i],lower,upper), m.tmp[i])

	lambda.new = c(tau, rho, m.new)
	
	###################################################
	# compute likelihood 
	###################################################
	# under current value of lambda
	loglik.curr = loglikHMM (data, dist, lambda.curr, m.denovo, m.maint, d.denovo, d.maint, b, c, DE.NOVO, PRINT.FLAG) + sum (dbeta(m.curr, alpha.m, beta.m, log=TRUE))
	# under proposed value of lambda
	loglik.new = loglikHMM (data, dist, lambda.new, m.denovo, m.maint, d.denovo, d.maint, b, c, DE.NOVO, PRINT.FLAG) + sum (dbeta(m.new, alpha.m, beta.m, log=TRUE))

	#################################################
	# compute Hastings ratio and decide to accept or 
	# reject the proposal
	#################################################
	log.ratio = loglik.new - loglik.curr
	if ((log.ratio!='NaN') & (log(runif(1))<log.ratio))
		result = list(m=m.new, accptm=accept.m+1)
	else
		result = list (m=m.curr, accptm=accept.m)

	return (result)
}

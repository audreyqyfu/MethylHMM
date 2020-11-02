################################################################################
# Sampler_tau_RPRD_ind
#
# 2011.10.10: Modified to estimate tau_RP=tau_RD.
# 2010.09.12: Add dist to argument and pass on to loglikHMM.
#
# Modified from Sampler_pi_ind.R.
#
# Function to update a single tau (Pr(0->1) in 1-step transition matrix)
# with a uniform prior.
#
# lambda={3 tau; 3 rho; m, vector of length S}
#
# input: i: the i-th tau to be updated
#		 tau.curr: current value of tau; vector of length 3
#		 n.proc: number of hidden processes (=3)
#        data: unprocessed double-stranded data; matrix of N by S
#		 dist: distances between adjacent CpG sites; vector of length S-1
#		 rho: jump rates; vector of length 3
#		 m.denovo: de novo rate of Dnmt1
#		 m.maint: rate of maintenance activity of Dnmt1
#		 d.denovo: de novo rate of Dnmt3 on daughter strand
#		 d.maint: rate of maintenance activity of Dnmt3 on daughter strand
#        b: failure of conversion rate
#        c: inappropriate conversion rate
#        r.m: parameter in the beta distn of m
#        g.m: parameter in the beta distn of m
#        sd.tau: sd in the normal transition kernel; 
#		 lower: lower bound of the uniform prior for tau
#		 upper: upper bound of the uniform prior for tau
#        DE.NOVO: 0 for assumption 1, 1 for assumption 2
#        PRINT.FLAG: 0 for no output, 1 for output
# output: result: tau vector with updated element
##################################################################################

Sampler_tau_RPRD_ind = function (tau.curr, n.proc, data, dist, rho, m, m.denovo, m.maint, d.denovo, d.maint, b, c, sd.tau, lower, upper, DE.NOVO=1, PRINT.FLAG=0)
{
	lambda.curr = c(as.vector (tau.curr), as.vector(rho), as.vector(m))
	lambda.new = lambda.curr
	par.new = rnorm (1, tau.curr[2], sd.tau)
	lambda.new[2:3] = ifelse (par.new>upper | par.new<lower, reflect(par.new, lower, upper), par.new)
	
	if (PRINT.FLAG)
	{
		print ("current lambda:")
		print (lambda.curr)
		print ("new lambda:")
		print (lambda.new)
	}
	
	###################################################
	# compute likelihood 
	###################################################
	# under current value of lambda (with theta instead of rho)
	loglik.curr = loglikHMM (data, dist, lambda.curr, m.denovo, m.maint, d.denovo, d.maint, b, c, DE.NOVO, PRINT.FLAG)
	# under proposed value of lambda (with theta instead of rho)
	loglik.new = loglikHMM (data, dist, lambda.new, m.denovo, m.maint, d.denovo, d.maint, b, c, DE.NOVO, PRINT.FLAG)
			
	#################################################
	# compute Hastings ratio and decide to accept or 
	# reject the proposal
	#################################################
	log.ratio = loglik.new - loglik.curr
	if ((log.ratio!='NaN') & (log(runif(1))<log.ratio))
		result = lambda.new[1:n.proc]
	else
		result = tau.curr
	
	if (PRINT.FLAG)
	{
		cat ("log lik using current lambda: ", loglik.curr, "\n")
		cat ("log lik using new lambda: ", loglik.new, "\n")
		cat ("difference: ", log.ratio, "\n")
		cat ("result: ", result, "\n")
	}
	
	return (result)
}

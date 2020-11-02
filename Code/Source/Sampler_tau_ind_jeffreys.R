#######################################################
# Sampler_tau_ind_jeffreys
#
# 2010.09.12: Add dist to argument and pass on to loglikHMM.
#
# Modified from Sampler_pi_ind_jeffreys.R.
#
# Function to update a single tau (Pr(0->1) in 1-step transition matrix)
# with a Jeffreys prior (beta (1/2, 1/2)).
#
# lambda={3 tau, 3 rho, m of S}
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
#        DE.NOVO: 0 for assumption 1, 1 for assumption 2
#        PRINT.FLAG: 0 for no output, 1 for output
# output: result: tau vector with updated element
#######################################################

Sampler_tau_ind_jeffreys = function (i, tau.curr, n.proc, data, dist, rho, m, m.denovo, m.maint, d.denovo, d.maint, b, c, sd.tau, DE.NOVO=1, PRINT.FLAG=0)
{
	lambda.curr = c(as.vector (tau.curr), as.vector(rho), as.vector(m))
	lambda.new = lambda.curr
	par.new = rnorm (1, tau.curr[i], sd.tau)
	taui.new = ifelse (par.new>1 | par.new<0, reflect (par.new,0,1), par.new)
	lambda.new[i] = taui.new
	
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
	loglik.curr = loglikHMM (data, dist, lambda.curr, m.denovo, m.maint, d.denovo, d.maint, b, c, DE.NOVO, PRINT.FLAG) + dbeta (tau.curr[i], 1/2,1/2,log=TRUE)
	# under proposed value of lambda (with theta instead of rho)
	loglik.new = loglikHMM (data, dist, lambda.new, m.denovo, m.maint, d.denovo, d.maint, b, c, DE.NOVO, PRINT.FLAG) + dbeta (taui.new, 1/2,1/2,log=TRUE)
			
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

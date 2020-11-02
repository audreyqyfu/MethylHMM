#######################################################
# Sampler_tau_ind_log
#
# 2010.09.30: Bug fix: should treat log10 transformed rho
#			  as a new random variable and update this r.v.
#			  for the uniform prior.  Should not transform 
#			  it back to the raw scale for updating!
# 2010.09.12: Add dist to argument and pass on to loglikHMM.
# 2010.09.10: Target log density modified because prior density 
#			  is log10 uniform rather than uniform.
#
# Modfified from Sampler_pi_ind_log.
# 
# Function to update a single tau (Pr(0->1) in 1-step transition matrix)
# with a log uniform prior.
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
#		 lower: lower bound of the uniform prior for tau
#		 upper: upper bound of the uniform prior for tau
#        DE.NOVO: 0 for assumption 1, 1 for assumption 2
#        PRINT.FLAG: 0 for no output, 1 for output
# output: result: tau vector with updated element
#######################################################

Sampler_tau_ind_log = function (i, tau.curr, n.proc, data, dist, rho, m, m.denovo, m.maint, d.denovo, d.maint, b, c, sd.tau, lower, upper, DE.NOVO=1, PRINT.FLAG=0)
{
	lower.log = log10 (lower)
	upper.log = log10 (upper)
	
	lambda.curr = c(as.vector (tau.curr), as.vector(rho), as.vector(m))
	lambda.new = lambda.curr
	log.new.tmp = rnorm (1, log10 (tau.curr[i]), sd.tau)
	log.new = ifelse (log.new.tmp<lower.log | log.new.tmp>upper.log, reflect (log.new.tmp, lower.log, upper.log), log.new.tmp)
	lambda.new[i] = 10^(log.new)
	
	if (PRINT.FLAG)
	{
		print ("current lambda (long):")
		print (lambda.curr)
		print ("new lambda (long):")
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
		cat ("log lik using current lambda (long): ", loglik.curr, "\n")
		cat ("log lik using new lambda (long): ", loglik.new, "\n")
		cat ("difference: ", log.ratio, "\n")
		cat ("result: ", result, "\n")
	}
	
	return (result)
}

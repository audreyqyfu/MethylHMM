#######################################################
# Sampler_rho_ind_unif
#
# 2010.09.12: Pass argument dist to loglikHMM.
# 2010.09.10: Discrete model representation.
# 
# Function to update a single rho using a 
# uniform prior.  
# 
# lambda={3 tau, 3 rho, m of S}
#
# input: i: the i-th pi to be updated
#		 rho.curr: current value of rho; vector of length 3
#		 n.proc: number of hidden processes (=3)
#        data: unprocessed double-stranded data; matrix of N by S
#		 dist: pairwise physical distance in bases; vector of S-1
#		 tau: Pr(0->1) in 1-step transition matrix; vector of length 3
#		 m.denovo: de novo rate of Dnmt1
#		 m.maint: rate of maintenance activity of Dnmt1
#		 d.denovo: de novo rate of Dnmt3 on daughter strand
#		 d.maint: rate of maintenance activity of Dnmt3 on daughter strand
#        b: failure of conversion rate
#        c: inappropriate conversion rate
#        r.m: parameter in the beta distn of m
#        g.m: parameter in the beta distn of m
#        sd.rho: sd in the normal proposal distribution for rho_i; scalar
#		 lower: lower bound in the uniform prior of rho; scalar
#		 upper: upper bound in the uniform prior of rho; scalar
#        DE.NOVO: 0 for assumption 1, 1 for assumption 2
#        PRINT.FLAG: 0 for no output, 1 for output
# output: result: rho vector with updated element
#######################################################

Sampler_rho_ind_unif = function (i, rho.curr, n.proc, data, dist, tau, m, m.denovo, m.maint, d.denovo, d.maint, b, c, sd.rho, lower, upper, DE.NOVO=1, PRINT.FLAG=0)
{
	lambda.curr = c(tau, rho.curr, m)
	
	##################################################
	# generate proposal for the i-th rho
	##################################################
	lambda.new = lambda.curr
	par.new = rnorm (1, rho.curr[i], sd.rho)
	lambda.new[n.proc+i] = ifelse ((par.new<lower) | (par.new>upper), reflect (par.new, lower, upper), par.new)
		
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
		result = lambda.new[n.proc+1:n.proc]
	else 
		result = rho.curr

	return (result)
}

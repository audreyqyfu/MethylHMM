#######################################################
# Sampler_rho_ind_log
#
# 2010.09.30: Bug fix: should treat log10 transformed rho
#			  as a new random variable and update this r.v.
#			  for the uniform prior.  Should not transform 
#			  it back to the raw scale for updating!
# 2010.09.12: Pass argument dist to loglikHMM.
#
# Modfified from Sampler_tau_ind_log.
# 
# Function to update a single rho with a log uniform prior.
#
# lambda={3 tau, 3 rho, m of S}
#
# input: i: the i-th rho to be updated
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
#		 lower: lower bound in the log uniform prior of rho; scalar, not log10 transformed
#		 upper: upper bound in the log uniform prior of rho; scalar, not log10 transformed
#        DE.NOVO: 0 for assumption 1, 1 for assumption 2
#        PRINT.FLAG: 0 for no output, 1 for output
# output: result: tau vector with updated element
#######################################################

Sampler_rho_ind_log = function (i, rho.curr, n.proc, data, dist, tau, m, m.denovo, m.maint, d.denovo, d.maint, b, c, sd.rho, lower, upper, DE.NOVO=1, PRINT.FLAG=0)
{
	lower.log = log10 (lower)
	upper.log = log10 (upper)
	
	lambda.curr = c(tau, rho.curr, m)
	lambda.new = lambda.curr
	log.new.tmp = rnorm (1, log10 (rho.curr[i]), sd.rho)
	log.new = ifelse (log.new.tmp<lower.log | log.new.tmp>upper.log, reflect (log.new.tmp, lower.log, upper.log), log.new.tmp)
	lambda.new[n.proc+i] = 10^(log.new)
	
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

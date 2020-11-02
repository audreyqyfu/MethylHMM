#######################################################
# Sampler_lambda
#
# 2011.10.10: Allow for tau_RP=tau_RD and rho_RP=rho_RD.
# 2011.09.24: Add argument M.EST, which may be set to 0 for in vitro data.
# 2010.07.22: allow different priors on the pi.
# 2010.04.16: allow rates of maintanence and de novo activities 
#			  for DNMT3 to be different
# 2009.09.08: modified for a single data file.
# 2009.08.11: call subroutines to separately update pi, rho and m.
# 2009.08.09: parametrization as in Note_10_EmissProbs
# 2008.09.08: under relaxed HMM on maintenance methyltransferases;
#				   see hmm_Aug2108 
# 
# Function to update elements in lambda using both data sets.
# lambda={3 tau, 3 rho, m of S}
#
# input: lambda.curr: current value of lambda; vector of 6+S
#        data: unprocessed double-stranded data; matrix of N by S
#		 dist: pairwise physical distance in bases; vector of length S-1
#		 m.denovo: de novo rate of Dnmt1
#		 m.maint: rate of maintenance activity of Dnmt1
#		 d.denovo: de novo rate of Dnmt3 on daughter strand
#		 d.maint: rate of maintenance activity of Dnmt3 on daughter strand
#        b: failure of conversion rate
#        c: inappropriate conversion rate
#        r.m: parameter in the beta distn of m
#        g.m: parameter in the beta distn of m
#        sd.tau: sd in the normal transition kernel; 
#               vector of length of the number of processes (3 here)
#		 sd.rho: sd in the proposal distribution for rho;
#				 vector of length 3
#        sd.m: s.d. in the normal transition kernel; vector of size 2
#		 tau.lower: vector of length 3; lower bound in the uniform prior for tau (M, RP, RD)
#		 tau.upper: vector of length 3; upper bound in the uniform prior for tau (M, RP, RD)
#		 rho.upper: upper bound in the uniform prior for rho
#		 m.lower: lower bound in the uniform prior for m
#		 m.upper: upper bound in the uniform prior for m
#		 TAU.PRIOR: which prior to assign to tau; uniform, uniform on log10 scale, Jeffreys prior
#		 RHO.PRIOR: which prior to assign to rho; uniform, uniform on log10 scale, Jeffreys prior
#        DE.NOVO: 0 for assumption 1, 1 for assumption 2
#        PRINT.FLAG: 0 for no output, 1 for output
#		 accept.m: acceptance rate of m; vector of size 2
#		 RPRD: 0 for not imposing rho_RP=rho_RD, 1 otherwise
# output: result: updated value of lambda and acceptance count of m
#######################################################

Sampler_lambda = function (lambda.curr, data, dist, m.denovo, m.maint, d.denovo, d.maint, b, c, r.m, g.m, sd.tau, sd.rho, sd.m, tau.lower, tau.upper, rho.lower, rho.upper, m.lower, m.upper, TAU.PRIOR = c("uniform", "logunif", "jeffreys"), RHO.PRIOR=c("uniform", "logunif", "Jeffreys"), DE.NOVO=1, M.EST=1, PRINT.FLAG=0, accept.m, RPRD=0)
{
	n.proc = length (sd.tau)
	n.par = length (lambda.curr)

	tau.curr = as.vector (lambda.curr[1:n.proc])
	rho.curr = as.vector (lambda.curr[(n.proc+1):(2*n.proc)])
	m.curr = lambda.curr[(2*n.proc+1):n.par]

	#######################################
	# update tau
	#######################################
	tau.new = Sampler_tau (tau.curr, data, dist, rho.curr, m.curr, m.denovo, m.maint, d.denovo, d.maint, b, c, sd.tau, tau.lower, tau.upper, TAU.PRIOR, DE.NOVO, PRINT.FLAG, RPRD)
	
	################################################
	# update rho
	################################################
	rho.new = Sampler_rho (rho.curr, data, dist, tau.new, m.curr, m.denovo, m.maint, d.denovo, d.maint, b, c, sd.rho, rho.lower, rho.upper, RHO.PRIOR, DE.NOVO, PRINT.FLAG, RPRD)
	
	########################################
	# update m
	########################################
	if (M.EST) {
		m.new = Sampler_m (m.curr, data, tau.new, rho.new, dist, m.denovo, m.maint, d.denovo, d.maint, b, c, r.m, g.m, sd.m, m.lower, m.upper, accept.m, DE.NOVO, PRINT.FLAG)
	}
	else {
		m.new = list (m=m.curr, accptm=0)
	}

	#########################################
	# combine tau, rho and m into lambda
	#########################################
	lambda.new = c(tau.new, rho.new, m.new$m)

	result=list(lambda=lambda.new, am=m.new$accptm)
	return (result)
	
}

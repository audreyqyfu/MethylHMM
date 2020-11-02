#######################################################
# Sampler_rho
#
# 2011.10.10: Allow for estimating rho_RP=rho_RD.
# 2010.09.30: Allow for Jeffreys prior.
# 2010.09.10: Discrete model representation.
# 
# updating 3 rho using both data sets
#
# lambda={3 tau, 3 rho, m of S}
#
# input: rho.curr: current value of rho; vector of length 3
#        data: unprocessed double-stranded data; matrix of N by S
#		 dist: pairwise physical distance in bases; vector of length S-1
#		 tau: stationary probabilities; vector of length 3
#		 m.denovo: de novo rate of Dnmt1
#		 m.maint: rate of maintenance activity of Dnmt1
#		 d.denovo: de novo rate of Dnmt3 on daughter strand
#		 d.maint: rate of maintenance activity of Dnmt3 on daughter strand
#        b: failure of conversion rate
#        c: inappropriate conversion rate
#        r.m: parameter in the beta distn of m
#        g.m: parameter in the beta distn of m
#        sd.rho: sd in the normal transition kernel; 
#               vector of length of the number of processes (3 here)
#		 rho.upper: upper bound in the uniform prior of rho; scalar
#		 rho.prior.rate: rate parameter in the exponential prior for rho; scalar
#		 RHO.PRIOR: which prior to assign to rho; uniform, uniform on log10 scale, Jeffreys prior
#        DE.NOVO: 0 for assumption 1, 1 for assumption 2
#        PRINT.FLAG: 0 for no output, 1 for output
#		 RPRD: 0 for not imposing rho_RP=rho_RD, 1 otherwise
# output: result: updated rho vector
#######################################################

Sampler_rho = function (rho.curr, data, dist, tau, m, m.denovo, m.maint, d.denovo, d.maint, b, c, sd.rho, lower, upper, RHO.PRIOR=c("uniform", "logunif", "jeffreys"), DE.NOVO=1, PRINT.FLAG=0, RPRD=0)
{
	RHO.PRIOR = match.arg (RHO.PRIOR)
	n.proc = length (rho.curr)
	
	if (RHO.PRIOR=="uniform")	# uniform prior
	{
		result = Sampler_rho_ind_unif (i=1, rho.curr, n.proc, data, dist, tau, m, m.denovo, m.maint, d.denovo, d.maint, b, c, sd.rho[1], lower[1], upper[1], DE.NOVO, PRINT.FLAG)
		rho.curr = result
		if (RPRD) {
			result = Sampler_rho_RPRD_ind_unif (rho.curr, n.proc, data, dist, tau, m, m.denovo, m.maint, d.denovo, d.maint, b, c, sd.rho[2], lower[2], upper[2], DE.NOVO, PRINT.FLAG)
			rho.curr = result
		}
		else {
			for (i in 2:n.proc)
			{
				result = Sampler_rho_ind_unif (i, rho.curr, n.proc, data, dist, tau, m, m.denovo, m.maint, d.denovo, d.maint, b, c, sd.rho[i], lower[i], upper[i], DE.NOVO, PRINT.FLAG)
				rho.curr = result
			}
		}
	}
	else if (RHO.PRIOR=="logunif")	# log uniform prior
	{
		for (i in 1:n.proc)
		{
			result = Sampler_rho_ind_log (i, rho.curr, n.proc, data, dist, tau, m, m.denovo, m.maint, d.denovo, d.maint, b, c, sd.rho[i], lower[i], upper[i], DE.NOVO, PRINT.FLAG)
			rho.curr = result
		}
	}
	else if (RHO.PRIOR=="jeffreys")
	{
		for (i in 1:n.proc)
		{
			result = Sampler_rho_ind_jeffreys (i, rho.curr, n.proc, data, dist, tau, m, m.denovo, m.maint, d.denovo, d.maint, b, c, sd.rho[i], DE.NOVO, PRINT.FLAG)
			rho.curr = result
		}		
	}
	
	return (result)
}

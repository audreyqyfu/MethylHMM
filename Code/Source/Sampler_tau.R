#######################################################
# Sampler_tau
#
# 2011.10.10: Allow for estimating tau_RP=tau_RD.
# 2011.09.24: Fix: log10 uniform prior is now applied to 
#			  all three taus.
# 2010.09.12: Pass argument dist on to loglikHMM.
#
# Update each of 3 tau consecutively.
#
# lambda={3 tau, 3 rho, m of S}
#
# input: tau.curr: current value of tau; vector of length 3
#        data: unprocessed double-stranded data; matrix of N by S
#		 dist: pairwise physical distance in bases; vector of S-1
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
#               vector of length of the number of processes (3 here)
#		 lower: vector of length 3; lower bounds of tau (M, RP, RD)
#		 upper: vector of length 3; upper bounds of tau (M, RP, RD)
#		 TAU.PRIOR: which prior to assign to tau; uniform, uniform on log10 scale, Jeffreys prior
#        DE.NOVO: 0 for assumption 1, 1 for assumption 2
#        PRINT.FLAG: 0 for no output, 1 for output
#		 RPRD: 0 for not imposing rho_RP=rho_RD, 1 otherwise
# output: result: updated tau vector
#######################################################

Sampler_tau = function (tau.curr, data, dist, rho, m, m.denovo, m.maint, d.denovo, d.maint, b, c, sd.tau, lower, upper, TAU.PRIOR = c("uniform", "logunif", "jeffreys"), DE.NOVO=1, PRINT.FLAG=0, RPRD=0)
{
	n.proc = length (tau.curr)

	#######################################
	# update tau one at a time
	#######################################
	if (TAU.PRIOR=="uniform")
	{
		result = Sampler_tau_ind (i=1, tau.curr, n.proc, data, dist, rho, m, m.denovo, m.maint, d.denovo, d.maint, b, c, sd.tau[1], lower[1], upper[1], DE.NOVO, PRINT.FLAG)
		tau.curr = result
		if (RPRD) {
			result = Sampler_tau_RPRD_ind (tau.curr, n.proc, data, dist, rho, m, m.denovo, m.maint, d.denovo, d.maint, b, c, sd.tau[2], lower[2], upper[2], DE.NOVO, PRINT.FLAG)
			tau.curr = result
		}
		else {
			for (i in 2:n.proc)
			{
				result = Sampler_tau_ind (i, tau.curr, n.proc, data, dist, rho, m, m.denovo, m.maint, d.denovo, d.maint, b, c, sd.tau[i], lower[i], upper[i], DE.NOVO, PRINT.FLAG)
				tau.curr = result
			}
		}
	}
	else if (TAU.PRIOR=="logunif")
	{
		for (i in 1:n.proc)
		{
#			if (i==1)
#				result = Sampler_tau_ind (i, tau.curr, n.proc, data, dist, rho, m, m.denovo, m.maint, d.denovo, d.maint, b, c, sd.tau[i], lower[i], upper[i], DE.NOVO, PRINT.FLAG)
#			else
				result = Sampler_tau_ind_log (i, tau.curr, n.proc, data, dist, rho, m, m.denovo, m.maint, d.denovo, d.maint, b, c, sd.tau[i], lower[i], upper[i], DE.NOVO, PRINT.FLAG)
			
			tau.curr = result
		}
	}
	else if (TAU.PRIOR=="jeffreys")
	{
		for (i in 1:n.proc)
		{
			result = Sampler_tau_ind_jeffreys (i, tau.curr, n.proc, data, dist, rho, m, m.denovo, m.maint, d.denovo, d.maint, b, c, sd.tau[i], DE.NOVO, PRINT.FLAG)
			tau.curr = result
		}		
	}
	
	return (result)
}

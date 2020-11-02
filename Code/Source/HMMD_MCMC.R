###########################################################################
# HMMD_MCMC
# 
# 2011.10.10: Add RPRD to allow for tau_RP=tau_RD and rho_RP=rho_RD.
# 2011.09.24: Add argument M.EST, which may be set to 0 for in vitro data.
# 2011.08.30: Typo in an message: when DNMT3.EST=0, shold print 
#			  "setting both d.maint and d.denovo to 1"
# 2010.09.12: Dist, replacing loc, is now an argument.
# 2010.09.12: Bug: argument dist is needed in all functions that call loglikHMM.
# 2010.09.09: MAJOR CHANGE: from a continuous-time model to a discrete-time model.
#			  Specify 1-step transition matrix and multiply it d times to 
#			  get transition matrix over d nt.
#			  Estimate in 1-step transition matrix tau=Pr(0->1) 
#			  and rho=Pr(a bound enzyme jumping off DNA) 
#			  such that rho(1-pi) = Pr(1->0).
# 2010.08.07: printing MCMC iteration numbers becomes an option.
# 2010.07.27: estimating mmaint and mdenovo becomes an option.
# 2010.07.21: log transform pi_RP and pi_RD; assign uniform prior.
# 2010.04.16: allow rates of maintanence and de novo activities 
#			  for DNMT3 to be different
#			  To not estimate these two rates, set d.denovo=0 and d.maint=1
# 2009.12.06: initial values for m as argument.
# 2009.09.08: modified for a single data file.
#			  fixed only Sampler_c.
# 2009.08.09: parametrization as in Note_10_EmissProbs
#			  need to finish sampler for gamma
# 2008.09.08: compute and output strand type probs
# 2008.08.19: under the relaxed model on maintenance methyltransferases;
#                  emission probs in hmm_July2808
#
# function to sample lambda in the hidden Markov model directly from 
# the data
# 
# as a discrete Markov process
# b fixed; estimate c, constant across sites
#
# input: file.mcmc: output of mcmc updates of parameters
#		 file.strandtype: output of strand type probs
#        data: unprocessed double-stranded binary data; matrix of N by S
#		 dist: distances between adjacent CpG sites; vector of length S-1
#        tau.curr: initial values of tau; vector of length 3
#		 rho.curr: jump rates; vector of length 3
#		 m.init: initial values for m; vector of length the number of sites
#		 rm.curr: initial value of parameter r in the beta prior of m
#		 gm.curr: initial value of parameter g in the beta prior of m
#		 c.curr: inappropriate conversion error rate
#		 m.denovo.curr: de novo rate of Dnmt1
#		 m.maint.curr: rate of maintenance activity of Dnmt1
#		 d.denovo.curr: de novo rate of Dnmt3 on daughter strand
#		 d.maint.curr: rate of maintenance activity of Dnmt3 on daughter strand
#        n.iter: number of iterations
#        burn.in: percentage of total length of mcmc run
#        step.size: number of iterations between recording mcmc samples
#        sd.tau: s.d. in the normal transition kernel; vector of length 
#               of the number of processes (3 here)
#		 sd.rho: standard deviation in the proposal distribution for rho;
#				 vector of length 3
#		 sd.m: standard deviation in the proposal distribution for m; scalar
#		 sd.rm: standard deviation in the proposal distribution for r.m; scalar
#		 sd.gm: standard deviation in the proposal distribution for log10(g.m); scalar
#		 sd.c: standard deviation in the proposal distribution for c; scalar
#		 sd.m.denovo: sd in the proposal distribution of m.denovo; scalar
#		 sd.m.maint:sd in the proposal distribution of m.maint; scalar
#		 sd.d.denovo: sd in the proposal distribution of d.denovo; scalar
#		 sd.d.maint:sd in the proposal distribution of d.maint; scalar
#		 rho.precision: precision in the estimates of rho; used to compute upper bound
#						in the uniform prior for rho
#		 rho.prior.rate: rate in the exponential prior for rho; scalar
#		 c.ub: upper bound in the uniform prior for c
#        c.lb: lower bound in the uniform prior for c
#		 tau.lower: vector of length 3; lower bounds of tau (M, RP, RD)
#		 tau.upper: vector of length 3; upper bounds of tau (M, RP, RD)
#        seed.value:scalar
#		 TAU.PRIOR: which prior to assign to tau; uniform, uniform on log10 scale, Jeffreys prior
#		 RHO.PRIOR: which prior to assign to rho: uniform, uniform on log10 scale 
#		 DNMT1.EST: 0 (default) for not estimating m.denovo (set to 0) and m.maint (set to 1); 1 for estimating them.
#		 DNMT3.EST: 0 (default) for not estimating d.denovo (set to 0) and d.maint (set to 1); 1 for estimating them.
#		 ERROR.c: 0 for not estimating it and value of c is set by c.curr; 1 (default) for estimating it
#        DE.NOVO: should always has value 1
#		 RPRD: 0 for not imposing rho_RP=rho_RD, 1 otherwise
#        PRINT.FLAG: 0 (default) for no output, 1 for output
#		 PRINT.ITER: 1 (default) for printing MCMC iteration numbers onto screen; 0 otherwise
# output: in file.mcmc: n.iter by 6+n.site+7
#	  columns are: 3 tau, 3 rho, m of size (n.site.1+n.site.2), r.m, g.m, c, m.denovo, m.maint, d.denovo, d.maint
#		  in file.strandtype: n.iter by s
############################################################################

HMMD_MCMC = function (file.mcmc, file.strandtype, data, dist, tau.curr, rho.curr, m.init, rm.curr, gm.curr, c.curr, m.denovo.curr, m.maint.curr, d.denovo.curr, d.maint.curr, error.b, n.iter, burn.in, step.size, sd.tau, sd.rho, sd.m, sd.rm, sd.gm, sd.c, sd.m.denovo, sd.m.maint, sd.d.denovo, sd.d.maint, tau.lower, tau.upper, rho.lower, rho.upper, c.ub, c.lb, seed.value, TAU.PRIOR=c("uniform", "logunif", "jeffreys"), RHO.PRIOR=c("uniform", "logunif", "jeffreys"), DNMT1.EST=0, DNMT3.EST=0, M.EST=1, ERROR.c=1, DE.NOVO=1, RPRD=0, PRINT.FLAG=0, PRINT.ITER=1)
{
    set.seed (seed.value)
	
	########################################
	# check conditions
	########################################
	TAU.PRIOR = match.arg (TAU.PRIOR)
	RHO.PRIOR = match.arg (RHO.PRIOR)
	if (DNMT1.EST)
	{
		cat ("estimating m.denovo or m.maint\n")
	}
	else
	{
		cat ("not estimating m.denovo or m.maint\n")
		cat ("setting m.denovo to be 0 and m.maint to be 1; input values are ignored\n")
		m.denovo.curr=0
		m.maint.curr=1
	}

	if (DNMT3.EST)
	{
		cat ("estimating d.denovo or d.maint\n")
#		if (DNMT1.EST)
#		{
#			if (TAU.PRIOR=="jeffreys")
#				stop ("only uniform or log10 uniform prior can be assigned to tau\n")
#			if ((tau.upper[1] - tau.lower[1] + tau.upper[3] - tau.lower[3]) > 1)
#				stop ("prior support for tau.M cannot overlap with prior support for tau.RD") 
#		}
	}
	else
	{
		cat ("not estimating d.denovo or d.maint\n")
		cat ("setting d.denovo to be 1 and d.maint to be 1; input values are ignored\n")
		d.denovo.curr=1
		d.maint.curr=1
		cat ("priors on tau are chosen to be:", TAU.PRIOR, "\n")
	}
	
	cat ("priors on rho are chosen to be:", RHO.PRIOR, "\n")

	if (M.EST) {
		cat ("estimating site-specific methylation probabilities m_j\n")
	}
	else {
		cat ("not estimating site-specific methylation probabilities m_j\n")
		cat ("m_j's are set to be:", m.init, "\n")
		cat ("r.m and g.m are set to 0 as well\n")
		rm.curr = 0
		gm.curr = 0
	}
	
	if (RPRD) {
		cat ("imposing tau_RP=tau_RD and rho_RP=rho_RD\n")
		cat ("using initial values for rho_RP and tau_RP, and sd values for RP\n")
		cat ("ignoring values set for the RD process\n")
	}
	
	if (ERROR.c)
	{
		cat ("estimating inapproripate conversion error rate c\n")
	}
	else
	{
		cat ("not estimating inappropriate conversion error rate c\n")
		cat ("c is set to be:", c.curr, "\n")
	}
	
	########################################
	# data processing
	########################################
	data = as.matrix (data)
	n.site = ncol(data)			# number of sites
	n = nrow(data)/2			# number of patterns

	########################################
	# initialization
	########################################
	n.proc = 3
	accept.tau = rep (0, n.proc)			# Pr(0->1) in 1-step transition matrix; (M,RP,RD)
	accept.rho = rep (0, n.proc)		# Pr(a bound enzyme jumping off DNA) in 1-step transition matrix; (M,RP,RD)
	accept.mpar = rep(0,2)			# r.m, g.m
	accept.c = 0						# c
	accept.m = 0
	accept.m.denovo = 0
	accept.m.maint = 0
	accept.d.denovo = 0
	accept.d.maint = 0

	need.append = 0

	# convert to alpha.m and beta.m
	alpha.m = rm.curr*(1-gm.curr)/gm.curr
	beta.m = (1-rm.curr)*(1-gm.curr)/gm.curr
	# initial values for m
	m.curr = m.init

	lambda.curr = c (tau.curr, rho.curr, m.curr)
	
	#####################################################
	# MCMC iterations
	#####################################################
	for (k in 1:n.iter)
	{
		if (PRINT.ITER)	cat (k)

		#####################################
		# compute strand type probs
		#####################################
		strandtype = StrandAssignmentProbs (data, lambda.curr, dist, m.denovo.curr, m.maint.curr, d.denovo.curr, d.maint.curr, error.b, c.curr, DE.NOVO, PRINT.FLAG)
	  
		#####################################
		# sample lambda
		#####################################
		result.lambda = Sampler_lambda (lambda.curr, data, dist, m.denovo.curr, m.maint.curr, d.denovo.curr, d.maint.curr, error.b, c.curr, rm.curr, gm.curr, sd.tau, sd.rho, sd.m, tau.lower, tau.upper, rho.lower, rho.upper, m.lower=0, m.upper=1, TAU.PRIOR, RHO.PRIOR, DE.NOVO, M.EST, PRINT.FLAG, accept.m, RPRD)
		lambda.new = result.lambda$lambda
		tau.new = as.vector (lambda.new[1:n.proc])
		rho.new = as.vector (lambda.new[(n.proc+1):(2*n.proc)])
		m.new = lambda.new[(2*n.proc+1):length(lambda.new)]
		lambda.new = c(tau.new, rho.new, m.new)
		accept.m = result.lambda$am

		######################################
		# sample r.m and g.m
		######################################
		if (M.EST) {
			rm.new = Sampler_r (rm.curr, m.new, gm.curr, sd.rm, 0, 1)
			gm.new = Sampler_g (gm.curr, m.new, rm.new, sd.gm, -4, 0)
		}
		else {
			rm.new = rm.curr
			gm.new = gm.curr
		}

		#####################################
		# sample c
		#####################################
		if (ERROR.c)
			c.new = Sampler_c (c.curr, data, dist, lambda.new, m.denovo.curr, m.maint.curr, d.denovo.curr, d.maint.curr, error.b, sd.c, c.lb, c.ub, DE.NOVO, PRINT.FLAG)
		else
			c.new=c.curr

		######################################
		# sample m.denovo and m.maint
		######################################
		if (DNMT1.EST)
		{
			m.denovo.new = Sampler_mdenovo (m.denovo.curr, data, dist, lambda.new, m.maint.curr, d.denovo.curr, d.maint.curr, error.b, c.new, sd.m.denovo, lower=0, upper=1, DE.NOVO, PRINT.FLAG)
			m.maint.new = Sampler_mmaint (m.maint.curr, data, dist, lambda.new, m.denovo.new, d.denovo.curr, d.maint.curr, error.b, c.new, sd.m.maint, lower=0, upper=1, DE.NOVO, PRINT.FLAG)
		}
		else
		{
			m.denovo.new = m.denovo.curr
			m.maint.new = m.maint.curr
		}
		
		######################################
		# sample d.denovo and d.maint
		######################################
		if (DNMT3.EST)
		{
			d.denovo.new = Sampler_ddenovo (d.denovo.curr, data, dist, lambda.new, d.maint.curr, m.denovo.new, m.maint.new, error.b, c.new, sd.d.denovo, lower=0, upper=1, DE.NOVO, PRINT.FLAG)
			d.maint.new = Sampler_dmaint (d.maint.curr, data, dist, lambda.new, d.denovo.new, m.denovo.new, m.maint.new, error.b, c.new, sd.d.maint, lower=0, upper=1, DE.NOVO, PRINT.FLAG)
		}
		else
		{
			d.denovo.new = d.denovo.curr
			d.maint.new = d.maint.curr
		}
		
		######################################
		# acceptance rate of MCMC samples
		######################################
		for (i in 1:n.proc)
			if (lambda.new[i]!=lambda.curr[i])
				accept.tau[i] = accept.tau[i]+1
		for (i in 1:n.proc)
			if (lambda.new[n.proc+i]!=lambda.curr[n.proc+i])
				accept.rho[i] = accept.rho[i]+1
		if (ERROR.c & c.new!=c.curr)
			accept.c = accept.c+1	
		if (DNMT1.EST)
		{
			if (m.denovo.new != m.denovo.curr)
				accept.m.denovo = accept.m.denovo + 1
			if (m.maint.new != m.maint.curr)
				accept.m.maint = accept.m.maint + 1
		}
		if (DNMT3.EST)
		{
			if (d.denovo.new != d.denovo.curr)
				accept.d.denovo = accept.d.denovo + 1
			if (d.maint.new != d.maint.curr)
				accept.d.maint = accept.d.maint + 1
		}
		if (rm.new!=rm.curr)
			accept.mpar[1] = accept.mpar[1]+1
		if (gm.new!=gm.curr)
			accept.mpar[2] = accept.mpar[2]+1

		######################################
		# output MCMC samples
		######################################
		result = c (lambda.new[1:(2*n.proc)], m.new, rm.new, gm.new, c.new, m.denovo.new, m.maint.new, d.denovo.new, d.maint.new)
		ncol.output = length (result)
		if ((k > (n.iter * burn.in)) && ((k %% step.size) == 0))
		{
			if (need.append == 0)
			{
				write (result, file.mcmc, ncol = ncol.output, sep = "\t")
				write (strandtype, file.strandtype, ncol=length (strandtype), sep="\t")
				need.append <- 1
			}
			else
			{
				write (result, file.mcmc, ncol = ncol.output, append = TRUE, sep = "\t")
				write (strandtype, file.strandtype, ncol=length (strandtype), append=TRUE, sep="\t")
			}
		}

		######################################
		# update current parameters
		######################################
		lambda.curr = lambda.new
		rm.curr = rm.new
		gm.curr = gm.new
		c.curr = c.new
		m.curr = m.new
		m.denovo.curr = m.denovo.new
		m.maint.curr = m.maint.new
		d.denovo.curr = d.denovo.new
		d.maint.curr = d.maint.new
    }
	if (PRINT.ITER)	cat ("\n")

    #################################
    # output acceptance rates
    #################################
	if (ERROR.c)
		cat ("acceptance rate of c is:", accept.c/n.iter, "\n")
    cat ("acceptance rate of tau is:", accept.tau/n.iter, "\n")
	cat ("aceeptance rate of rho is:", accept.rho/n.iter, "\n")
	if (M.EST) {
		cat ("acceptance rate of m is:", accept.m/n.iter, "\n")
		cat ("acceptance rates of r.m and g.m are:", accept.mpar/n.iter, "\n")
	}
	if (DNMT1.EST)
	{
		cat ("acceptance rate of m.denovo is:", accept.m.denovo/n.iter, "\n")
		cat ("acceptance rate of m.maint is:", accept.m.maint/n.iter, "\n")
	}
	if (DNMT3.EST)
	{
		cat ("acceptance rate of d.denovo is:", accept.d.denovo/n.iter, "\n")
		cat ("acceptance rate of d.maint is:", accept.d.maint/n.iter, "\n")
	}
}



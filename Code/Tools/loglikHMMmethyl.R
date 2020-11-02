#######################################################################
# loglikHMMmethyl
#
# 2011.09.11: Incorporate dm and dd.
# 2010.10.11: Modified for the discrete HMM.
# 
# Compute log likelihood for all ds sequences under the HMM
# with x hidden states and observed states of Q',D'.
# If x=8, hidden states are M,RP,RD
# If x=16, hidden states are P, M, RP, RD.
# If x=64, hidden states are Q,D,P,M,RP,RD.
#
# Input: ds.data: matrix of 2N by S
#		 dist: vector of length S-1
#		 pi.mcmc: vector of length 3
#		 rho.mcmc: vector of length 3
#		 m.mcmc: vector of length S
#		 mm.mcmc: scalar
#		 md.mcmc: scalar
#		 dm.mcmc: scalar
#		 dd.mcmc: scalar
#		 c.mcmc: scalar
#		 b: scalar
# Output: vector of length N
#######################################################################
loglikHMMmethyl = function (ds.data, dist, tau.mcmc, rho.mcmc, m.mcmc, mm.mcmc=1, md.mcmc=0, dm.mcmc=1, dd.mcmc=1, c.mcmc, b, model=c("h8", "h16", "h64"))
{
	model = match.arg (model)
		
	# convert (Q',D') to value in {1, 2, 3, 4}
	N = nrow(ds.data)/2 #number of individuals
	S = ncol (ds.data)
	top = as.matrix(ds.data[2*(0:(N-1))+1,]) # odd rows
	bot = as.matrix(ds.data[2*(0:(N-1))+2,]) # even rows
	data.t = 2*top+bot+1
	data.b = 2*bot+top+1
	
	loglik = rep (0, N)

	# compute stationary probability pi using tau and rho
	pi.mcmc = tau.mcmc / (tau.mcmc + rho.mcmc * (1-tau.mcmc))
	
	# compute marginal distribution at first site
	latentMarginal = computeLatentMarginal (pi=pi.mcmc, m=m.mcmc[1], mm=mm.mcmc, md=md.mcmc, dm=dm.mcmc, dd=dd.mcmc, model=model)
	# compute emission matrix
	trans = computeLatentTrans (tau=tau.mcmc, rho=rho.mcmc, dist=dist, m=m.mcmc[-1], mm=mm.mcmc, md=md.mcmc, dm=dm.mcmc, dd=dd.mcmc, model=model)
	# compute transition matrix
	emiss.error = computeLatentEmissError (m=m.mcmc, mm=mm.mcmc, md=md.mcmc, dm=dm.mcmc, dd=dd.mcmc, error.b=b, error.c=c.mcmc, model=model)
	
	for (i in 1:N)
	{
		loglik[i] = forward_backward (data.t[i,], latentMarginal, trans, emiss.error)
		if (sum (top[i,]-bot[i,]==0) < S)
		{
			tmp = forward_backward (data.b[i,], latentMarginal, trans, emiss.error)
			loglik[i] = log (exp(loglik[i]) + exp(tmp))
		}
	}
	
	return (loglik)
	
}

###############################################################
# ViterbiHMMmethyl
#
# 2011.09.11: Incorporate dm and dd.
# 2010.10.11: Modified for discrete HMM. Add comments.
#
# Function to apply the Viterbi algorithm to the ds seq data
# under a hidden Markov model with x hidden states and observed 
# states of Q',D'.
# If x=8, hidden states are M,RP,RD
# If x=16, hidden states are P, M, RP, RD.
# If x=64, hidden states are Q,D,P,M,RP,RD.
# 
# Input: ds.data: matrix of 2N by S; rows (1,2), (3,4), etc are 
#				  top and bottom strands of a pattern
#		 dist: vector of length S-1; distances in bp between adjacent sites
#		 tau.mcmc: vector of 3; estimated associating probability for M, RP, RD
#		 rho.mcmc: vector of 3; estimated dissociating probability for M, RP, RD
#		 m.mcmc: vector of length S; estimated methylation probability at S sites
#		 mm.mcmc: scalar; estimated maintenance activity rate of DNMT1
#		 md.mcmc: scalar; estimated de novo activity rate of DNMT1
#		 dm.mcmc: scalar; estimated maintenance activity rate of DNMT3
#		 dd.mcmc: scalar; estimated de novo activity rate of DNMT3
#		 c.mcmc: scalar; estimated inappropriate bisulfite conversion error rate
#		 b: scalar; failured bisulfite conversion error rate
#		 loglik: vector of length N; each log Pr(pattern i)
#		 model: choose from h8, h16 and h64; specifies number of hidden states
# Output: list of 4 objects:
#		  condprob: vector of length N; each Pr(Viterbi path|a pattern)
#		  path: matrix of N by S; each row the Viterbi path for a pattern
#		  logprobmax: vector of length N; each log Pr(Viterbi path of pattern i, and pattern i)
#		  loglik: vector of length N; each log Pr(pattern i)
################################################################
ViterbiHMMmethyl = function (ds.data, dist, tau.mcmc, rho.mcmc, m.mcmc, mm.mcmc=1, md.mcmc=0, dm.mcmc=1, dd.mcmc=1, c.mcmc, b=0.003, loglik, model=c("h8", "h16", "h64"))
{
	model = match.arg (model)
	
	# convert (Q',D') to value in {1, 2, 3, 4}
	N = nrow(ds.data)/2 #number of individuals
	S = ncol (ds.data)
	top = as.matrix(ds.data[2*(0:(N-1))+1,]) # odd rows
	bot = as.matrix(ds.data[2*(0:(N-1))+2,]) # even rows
	data.t = 2*top+bot+1
	data.b = 2*bot+top+1
	
	logprobmax = rep (0, N)
	path = matrix (0, nrow=N, ncol=S)
	
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
		cat (i)
		result.t = Viterbi (data.t[i,], latentMarginal, emiss.error, trans)
		result.b = Viterbi (data.b[i,], latentMarginal, emiss.error, trans)
		
		if (result.t$logmaxprob > result.b$logmaxprob)
		{
			logprobmax[i] = result.t$logmaxprob
			path[i,] = result.t$path
		}
		else
		{
			logprobmax[i] = result.b$logmaxprob
			path[i,] = result.b$path
		}
	}
	cat ("\n")
	
	# compute Pr (Viterbi path | Data) for each pattern
	condprob = exp (logprobmax - loglik)
	return (list (condprob = condprob, path = path, logprobmax = logprobmax, loglik = loglik))
}

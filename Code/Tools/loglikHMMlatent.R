########################################################################
# loglikHMMlatent
# 
# Computes log of probability Pr(obs, latent) for an observed string.
#
# Input: obs: vector of length T
#		 pi: vector of length nHidden
#		 emiss: array of T by nHidden by nObs
#		 trans: array of T by nHidden by nHidden
# Output:
########################################################################
loglikHMMlatent = function (obs, latent, pi, emiss, trans)
{
	T = length (obs)
	if (is.matrix(emiss))
	{
		emiss.tmp = emiss
		emiss = array (0, dim=c (T, nrow(emiss), ncol(emiss)))
		for (t in 1:T)
		emiss[t,,] = emiss.tmp
	}
	
	logL = log (pi[latent[1]]) + log (emiss[1,latent[1], obs[1]])
	for (t in 2:T)
	logL = logL + log (trans[t-1,latent[t-1], latent[t]]) + log (emiss[t,latent[t], obs[t]])
	
	return (logL)
}


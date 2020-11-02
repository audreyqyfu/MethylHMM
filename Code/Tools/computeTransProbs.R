############################################################################
# computeTransProbs
#
# Computes Pr(0->1) and Pr(1->0) over a set of distances.
#
# Input: tau: Pr(0->1) in 1-step transition matrix; vector of length 3
#        rho: Pr(single bound enzyme becomes unassociated) in 
#			  1-step transition matrix; vector of length 3
#		 dist: distances; integer vector of length S
# Output: list of two components: 
#		  probs0to1: Pr(0->1) for dist; matrix of S by 3
#		  probs1to0: Pr(1->0) for dist; matrix of S by 3
#############################################################################
computeTransProbs = function (tau, rho, dist, n.proc)
{
	S = length (dist)
	probs0to1 = matrix (0, nrow=S, ncol=n.proc)
	probs1to0 = probs0to1
	for (j in 1:n.proc)
	{
		elements = cbind (rep (tau[j], S), rep(rho[j], S), dist)
		probs = apply (elements, 1, computeTransProbsInd)
		probs0to1[,j] = probs[1,]
		probs1to0[,j] = probs[2,]
	}

	return (list (probs0to1=probs0to1, probs1to0=probs1to0))
}
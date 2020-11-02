############################################################################
# computeTransProbsInd
#
# Computes Pr(0->1) and Pr(1->0) over one distance.
#
# Input: x: vector of length 3; containing 
#		    tau: Pr(0->1) in 1-step transition matrix; scalar
#           rho: Pr(single bound enzyme becomes unassociated) in 
#			     1-step transition matrix; scalar
#		    dist: distances; integer
# Output: vector of length 2:
#		  probs0to1: Pr(0->1) for dist
#		  probs1to0: Pr(1->0) for dist
#############################################################################
computeTransProbsInd = function (x)
{
	tpm = matrix (c(1-x[1], x[1], x[2]*(1-x[1]), 1-x[2]*(1-x[1])), byrow=TRUE, ncol=2)
	tpm.d = mtx.exp (tpm,x[3])
	result = c(tpm.d[1,2], tpm.d[2,1])
	return (result)
}
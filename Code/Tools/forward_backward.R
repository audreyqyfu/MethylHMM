#######################################################################
# forward_backward
#
# Computes the log likelihood of an observed vector.
#
# Input: obs: vector of length T
#		 pi: vector of length nHidden
#		 trands: array of T by nHidden by nHidden
#		 emiss: array of T by nHidden by nObs
# Output: scalar
#######################################################################
forward_backward = function (obs, pi, trans, emiss, PRINT=FALSE, FORWARD.ONLY=TRUE)
{
	T = length (obs)
	nHidden = length (pi)
	Alpha = matrix (0, nrow=T, ncol=nHidden)
	Beta = Alpha
	
	# forward pass
	Alpha[1,] = pi * as.vector (emiss[1,,obs[1]])
	for (t in 2:T)
	{
		for (i in 1:nHidden)
		{
			Alpha[t,i] = sum (Alpha[t-1,] * trans[t-1,,i]) * emiss[t,i,obs[t]]
		}
	}
	
	# backward pass
	if (!FORWARD.ONLY)
	{
		Beta[T,] = 1
		for (t in (T-1):1)
			for (i in 1:nHidden)
				Beta[t,i] = sum (trans[t,i,] * emiss[t+1,,obs[t+1]] * Beta[t+1,])
	}
	
	lik = sum (Alpha[T,])
	if (PRINT & !FORWARD.ONLY)
	{
		for (t in 1:T)
		{
			lik.each = sum (Alpha[t,] * Beta[t,])
			cat ("site ", t, ":", lik.each, "\n")
		}
	}
	
	return (log (lik))
	
}


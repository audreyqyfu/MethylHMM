#############################################################
# Viterbi
#
# Viterbi algorithm for a hidden Markov model
#
# Input: obs: vector of length T
#		 hidMarginal: vector of length nHidden
#		 emiss: array of T by nHidden by nObs
#		 trans: array of T by nHidden by nHidden
# Output: list containing two objects:
#		  logmaxprob: joint prob of Viterbi path and data on log scale
#		  path: vector
#############################################################
Viterbi = function (obs, hidMarginal, emiss, trans)
{
	T = length (obs)
	N = length (hidMarginal)
	
	prob.matrix = matrix (0, nrow=N, ncol=T)
	path.matrix = prob.matrix
	
	###################################
	# Initialization
	###################################
	prob.matrix[,1] = hidMarginal * emiss[1,, obs[1]]
	
	###################################
	# Recursion
	###################################
	for (t in 2:T)
	{
		for (i in 1:N)
		{
			tmp = prob.matrix[,t-1] * trans[t-1, ,i]
			prob.matrix[i,t] = emiss[t,i, obs[t]] * max (tmp)
			if (max(tmp)==0)
			path.matrix[i,t] = 0
			else
			path.matrix[i,t] = which (tmp==max(tmp))
		}
	}
	
	###################################
	# Termination
	###################################
	prob.max = max (prob.matrix[,T])
	path.find = rep (0, T)
	path.find[T] = which (prob.matrix[,T] == prob.max)
	for (t in (T-1):1)
	path.find[t] = path.matrix[path.find[t+1], t+1]
	
	return (list(logmaxprob = log(prob.max), path = path.find))
}


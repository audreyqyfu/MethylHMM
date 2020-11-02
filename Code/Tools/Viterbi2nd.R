#############################################################
# Viterbi2nd
#
# Finds candidates of the second most likely path for a hidden Markov model.
#
# Input: obs: vector of length T
#		 pi: vector of length nHidden; initial marginal distribution
#			 of hidden states
#		 emiss: array of T by nHidden by nObs
#		 trans: array of T by nHidden by nHidden
# Output: path candidates: matrix of T by T; i-th row is the Viterbi path
#		  among a set of paths that have the first i-1 elements identical to 
#		  the Viterbi path of data.
#############################################################
Viterbi2nd = function (obs, pi, emiss, trans, vpath)
{
	T = length (obs)
	N = length (pi)
	
	prob.max=rep(0, T)
	path.find = matrix (0, nrow=T, ncol=T)
	
	prob.matrix = matrix (0, nrow=N, ncol=T)
	path.matrix = prob.matrix
	
	prob.matrix[,1] = pi * emiss[1,, obs[1]]
	
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

	for (t in 1:(T-1))
	{
		Tnew = T-t+1
		prob.matrix.tmp = prob.matrix[,t:T]
		path.matrix.tmp = matrix (0, nrow=N, ncol=Tnew)
		prob.matrix.tmp[vpath[t], 1] = 0
		
		if (t<T)
		{
			for (k in 2:Tnew)
			{
				for (i in 1:N)
				{
					tmp = prob.matrix.tmp[,k-1] * trans[t+k-2,,i]
					prob.matrix.tmp[i,k] = emiss[t+k-1, i, obs[t+k-1]] * max (tmp)
					if (max(tmp)==0)
						path.matrix.tmp[i,k] = 0
					else
						path.matrix.tmp[i,k] = which (tmp==max(tmp))
				}
			}
		}
		
		
		path.find.tmp = rep (0, Tnew)
		prob.max.tmp = max (prob.matrix.tmp[,Tnew])
		
		
		path.find.tmp[Tnew] = which (prob.matrix.tmp[,Tnew] == prob.max.tmp)
		for (j in (Tnew-1):1)
		path.find.tmp[j] = path.matrix.tmp[path.find.tmp[j+1], j+1]
		path.find[t,1:(t-1)] = vpath[1:(t-1)]
		path.find[t,t:T] = path.find.tmp
	}
	
	path.find[T,1:(T-1)] = vpath[1:(T-1)]
	prob.matrix.tmp = prob.matrix[,T]
	prob.matrix.tmp[vpath[T]] = 0
	path.find[T,T] = which (prob.matrix.tmp == max (prob.matrix.tmp))
	
	# compute log likelihood of each path
#	loglikpath = rep (0, T)
#	for (t in 1:T)
#		loglikpath[t] = loglikHMMlatent (obs, path.find[t,], pi, emiss, trans)
	
#	logprob = max (loglikpath)
#	path = path.find[which (loglikpath == logprob),]
	
#	return (list (logprob=logprob, path = path, path.all = path.find))
	return (path.find)
}

#############################################################
# Viterbi2ndPathProb
#
# Finds the second most likely path and computes its probability Pr(path|data)
# for a hidden Markov model.
#
# Input: obs: vector of length T
#		 pi: vector of length nHidden; initial marginal distribution
#			 of hidden states
#		 emiss: array of T by nHidden by nObs
#		 trans: array of T by nHidden by nHidden
#		 vpath: vector of length T; Viterbi path of obs
#		 loglik: scalar; Pr(obs)
# Output: list of two objects:
#		  path: matrix, each row a 2nd most likely path
#		  condprob: scalar, Pr(path|obs)
#############################################################
Viterbi2ndPathProb = function (obs, pi, emiss, trans, vpath, loglik)
{
	T = length (obs)
	# generate path candidates in matrix of T by T
	candidates = Viterbi2nd (obs, pi, emiss, trans, vpath)
	# compute log of joint probability of path and obs for each candidate path
	# generate a vector of length T
	logprob = 0
	for (t in 1:T)
		logprob[t] = loglikHMMlatent (obs, candidates[t,], pi, emiss, trans)
	# find max log joint probability of path
	maxlogprob = max (logprob)
	# find indices corresponding to paths with max log joint probability
	maxind = which (logprob==maxlogprob)
	
	# find 2nd most likely paths
	# could be multiple
	# generate matrix with each row being a path
	path = candidates[maxind,]
	# compute Pr(path|obs)
	condprob = exp (maxlogprob - loglik)
	
	return (list (path=path, condprob=condprob))
}
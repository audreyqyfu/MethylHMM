##########################################################################
# loglikHMM 
#
# 2010.09.12: Bug: forgot to have dist as input argument.
# 2010.09.10: Discrete model representation.
# 2010.04.16: allow rates of maintanence and de novo activities 
#			  for DNMT3 to be different
# 2009.08.09: parametrization as in Note_10_EmissProbs
# 2009.07.26: file name and function name changed.
#
# Wrapper function that calls the C subroutine "loglikHMM" 
# to compute log likelihood of all double-stranded sequences
# under the hidden Markov model.
#
# input: data
#        lambda: tau, vector of length 3; rho, vector of length 3; m of length n.site
#		 m.denovo:  de novo rate of Dnmt1
#		 m.maint: rate of maintenance activity of Dnmt1
#		 d.denovo: de novo rate of Dnmt3 on daughter strand
#		 d.maint: rate of maintenance activity of Dnmt3 on daughter strand
#        error.b: failure of conversion
#        error.c: inappropriate conversion
# output: log likelihood given the data and current parameter values
##########################################################################

loglikHMM = function (data, dist, lambda, m.denovo, m.maint, d.denovo, d.maint, error.b, error.c, DE.NOVO=1, PRINT.FLAG=0)
{
	N=as.integer(nrow(data)/2)
	S=as.integer(ncol(data))
	n.hidden=as.integer(8)
	n.obs=as.integer(4)
	n.proc=as.integer(3)
	data.tmp=as.integer(as.vector(t(data)))
	
	tau = as.vector (lambda[1:n.proc])
	rho = as.vector (lambda[(n.proc+1):(2*n.proc)])
	m = as.vector (lambda[(2*n.proc+1):length(lambda)])

	# compute transition probabilities over nucleotide distances
	probs = computeTransProbs (tau, rho, dist, n.proc)
	prob.0to1 = probs$probs0to1
	prob.1to0 = probs$probs1to0
	# compute stationary probabilities pi=Pr(in state 1)
	pi = tau / (tau + rho *(1-tau))
	# compile parameter vector lambda
	lambda.tmp = c(pi, as.vector (prob.0to1), as.vector (prob.1to0), m)
	
	# assign right data types
	lambda.tmp=as.double(lambda.tmp)
	m.denovo.tmp = as.double (m.denovo)
	m.maint.tmp = as.double (m.maint)
	d.denovo.tmp = as.double (d.denovo)
	d.maint.tmp = as.double (d.maint)
	error.b.tmp=as.double(error.b)
	error.c.tmp=as.double(error.c)

	out=.C("loglikHMM",loglik=as.double(0), data.tmp, lambda.tmp, m.denovo.tmp, m.maint.tmp, d.denovo.tmp, d.maint.tmp, error.b.tmp, error.c.tmp, N, S, n.hidden, n.obs, n.proc, as.integer(DE.NOVO), as.integer(PRINT.FLAG))
	return(out$loglik)
}


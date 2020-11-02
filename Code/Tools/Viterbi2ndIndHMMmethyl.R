#############################################################
# Viterbi2ndIndHMMmethyl
#
# 2011.09.11: Incorporate dm and dd.
# 2010.10.11: Modified for discrete HMM.
#
# Finds the second most likely path for a single sequence.
#
# Input: ds.data: matrix of 2 by S; top and bottom strands of a methylation pattern
#		 vpath: vector of length S; Viterbi path of the pattern
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
#		 loglik: scalar; log Pr(pattern i)
#		 stprob: scalar; Pr(top parent and bottom daughter)
#		 model: choose from h8, h16 and h64; specifies number of hidden states
# Output: list containing log Pr(O,Q) and Viterbi path
#############################################################
Viterbi2ndIndHMMmethyl = function (ds.data, vpath, dist, tau.mcmc, rho.mcmc, m.mcmc, mm.mcmc=1, md.mcmc=0, dm.mcmc=1, dd.mcmc=1, c.mcmc, b=0.003, loglik, stprob, model=c("h8", "h16", "h64"))
{
	# use strand assignment probability to determine which strand type to go with
	parent="t"
	# convert (Q',D') to value in {1, 2, 3, 4}
	data = 2*ds.data[1,]+ds.data[2,]+1
	if (stprob<0.5)
	{
		parent="b"
		# convert (Q',D') to value in {1, 2, 3, 4}
		data = 2*ds.data[2,]+ds.data[1,]+1
	}
	
	# compute stationary probability pi using tau and rho
	pi.mcmc = tau.mcmc / (tau.mcmc + rho.mcmc * (1-tau.mcmc))
	
	# compute marginal distribution at first site
	latentMarginal = computeLatentMarginal (pi=pi.mcmc, m=m.mcmc[1], mm=mm.mcmc, md=md.mcmc, dm=dm.mcmc, dd=dd.mcmc, model=model)
	# compute emission matrix
	trans = computeLatentTrans (tau=tau.mcmc, rho=rho.mcmc, dist=dist, m=m.mcmc[-1], mm=mm.mcmc, md=md.mcmc, dm=dm.mcmc, dd=dd.mcmc, model=model)
	# compute transition matrix
	emiss.error = computeLatentEmissError (m=m.mcmc, mm=mm.mcmc, md=md.mcmc, dm=dm.mcmc, dd=dd.mcmc, error.b=b, error.c=c.mcmc, model=model)

	# find 2nd most likely paths for each strand assignment
	v2nd = Viterbi2ndPathProb (data, latentMarginal, emiss.error, trans, vpath, loglik)
#	if (is.matrix(v2nd$path))
#	{
#		v2nd.convert = list (0)
#		for	(i in nrow(v2nd$path))
#			v2nd.convert[[i]] = summaryOnePath (ds.data, v2nd$path[i,])
#	}
#	else
#		v2nd.convert = summaryOnePath (ds.data, v2nd$path)

#	return (list (parent=parent, stprob=stprob, condprob=v2nd$condprob, path=v2nd$path, summary=v2nd.convert))
	return (list (parent=parent, stprob=stprob, condprob=v2nd$condprob, path=v2nd$path))
}
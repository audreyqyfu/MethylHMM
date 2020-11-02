###############################################################
# computePostProbEventHemi
# 
# Compute Pr(methylation event | a hemi) for all hemis in data.
# Four types of methylation events:
# failure of maintenance
# de novo on parent
# de novo on daughter
# measurement error
# Represented by Pr (P,Q,D|a hemi)
#
# Input: ds.data: matrix of 2N by S
#		 assign.pr: matrix of mcmc.length by N (number of patterns)
#		 m: matrix of mcmc.length by S (number of sites)
#	     c: vector of mcmc.length
#		 four parameters below are from the HMM
#		 tau: matrix of mcmc.length by 3
#		 rho: matrix of mcmc.length by 3
#		 mm: vector of mcmc.length
#		 md: vector of mcmc.length
#		 three parameters below are from the IEM
#		 mu: matrix of mcmc.length by S
#		 dp: matrix of mcmc.length by S
#		 dd: matrix of mcmc.length by S
#		 b: scalar; failed bisulfite conversion error rate
# Output: array of mcmc.length by nHemi by 4 (events)
###############################################################
computePostProbEventHemi = function (ds.data, assign.pr, m, c, tau, rho, md, mm, mu, dp, dd, b=0.003, model=c("HMM", "IEM"))
{
	model=match.arg(model)
	mcmc.length = nrow(assign.pr)
	S = ncol (ds.data)
	
	########################################
	# find coordinates of hemis
	########################################
	h.coord = findHemiCoord (ds.data)
	n.hemi = nrow (h.coord)
	
	if (model=="HMM")
	{
		########################################
		# compute stationary prob pi
		########################################
		pi = tau / (tau + rho*(1-tau))
	}
	
	post = array (0, dim = c(mcmc.length, n.hemi, 4))
	for (k in 1:mcmc.length)
	{
		print (k)
		
		# Pr (top strand being parent | a pattern)
		prk = as.vector (assign.pr[k,])
		h.prk = prk[h.coord[,1]]
		# Pr((Q',D')=(1,0)|hemi)
		post.tp = h.prk * h.coord[,3] + (1-h.prk) * (1-h.coord[,3])
		
		# compute Pr((Q',D')|(Q,D))
		error.QD = array (0, dim=c(S,4,4))
		if (is.vector (c))
		{
			error = matrix (c(1-b,b,c[k],1-c[k]), byrow=TRUE, ncol=2)
			error.QD.tmp = kronecker (error, error)
			for (s in 1:S)
				error.QD[s,,] = error.QD.tmp
		}
		else
		{
			for (s in 1:S)
			{
				error = matrix (c(1-b,b,c[k,s],1-c[k,s]), byrow=TRUE, ncol=2)
				error.QD[s,,] = kronecker (error, error)
			}
		}
		
		emiss.PQD = array (0, dim=c(S, 2, 4))
		if (model=="HMM")
		{
			#################################################
			# compute Pr((Q,D)|P) by summing over (M, RP, RD)
			#################################################
			# compute emission probs Pr(Q,D|P,M,RP,RD)
			# produce 16-by-4 matrix
			emiss.h16 = computeLatentEmiss (mm=mm[k], md=md[k], model="h16")
			# compute marginals Pr(M,RP,RD)
			marg.h8 = computeLatentMarginal (pi=pi[k,], model="h8")
			marg.h16.tmp = c(marg.h8, marg.h8)
			# compute Pr(P,M,RP,RD, Q,D) 
			# matrix of 16 by 4; each entry a joint probability
			joint.h64 = emiss.h16 * marg.h16.tmp
			# compute Pr(Q,D|P)
			# matrix of 2 by 4; each entry a conditional probability
			emiss.PQD.tmp = rbind (apply (joint.h64[1:8,], 2, sum), apply (joint.h64[9:16,], 2, sum))
			
			for (s in 1:S)
				emiss.PQD[s,,] = emiss.PQD.tmp
		}
		if (model=="IEM")
		{
			for (s in 1:S)
			{
				emiss.PQD[s,1,] = c((1-dp[k,s])*(1-dd[k,s]), (1-dp[k,s])*dd[k,s], dp[k,s]*(1-dd[k,s]), dp[k,s]*dd[k,s])
				emiss.PQD[s,2,3:4] = c(1-mu[k,s], mu[k,s])
			}
		}
		
		#################################################
		# compute Pr(P,Q,D | Q',D')
		# array of S by 4 by 8
		#################################################
		post.PQD = array (0, dim=c(S, 4, 8))
		for (s in 1:S)
		{
			# Pr (P, Q,D)
			joint.PQD = emiss.PQD[s,,] * c(1-m[k,s], m[k,s])
			# Pr (P,Q,D, Q',D')
			joint.PQD.error = rbind (error.QD[s,,] * joint.PQD[1,], error.QD[s,,]*joint.PQD[2,])
			# Pr(Q',D')
			colsums = apply (joint.PQD.error, 2, sum)
			# Pr(P,Q,D | Q',D')
			post.PQD[s,,] = t(joint.PQD.error) / colsums
		}
		
#		cat ("finish computing post PQD\n")
#		print (dim (post.PQD))
		
		#################################################
		# compute Pr(event| hemi)
		# deriving from Pr (P,Q,D | hemi)
		# array of mcmc.length by nHemi by 4
		#################################################
		for (t in 1:n.hemi)
		{
#			print (t)
#			print (h.coord[t,1])
#			print (post.PQD[h.coord[t,2],2,])
#			print (post.tp[t])
			# Pr(P,Q,D, (Q',D')=(0,1)|H)
			post.01 = post.PQD[h.coord[t,2],2,] * (1-post.tp[t])
#			print (post.01)
			# Pr(P,Q,D, (Q',D')=(1,0)|H)
			post.10 = post.PQD[h.coord[t,2],3,] * post.tp[t]
#			print (post.10)
			# failure of maintenance event
			post[k,t,1] = post.01[7] + post.10[7]
#			print (post[k,t,1])
			# de novo on parent 
			post[k,t,2] = post.01[3] + post.10[3]
#			print (post[k,t,2])
			# de novo on daughter
			post[k,t,3] = post.01[2] + post.10[2]
#			print (post[k,t,3])
			# measurement error
			post[k,t,4] = 1-sum (post[k,t,1:3])
#			print (post[k,t,4])
		}
#		cat ("finish computing post event\n")
	}
	return (post)
}

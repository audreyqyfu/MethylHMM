##############################################################################
# computeLatentTrans
#
# 2011.09.11: Incorporate dm and dd.
# 2010.10.11: Modify for the discrete HMM.
#
# Computes the transition probability matrix for the hidden Markov chain.
#
# Input: tau: vector of length 3
#		 rho: vector of length 3
#		 dist: distances between sites; vector of length S
#		 m: methylation densities at sites excluding the first one; vector of length S
# Output: array of S by nHidden by nHidden
##############################################################################
computeLatentTrans = function (tau, rho, dist, m, mm, md, dm, dd, n.proc=3, n.obs=4, model=c("h8", "h16", "h64"))
{
	model = match.arg (model)
	S = length (dist)
	
	probs = computeTransProbs (tau, rho, dist, n.proc)
	# matrix of S by 3 (M, RP, RD)
	prob.0to1 = probs$probs0to1
	# matrix of S by 3 (M, RP, RD)
	prob.1to0 = probs$probs1to0

	if (model == "h8")
	{
		n.state = 8
		trans = array (0, dim=c(S,n.state, n.state))
		for (s in 1:S)
		{
			single = array (0, dim=c(n.proc, 2, 2))
			for (k in 1:n.proc)
			{
				single[k,1,1] = 1 - prob.0to1[s,k]
				single[k,2,1] = prob.1to0[s,k]
				single[k,,2] = 1 - single[k,,1]
				
				if (k==1)
					joint = single[k,,]
				else
					joint = kronecker (joint, single[k,,])
			}
			
			trans[s,,] = joint
		}		
	}
	
	if (model == "h16")
	{
		n.state = 16
		trans = array (0, dim=c(S,n.state, n.state))
		for (s in 1:S)
		{
			m.single = matrix (0, nrow=2, ncol=2)
			m.single[,1] = 1-m[s]
			m.single[,2] = m[s]
		
			single = array (0, dim=c(n.proc, 2, 2))
			joint = m.single
			for (k in 1:n.proc)
			{
				single[k,1,1] = 1 - prob.0to1[s,k]
				single[k,2,1] = prob.1to0[s,k]
				single[k,,2] = 1 - single[k,,1]
			
				joint = kronecker (joint, single[k,,])
			}
		
			trans[s,,] = joint
		}
	}

	if (model == "h64")
	{
		# Pr(Q2,D2,P2,M2,RP2,RD2|Q1,D1,P1,M1,RP1,RD1)
		# = Pr(Q2,D2|P2,M2,RP2,RD2)Pr(P2,M2,RP2,RD2|Q1,D1,P1,M1,RP1,RD1)
		# = Pr(Q2,D2|P2,M2,RP2,RD2)Pr(P2,M2,RP2,RD2|P1,M1,RP1,RD1)
		# = emiss for h16 * trans for h16
		n.state = 64
		n.state.tmp = 16
		trans = array (0, dim=c(S,n.state, n.state))
		for (s in 1:S)
		{
			m.single = matrix (0, nrow=2, ncol=2)
			m.single[,1] = 1-m[s]
			m.single[,2] = m[s]
			
			single = array (0, dim=c(n.proc, 2, 2))
			joint = m.single
			for (k in 1:n.proc)
			{
				single[k,1,1] = 1 - prob.0to1[s,k]
				single[k,2,1] = prob.1to0[s,k]
				single[k,,2] = 1 - single[k,,1]
				
				joint = kronecker (joint, single[k,,])
			}
			
			emiss = computeLatentEmiss (mm=mm, md=md, dm=dm, dd=dd, n.obs=n.obs, model="h16")	
			
			for (i in 1:n.obs)
			{
				for (j in 1:n.obs)
				{
					trans[s,n.state.tmp*(i-1)+1:n.state.tmp,n.state.tmp*(j-1)+1:n.state.tmp] = t(t(joint) * emiss[,j])
				}
			}
		}
	}
		
	
	return (trans)
}


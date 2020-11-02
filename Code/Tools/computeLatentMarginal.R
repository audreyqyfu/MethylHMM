###############################################################
# computeLatentMarginal
# 
# 2011.09.11: Incorporate dm and dd.
# 2010.10.11: Simplify code and add comments.
# 
# Computes initial marginal distribution for the HMM with x hidden states.
# If x=8, hidden states are M,RP,RD
# If x=16, hidden states are P, M, RP, RD.
# If x=64, hidden states are Q,D,P,M,RP,RD.
# 
# Input: pi: vector of length 3
#		 m: scalar; methylation density at site 1
#		 mm: scalar
#		 md: scalar
# Output: vector of length nHidden
###############################################################
computeLatentMarginal = function (pi, m, mm=1, md=0, dm=1, dd=1, model=c("h8", "h16", "h64"))
{
	model = match.arg (model)
	marginal = 0
	
	if (model == "h8")
	{
		# Pr (M,RP,RD)
		marginal[1] = (1-pi[1]) * (1-pi[2]) * (1-pi[3])
		marginal[2] = (1-pi[1]) * (1-pi[2]) * (pi[3])
		marginal[3] = (1-pi[1]) * (pi[2]) * (1-pi[3])
		marginal[4] = (1-pi[1]) * (pi[2]) * (pi[3])
		
		marginal[5] = (pi[1]) * (1-pi[2]) * (1-pi[3])
		marginal[6] = (pi[1]) * (1-pi[2]) * (pi[3])
		marginal[7] = (pi[1]) * (pi[2]) * (1-pi[3])
		marginal[8] = (pi[1]) * (pi[2]) * (pi[3])
	}
	else
	{
		# h16
		# Pr (P,M,RP,RD)
		marginal[1] = (1-pi[1]) * (1-pi[2]) * (1-pi[3])
		marginal[2] = (1-pi[1]) * (1-pi[2]) * (pi[3])
		marginal[3] = (1-pi[1]) * (pi[2]) * (1-pi[3])
		marginal[4] = (1-pi[1]) * (pi[2]) * (pi[3])
		
		marginal[5] = (pi[1]) * (1-pi[2]) * (1-pi[3])
		marginal[6] = (pi[1]) * (1-pi[2]) * (pi[3])
		marginal[7] = (pi[1]) * (pi[2]) * (1-pi[3])
		marginal[8] = (pi[1]) * (pi[2]) * (pi[3])
		
		marginal[9:16] = marginal[1:8]
		marginal[1:8] = marginal[1:8] * (1-m)
		marginal[9:16] = marginal[9:16] * m
		
		if (model=="h64")
		{
			# Pr (Q,D,P,M,RP,RD) = Pr(Q,D|P,M,RP,RD)Pr(P,M,RP,RD)
			emiss = computeLatentEmiss (mm=mm, md=md, dm=dm, dd=dd, model="h16")	
			marginalh64 = as.vector (emiss * marginal)
			marginal = marginalh64			
		}
	}
	
	return (marginal)
}
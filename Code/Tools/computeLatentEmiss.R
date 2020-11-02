##############################################################################
# computeLatentEmiss
#
# 2011.09.11: Incorporate dm and dd.
# 2010.10.11: Add comments.
#
# Computes the emission matrix (without error) for the hidden Markov chain.
#
# Input: m: vector of length S; needed only for h8
#		 mm: scalar
#		 md: scalar
#		 dm: scalar
#		 dd: scalar
# Output: if h8, array of S by nHidden by nObs;
#		  if h16, matrix of nHidden by nObs
##############################################################################
computeLatentEmiss = function (m, mm=1, md=0, dm=1, dd=1, n.obs=4, model=c("h8", "h16"))
{
	model = match.arg (model)

	if (model == "h8")
	{
		# Pr ((Q,D)|(M,RP,RD))
		# emission matrix different at different positions
		# produce an array of S by n.state by n.obs
		n.state = 8
		S = length (m)
		emiss = array (0, dim = c(S, n.state, n.obs))
		
		for (s in 1:S)
		{
			emiss[s,1,1] = 1-m[s]
			emiss[s,1,3] = m[s]
			emiss[s,2,1] = (1.0 - m[s]) * (1.0 - dd)
			emiss[s,2,2] = (1-m[s]) * dd
			emiss[s,2,3] = m[s] * (1-dm)
			emiss[s,2,4] = m[s] * dm
			emiss[s,3,3] = 1
			emiss[s,4,3] = (1-m[s])*(1-dd) + m[s]*(1-dm)
			emiss[s,4,4] = (1-m[s])*dd + m[s]*dm
			
			emiss[s,5,1] = (1-m[s])*(1 - md)
			emiss[s,5,2] = (1-m[s])*md
			emiss[s,5,3] = m[s]*(1-mm)
			emiss[s,5,4] = m[s]*mm
			emiss[s,6,1] = (1-m[s]) * (1-md)*(1-dd)
			emiss[s,6,2] = (1-m[s]) * (md+dd-md*dd)
			emiss[s,6,3] = m[s]*(1-mm)*(1-dm)
			emiss[s,6,4] = m[s]*(mm+dm-mm*dm)
			emiss[s,7,3] = (1-m[s])*(1 - md) + m[s]*(1-mm)
			emiss[s,7,4] = (1-m[s])*md + m[s]*mm
			emiss[s,8,3] = m[s]*(1-mm)*(1-dm) + (1-m[s])*(1-md)*(1-dd)
			emiss[s,8,4] = m[s]*(mm+dm-mm*dm) + (1-m[s])*(md+dd-md*dd)
		}
	}
	
	if (model == "h16")
	{
		# Pr ((Q,D)|(P,M,RP,RD))
		# emission matrix is the same at all positions
		# hence calculate only at one position here
		n.state = 16
		emiss = matrix (0, nrow=n.state, ncol=n.obs)
		emiss[1,1] = 1
		emiss[2,1] = 1-dd
		emiss[2,2] = dd
		emiss[3,3] = 1
		emiss[4,3] = 1-dd
		emiss[4,4] = dd
		
		emiss[5,1] = 1 - md
		emiss[5,2] = md
		emiss[6,1] = (1-md)*(1-dd)
		emiss[6,2] = md + dd - md*dd
		emiss[7,3] = 1 - md
		emiss[7,4] = md
		emiss[8,3] = (1-md)*(1-dd)
		emiss[8,4] = md+dd-md*dd
		
		emiss[9,3] = 1
		emiss[10,3] = 1-dm
		emiss[10,4] = dm
		emiss[11,3] = 1
		emiss[12,3] = 1-dm
		emiss[12,4] = dm
		
		emiss[13,3] = 1 - mm
		emiss[13,4] = mm
		emiss[14,3] = (1-mm)*(1-dm)
		emiss[14,4] = mm+dm-mm*dm
		emiss[15,3] = 1 - mm
		emiss[15,4] = mm
		emiss[16,3] = (1-mm)*(1-dm)
		emiss[16,4] = mm+dm-mm*dm
	}
	
	return (emiss)
}


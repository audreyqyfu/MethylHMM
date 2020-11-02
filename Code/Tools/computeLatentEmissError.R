###############################################################
# computeLatentEmissError
# 
# 2011.09.11: Incorporate dm and dd.
# 2010.10.11: Simplify calculation under the h64 model.
#			  Add comments.
#
# Computes emission matrix for the HMM with x hidden states.
# If x=8, hidden states are M,RP,RD
# If x=16, hidden states are P, M, RP, RD.
# If x=64, hidden states are Q,D,P,M,RP,RD.
# 
# Input: m: vector of length S
#		 mm: scalar
#		 md: scalar
#		 error.b: scalar
#		 error.c: scalar
# Output: array of S by nHidden by nObs
###############################################################
computeLatentEmissError = function (m, mm=1, md=0, dm=1, dd=1, error.b, error.c, n.obs=4, model=c("h8", "h16", "h64"))
{
	model = match.arg (model)
	
	# compute error matrix Pr ((Q', D')|(Q,D))
	# produce 4-by-4 matrix
	# with columns corresponding to (Q',D')
	# and rows corresponding to (Q,D)
	# in following order: (0,0), (0,1), (1,0), (1,1)
	error.single = matrix (0, nrow=2, ncol=2)
	error.single[,1] = c(1-error.b, error.c)
	error.single[,2] = 1-error.single[,1]
	error = kronecker (error.single, error.single)

	# for this model, 
	# emission matrix with error is not the same at all positions
	# produce an array of S by n.state by n.obs
	if (model == "h8")
	{
		# produce an array of S by 8 by 4
		n.state = 8
		emiss = computeLatentEmiss (m=m, mm=mm, md=md, dm=dm, dd=dd, model="h8")

		S = dim (emiss)[1]
		emiss.error = emiss
		for (s in 1:S)
		emiss.error[s,,] = emiss[s,,] %*% error
	}		

	# for two models below
	# emission matrix with error is the same at all positions
	# compute this matrix at a single position
	# produce n.state-by-n.obs matrix
	if (model == "h16")
	{
		# Pr((Q',D')|(P,M,RP,RD)) = (sum over Q,D) Pr((Q',D')|(Q,D)) Pr ((Q,D)|(P,M,RP,RD))
		n.state = 16
		emiss = computeLatentEmiss (m=m, mm=mm, md=md, dm=dm, dd=dd, model="h16")
		
		# produce a matrix of 16 by 4 (n.obs)
		emiss.error = emiss %*% error
	}		
	
	if (model == "h64")
	{
		# Pr((Q',D')|(Q,D,P,M,RP,RD)) = Pr((Q',D')|(Q,D))
		# in emiss.error, a matrix of 64 by 4 (n.obs)
		# each of rows 1-16 is Pr(.|(0,0)), row 1 of error
		# each of rows 17-32 is Pr(.|(0,1)), row 2 of error
		# each of rows 33-48 is Pr(.|(1,0)), row 3 of error
		# each of rows 49-64 is Pr(.|(1,1)), row 4 of error
		# expand matrix error
		# such that each row is repeated consecutively 16 times
		n.state = 64
		emiss.error = apply (error, 2, rep, each=n.state/nrow(error))
#		n.rep = n.state / nrow (error)
#		emiss.error = matrix (0, nrow=n.state, ncol=n.obs)
#		for (i in 1:nrow (error))
#		{
#			tmp = rep (error[i,], n.rep)
#			tmp.mtx = matrix (tmp, ncol=n.obs, byrow=TRUE)
#			emiss.error[n.rep * (i-1) + 1:n.rep,] = tmp.mtx
#		}
	}		
	
	# for the two models above, 
	# produce an array of S by n.state by n.obs
	# matrix at each S is the same
#	if (is.matrix (emiss.error))
	if (model!="h8")
	{
		emiss.tmp = emiss.error
		S = length (m)
		emiss.error = array (0, dim=c(S, nrow (emiss.tmp), ncol(emiss.tmp)))
		for (s in 1:S)
			emiss.error[s,,] = emiss.tmp
	}
	
	return (emiss.error)
}
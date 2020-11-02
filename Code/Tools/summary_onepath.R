#######################################################################
# summaryOnePath
# 
# 2011.09.11: Bug fix: use new function integer.base.b2 to fill in 
#			  extra digits with 0
#
# Convert object from Viterbi algorithms to readable format.
#
# Input: ds.data: one sequence; matrix of 2 by S
#		 vpath: vector of length S
#		 path.object: list produced from calling Viterbi2ndIndHMMmethyl
# Output: list
#######################################################################
summaryOnePath = function (ds.data, path)
{
	n.col=6
	n.site = ncol (ds.data)
	converted = integer.base.b2 (path-1, L=6)
	colnames (converted) = c("Q","D","P","M","RP","RD")
	converted = cbind (t(ds.data), converted)
	colnames (converted)[1:2] = c("x","y")

	return (t(converted))
}
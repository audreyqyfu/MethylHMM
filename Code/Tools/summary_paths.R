#######################################################################
# summaryPaths
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
summaryPaths = function (ds.data, vpath, vpath.prob, path.object)
{
	n.col=6
	n.site = ncol (ds.data)
	paths = rbind (vpath, path.object$path)
	rownames (paths) = 1:nrow(paths)
	converted = matrix (0, nrow=n.site, ncol=nrow(paths)*n.col)
	for (i in 1:nrow(paths))
		converted[,(i-1)*n.col+1:n.col] = integer.base.b2 (paths[i,]-1, L=6)
	colnames (converted) = rep (c("Q","D","P","M","RP","RD"),nrow(paths))
	converted = cbind (t(ds.data), converted)
	colnames (converted)[1:2] = c("x","y")

	return (list (condprob = c(vpath.prob, path.object$condprob[1]), paths = paths, converted = t(converted)))
}
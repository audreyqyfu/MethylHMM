###########################################
# findHemiCoord
#
# Find coordinates of hemimethylated dyads 
# in double-stranded data. 
#
# Input: ds.data: matrix of 2N by S
# Output: matrix of nhemi by 3 with 
#		  col 1: row index (which pattern)
#		  col 2: column index (which site)
#		  col 3: top methylated (1), bottom 
#				 methylated (0)
###########################################
findHemiCoord = function (ds.data)
{
	N = nrow(ds.data)/2 #number of individuals
	S = ncol (ds.data)
	top = as.matrix(ds.data[2*(0:(N-1))+1,]) # odd rows
	bot = as.matrix(ds.data[2*(0:(N-1))+2,]) # even rows
	h = top-bot
	h.top = which (h==1, arr.ind=TRUE)
	h.bot = which (h==(-1), arr.ind=TRUE)
	status = c(rep (1, nrow (h.top)), rep (0, nrow (h.bot)))
	h.coord = as.matrix (rbind (h.top, h.bot))
	h.coord = cbind (h.coord, status)
	rownames (h.coord) = 1:nrow (h.coord)
	
	return (h.coord)
}

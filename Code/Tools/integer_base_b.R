#########################################
# integer.base.b
#
# author: Spencer Graves 
# http://tolstoy.newcastle.edu.au/R/help/03b/3251.html
#
# Converts decimal to binary
#########################################
integer.base.b <- function(x, b=2)
{ 
	xi <- as.integer(x) 
	if(any(is.na(xi) | ((x-xi)!=0))) 
	print(list(ERROR="x not integer", x=x)) 
	N <- length(x) 
	xMax <- max(x)	
	ndigits <- (floor(logb(xMax, base=2))+1) 
	Base.b <- array(NA, dim=c(N, ndigits)) 
	for(i in 1:ndigits){#i <- 1 
		Base.b[, ndigits-i+1] <- (x %% b) 
		x <- (x %/% b) 
	} 
	if(N ==1) Base.b[1, ] else Base.b 
} 


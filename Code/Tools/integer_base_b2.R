#########################################
# integer.base.b2
#
# based on integer.base.b by Spencer Graves from the following URL
# http://tolstoy.newcastle.edu.au/R/help/03b/3251.html
#
# Converts decimal to binary and fill in 
# additional digits with 0
#########################################
integer.base.b2 <- function(x, b=2, L)
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
	if(ndigits < L) {
		result = array (0, dim=c(N, L))
		result[,(L-ndigits+1):L] = Base.b
	}
	else
		result = Base.b
	if(N ==1) result[1, ] else result
} 


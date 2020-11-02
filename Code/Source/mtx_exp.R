##################################################
# mtx.exp
#
# Computes matrix raised to power n.
# Authors: Vincente Canto Cassola and Martin Maechler
# Code taken from package Biodem.
##################################################
mtx.exp = function (X, n) 
{
    if (n != round(n)) {
        n <- round(n)
        warning("rounding exponent `n' to", n)
    }
    phi <- diag(nrow = nrow(X))
    pot <- X
    while (n > 0) {
        if (n%%2) 
		phi <- phi %*% pot
        n <- n%/%2
        pot <- pot %*% pot
    }
    return(phi)
}

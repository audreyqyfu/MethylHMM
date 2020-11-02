#############################################################################
# Sampler_r 
#
# 2009.09.09: modified for a single input data file; lower and upper bounds 
#			  of the support for r are now arguments.
#
# function to update r via random walk Metropolis. 
# using both data sets.
#
# The prior on r is uniform (0, 1)
# 
# input: current: scalar; current value of r
#        rate: list of 2 vectors; each vector of length s; current values of rate
#        g: the current value of the other hyperparameter
# output: next: scalar; next value of r
##############################################################################

Sampler_r = function (current, rate, g, sd, lower, upper)
{
    # generate a new value from a normal distribution
    # with current as the mean
    new.tmp <- rnorm (1, mean=current, sd=sd)
	new = ifelse (new.tmp>1 | new.tmp<0, reflect (new.tmp, 0, 1), new.tmp)

    alpha.new <- new * (1 - g) / g
    beta.new <- (1 - new) * (1 - g) / g
    alpha.curr <- current * (1 - g) / g
    beta.curr <- (1 - current) * (1 - g) / g

    # compute the Hastings ratio on the log scale
    log.ratio <- sum (dbeta (rate, shape1=alpha.new, shape2=beta.new, log = TRUE)) - sum (dbeta (rate, shape1=alpha.curr, shape2=beta.curr, log = TRUE))

    # choose the value that gives a smaller probability
	result = ifelse ((log.ratio != 'NaN') & (log (runif (1, min=0, max=1)) < log.ratio), new, current) 

    return (result)
}


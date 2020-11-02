#############################################################################
# Sampler_g
#
# 2009.09.09: modified to take a single input data file; reflect proposals 
#			  from outside the support to within the support.
#
# function to update the rate via random walk Metropolis. 
# using both data sets
#
# The prior on g is log uniform (lower, upper)
# 
# input: current: integer; current value of g
#        rate: list of 2 vectors; each vector of length s; current values of rate
#        r: the current value of the other hyperparameter
# output: next: integer; next value of g
##############################################################################

Sampler_g <- function (current, rate, r, sd, lower, upper)
{
    # generate a new value from a normal distribution
    # with current as the mean
    log.new.tmp <- rnorm (1, log10 (current), sd)
	log.new = ifelse (log.new.tmp<lower | log.new.tmp>upper, reflect (log.new.tmp, lower, upper), log.new.tmp)

    alpha.new <- r * (1 - 10^(log.new)) / 10^(log.new)
    beta.new <- (1 - r) * (1 - 10^(log.new)) / 10^(log.new)
    alpha.curr <- r * (1 - current) / current
    beta.curr <- (1 - r) * (1 - current) / current

    log.ratio <- sum (dbeta (rate, alpha.new, beta.new, log = TRUE)) - sum (dbeta (rate, alpha.curr, beta.curr, log = TRUE)) 
	result = ifelse ((log.ratio != 'NaN') & (log(runif(1)) < log.ratio), 10^(log.new), current) 

    return (result)
}


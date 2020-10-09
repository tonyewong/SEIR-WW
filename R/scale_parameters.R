##==============================================================================
## scale_parameters.R
##
## The R function to generate a Latin Hypercube sample only draws from the
## uniform marginal distributions on [0,1]. This function treats those draws as
## samples from the CDF of each parameter, and uses the quantile function to
## map them back to the values for each parameter.
##
##  X = 2d array of parameters. #rows = #samples, #columns = #parameters
##  parameter_names = 1d array of parameter names
##  priors = prior distribution object returned from `priors_func`
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

scale_parameters <- function(X, parameter_names, priors) {

  # scale parameters X (from Unif(0,1)) to whatever distributions are
  # specified below. Assumed that each column of X is a different parameter,
  # which are in the same order as `parameter_names`, and each row is a
  # different concomitant set of parameters.

  X_scaled <- X

  for (pp in names(priors)) {
    ip <- match(pp, parameter_names)
    if (priors[[pp]]$distribution == "gaussian") {
      X_scaled[,ip] <- qnorm(p=X[,ip], mean=priors[[pp]]$params[match("mean", priors[[pp]]$parnames)], sd=priors[[pp]]$params[match("sd", priors[[pp]]$parnames)])

    } else if (priors[[pp]]$distribution == "binomial") {
      X_scaled[,ip] <- qbinom(p=X[,ip], size=priors[[pp]]$params[match("size", priors[[pp]]$parnames)], prob=priors[[pp]]$params[match("p", priors[[pp]]$parnames)])

    } else if (priors[[pp]]$distribution == "neg_binomial") {
      X_scaled[,ip] <- qnbinom(p=X[,ip], size=priors[[pp]]$params[match("size", priors[[pp]]$parnames)], prob=priors[[pp]]$params[match("p", priors[[pp]]$parnames)])

    } else if (priors[[pp]]$distribution == "triangle") {
      X_scaled[,ip] <- qtriangle(p=X[,ip], a=priors[[pp]]$params[match("a", priors[[pp]]$parnames)], b=priors[[pp]]$params[match("b", priors[[pp]]$parnames)], c=priors[[pp]]$params[match("c", priors[[pp]]$parnames)])

    } else if (priors[[pp]]$distribution == "disc_uniform") {
      X_scaled[,ip] <- qdunif(p=X[,ip], min=priors[[pp]]$params[match("min", priors[[pp]]$parnames)], max=priors[[pp]]$params[match("max", priors[[pp]]$parnames)])

    } else if (priors[[pp]]$distribution == "cont_uniform") {
      X_scaled[,ip] <- qunif(p=X[,ip], min=priors[[pp]]$params[match("min", priors[[pp]]$parnames)], max=priors[[pp]]$params[match("max", priors[[pp]]$parnames)])

    } else {
      print(paste("ERROR: unrecognized `priors$distribution` value of",priors[[pp]]$distribution))

    }
  }

  return(X_scaled)

}

##==============================================================================
## End
##==============================================================================
